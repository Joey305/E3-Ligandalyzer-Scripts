#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
10_MCS_Matching.py
===============================================================================
Purpose:
    Build atom-level 2D SMILES <-> 3D PDB ligand mappings using RDKit MCS.

This script creates unique recruiter IDs for every ligand instance:
    L00001, L00002, L00003, ...

A ligand instance is defined as:
    Ligase + pdb_id + Ligand + Variant

Inputs:
    Ligase_Table/Ligase_Ligand_SASA_atoms.csv
        Atom-level 3D SASA table.

    Ligase_Table/Ligase_SMILE_Codes.csv
        Validated SMILES table from metadata generation.
        Used as the preferred SMILES source by Ligase + Ligand.

    Ligases/<Ligase>/PDB/<pdb_id>_<Ligand>.pdb
    Ligases/<Ligase>/PDB/<pdb_id>_<Ligand>_<Variant>.pdb
        Active cleaned PDB files.

Outputs:
    Ligase_Table/Ligase_Ligands_Smiles_3DMapped.csv
        Main atom-level SMILES-to-3D mapping table.

    Ligase_Table/Recruiter_SMILES_Map.csv
        Columns: SMILES, RECRUITER_CODE

    Ligase_Table/Ligase_SMILE_Codes.csv
        Columns: SMILES, RECRUITER_CODE
        Existing file is backed up before overwrite.

    Ligase_Table/Ligase_SMILE_Codes_Atoms.csv
        Columns: RECRUITER_CODE, smiles_atom_index, smile_atom

    Ligase_Table/Ligand_Instance_Recruiter_Codes.csv
        Instance mapping table from biological identity to L-code.

    Ligase_Table/MCS_Matching_Audit.csv
        Per-instance MCS audit.

Usage:
    python 10_MCS_Matching.py --dry-run
    python 10_MCS_Matching.py --overwrite

Notes:
    - SMILES atom indices start at zero.
    - Hydrogens are removed before MCS so mapping is heavy-atom based.
    - If the MCS does not cover all SMILES/PDB heavy atoms, the row is still
      written as partial and flagged in MCS_Matching_Audit.csv.
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


# =============================================================================
# Defaults
# =============================================================================

DEFAULT_ATOMS_CSV = Path("Ligase_Table/Ligase_Ligand_SASA_atoms.csv")
DEFAULT_SMILES_CSV = Path("Ligase_Table/Ligase_SMILE_Codes.csv")
DEFAULT_LIGASE_ROOT = Path("Ligases")
DEFAULT_OUTDIR = Path("Ligase_Table")

MAPPED_COLUMNS = [
    "Ligase",
    "pdb_id",
    "Ligand",
    "Variant",
    "RECRUITER_CODE",
    "Chain",
    "atom_id",
    "exact_atom",
    "atom_type",
    "x",
    "y",
    "z",
    "smiles_atom_index",
    "smile_atom",
]

RECRUITER_MAP_COLUMNS = [
    "SMILES",
    "RECRUITER_CODE",
]

SMILE_CODES_COLUMNS = [
    "SMILES",
    "RECRUITER_CODE",
]

SMILE_ATOMS_COLUMNS = [
    "RECRUITER_CODE",
    "smiles_atom_index",
    "smile_atom",
]

INSTANCE_COLUMNS = [
    "RECRUITER_CODE",
    "Ligase",
    "pdb_id",
    "Ligand",
    "Variant",
    "Original_Ligase_Ligand_Code",
    "SMILES",
    "Canonical_SMILES",
    "SMILES_Source",
    "PDB_File",
    "Instance_Key",
]

AUDIT_COLUMNS = [
    "RECRUITER_CODE",
    "Ligase",
    "pdb_id",
    "Ligand",
    "Variant",
    "PDB_File",
    "Status",
    "SMILES_Source",
    "SMILES_Heavy_Atoms",
    "PDB_Heavy_Atoms",
    "MCS_Atoms",
    "Mapped_Rows",
    "Notes",
    "Error",
]


# =============================================================================
# Small helpers
# =============================================================================

def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def clean_str(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and math.isnan(x):
        return ""
    return str(x).strip()


def normalize_variant(x) -> str:
    s = clean_str(x)
    if not s or s.lower() == "nan":
        return "1"
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except Exception:
        pass
    return s


def backup_if_exists(path: Path) -> Optional[Path]:
    if not path.exists():
        return None
    b = path.with_name(f"{path.stem}.backup_{stamp()}{path.suffix}")
    shutil.copy2(path, b)
    return b


def write_csv(path: Path, columns: List[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# =============================================================================
# Input loading
# =============================================================================

def load_atoms(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing SASA atoms file: {path}")

    df = pd.read_csv(path)

    required = {
        "Ligase", "pdb_id", "Ligand", "Variant", "Chain",
        "atom_id", "exact_atom", "atom_type", "x", "y", "z",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")

    df = df.copy()
    df["Ligase"] = df["Ligase"].astype(str).str.strip()
    df["pdb_id"] = df["pdb_id"].astype(str).str.strip().str.upper()
    df["Ligand"] = df["Ligand"].astype(str).str.strip().str.upper()
    df["Variant"] = df["Variant"].apply(normalize_variant)
    df["Chain"] = df["Chain"].astype(str).str.strip()

    # Keep heavy atoms only for mapping.
    if "is_heavy_atom" in df.columns:
        mask = df["is_heavy_atom"].astype(str).str.lower().isin(["true", "1", "yes"])
        df = df[mask].copy()
    else:
        df = df[df["atom_type"].astype(str).str.upper() != "H"].copy()

    df["atom_id"] = pd.to_numeric(df["atom_id"], errors="coerce").astype("Int64")

    return df


def load_smiles(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing SMILES table: {path}")

    df = pd.read_csv(path)

    if "SMILES" not in df.columns:
        raise ValueError(f"{path} must contain SMILES column.")

    df = df.copy()

    if "Ligase" in df.columns:
        df["Ligase"] = df["Ligase"].astype(str).str.strip()
    if "Ligand" in df.columns:
        df["Ligand"] = df["Ligand"].astype(str).str.strip().str.upper()

    if "Canonical_SMILES" not in df.columns:
        df["Canonical_SMILES"] = df["SMILES"]

    # Key by Ligase + Ligand when available.
    if {"Ligase", "Ligand"}.issubset(df.columns):
        df["_SMILES_KEY"] = df["Ligase"] + "|" + df["Ligand"]
    elif "RECRUITER_CODE" in df.columns:
        # Accept old Ligase_Ligand recruiter code format if that is all we have.
        parts = df["RECRUITER_CODE"].astype(str).str.split("_", n=1, expand=True)
        if parts.shape[1] == 2:
            df["Ligase"] = parts[0]
            df["Ligand"] = parts[1].str.upper()
            df["_SMILES_KEY"] = df["Ligase"] + "|" + df["Ligand"]
        else:
            raise ValueError(
                f"{path} lacks Ligase/Ligand columns and RECRUITER_CODE is not Ligase_Ligand format."
            )
    else:
        raise ValueError(f"{path} must contain Ligase + Ligand or RECRUITER_CODE.")

    return df


def build_smiles_lookup(smiles_df: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    lookup = {}
    for _, r in smiles_df.iterrows():
        key = clean_str(r.get("_SMILES_KEY"))
        smi = clean_str(r.get("SMILES"))
        if not key or not smi:
            continue
        lookup[key] = {
            "SMILES": smi,
            "Canonical_SMILES": clean_str(r.get("Canonical_SMILES")) or smi,
        }
    return lookup


# =============================================================================
# PDB ligand extraction
# =============================================================================

def find_pdb_file(ligase_root: Path, ligase: str, pdb_id: str, ligand: str, variant: str) -> Optional[Path]:
    pdb_dir = ligase_root / ligase / "PDB"

    candidates = []
    if variant and variant != "1":
        candidates.append(pdb_dir / f"{pdb_id}_{ligand}_{variant}.pdb")
    candidates.append(pdb_dir / f"{pdb_id}_{ligand}.pdb")
    candidates.append(pdb_dir / f"{pdb_id}_{ligand}_{variant}.pdb")

    for p in candidates:
        if p.exists():
            return p

    # Fallback glob.
    hits = sorted(pdb_dir.glob(f"{pdb_id}_{ligand}*.pdb"))
    return hits[0] if hits else None


def residue_name_from_pdb_line(line: str) -> str:
    return line[17:20].strip().upper()


def atom_serial_from_pdb_line(line: str) -> Optional[int]:
    try:
        return int(line[6:11])
    except Exception:
        return None


def extract_ligand_pdb_block(pdb_path: Path, ligand: str) -> Tuple[str, int]:
    hetatm_lines: List[str] = []
    conect_lines: List[str] = []
    ligand_serials: List[int] = []

    with pdb_path.open(errors="replace") as f:
        all_lines = f.readlines()

    for line in all_lines:
        if not line.startswith("HETATM"):
            continue
        if residue_name_from_pdb_line(line) == ligand.upper():
            hetatm_lines.append(line)
            serial = atom_serial_from_pdb_line(line)
            if serial is not None:
                ligand_serials.append(serial)

    ligand_serial_set = set(ligand_serials)

    for line in all_lines:
        if not line.startswith("CONECT"):
            continue
        nums = []
        for i in range(6, len(line), 5):
            chunk = line[i:i + 5].strip()
            if chunk.isdigit():
                nums.append(int(chunk))
        if nums and nums[0] in ligand_serial_set:
            conect_lines.append(line)

    return "".join(hetatm_lines + conect_lines), len(hetatm_lines)


def mol_from_pdb_block(pdb_block: str) -> Tuple[Optional[Chem.Mol], str]:
    if not pdb_block.strip():
        return None, "Empty ligand PDB block"

    errors = []

    try:
        mol = Chem.MolFromPDBBlock(
            pdb_block,
            sanitize=True,
            removeHs=False,
            proximityBonding=True,
        )
        if mol is not None:
            return mol, ""
    except Exception as e:
        errors.append(f"sanitize=True failed: {type(e).__name__}: {e}")

    try:
        mol = Chem.MolFromPDBBlock(
            pdb_block,
            sanitize=False,
            removeHs=False,
            proximityBonding=True,
        )
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except Exception as e:
                errors.append(f"post-sanitize failed: {type(e).__name__}: {e}")
            return mol, "; ".join(errors)
    except Exception as e:
        errors.append(f"sanitize=False failed: {type(e).__name__}: {e}")

    return None, "; ".join(errors) if errors else "RDKit returned None"


def remove_hs_keep_info(mol: Chem.Mol) -> Chem.Mol:
    try:
        return Chem.RemoveHs(mol, sanitize=False)
    except Exception:
        return mol


# =============================================================================
# RDKit / MCS
# =============================================================================

def smiles_mol_from_smiles(smiles: str) -> Tuple[Optional[Chem.Mol], str]:
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is not None:
            return mol, ""
    except Exception as e:
        err1 = f"sanitize=True failed: {type(e).__name__}: {e}"
    else:
        err1 = "MolFromSmiles returned None"

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is not None:
            return mol, err1 + "; parsed sanitize=False"
    except Exception as e:
        return None, err1 + f"; sanitize=False failed: {type(e).__name__}: {e}"

    return None, err1


def canonicalize_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def generate_smiles_from_pdb_mol(pdb_mol_no_h: Chem.Mol) -> str:
    try:
        return Chem.MolToSmiles(pdb_mol_no_h, canonical=True, isomericSmiles=True)
    except Exception:
        return ""


def run_mcs(
    smiles_mol: Chem.Mol,
    pdb_mol: Chem.Mol,
    timeout: int,
    bond_compare: str = "any",
) -> Tuple[Optional[Tuple[int, ...]], Optional[Tuple[int, ...]], str]:
    notes = []

    # First try direct substructure, both directions.
    try:
        match_pdb = pdb_mol.GetSubstructMatch(smiles_mol)
        if match_pdb:
            match_smiles = tuple(range(smiles_mol.GetNumAtoms()))
            return match_smiles, match_pdb, "DIRECT_SMILES_IN_PDB"
    except Exception as e:
        notes.append(f"direct smiles-in-pdb failed: {type(e).__name__}")

    try:
        match_smiles = smiles_mol.GetSubstructMatch(pdb_mol)
        if match_smiles:
            match_pdb = tuple(range(pdb_mol.GetNumAtoms()))
            return match_smiles, match_pdb, "DIRECT_PDB_IN_SMILES"
    except Exception as e:
        notes.append(f"direct pdb-in-smiles failed: {type(e).__name__}")

    if bond_compare.lower() == "any":
        bond_cmp = rdFMCS.BondCompare.CompareAny
    else:
        bond_cmp = rdFMCS.BondCompare.CompareOrder

    try:
        mcs = rdFMCS.FindMCS(
            [smiles_mol, pdb_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=bond_cmp,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            matchValences=False,
            timeout=timeout,
        )
    except Exception as e:
        return None, None, "; ".join(notes + [f"FindMCS exception: {type(e).__name__}: {e}"])

    if not mcs or not mcs.smartsString:
        return None, None, "; ".join(notes + ["No MCS SMARTS produced"])

    patt = Chem.MolFromSmarts(mcs.smartsString)
    if patt is None:
        return None, None, "; ".join(notes + ["MCS SMARTS could not be parsed"])

    match_smiles = smiles_mol.GetSubstructMatch(patt)
    match_pdb = pdb_mol.GetSubstructMatch(patt)

    if not match_smiles or not match_pdb:
        return None, None, "; ".join(notes + ["MCS pattern did not match both molecules"])

    return match_smiles, match_pdb, f"MCS:{mcs.smartsString}"


# =============================================================================
# Main mapping
# =============================================================================

def build_instances(atoms_df: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["Ligase", "pdb_id", "Ligand", "Variant"]
    instances = (
        atoms_df[group_cols]
        .drop_duplicates()
        .sort_values(group_cols)
        .reset_index(drop=True)
    )
    instances["RECRUITER_CODE"] = [f"LR{i:05d}" for i in range(1, len(instances) + 1)]
    instances["Instance_Key"] = (
        instances["Ligase"] + "|" +
        instances["pdb_id"] + "|" +
        instances["Ligand"] + "|" +
        instances["Variant"].astype(str)
    )
    instances["Original_Ligase_Ligand_Code"] = instances["Ligase"] + "_" + instances["Ligand"]
    return instances


def process_instance(
    instance: Dict[str, object],
    atoms_group: pd.DataFrame,
    smiles_lookup: Dict[str, Dict[str, str]],
    ligase_root: Path,
    timeout: int,
    bond_compare: str,
    prefer_pdb_smiles: bool,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], Optional[Dict[str, object]], Dict[str, object]]:
    ligase = clean_str(instance["Ligase"])
    pdb_id = clean_str(instance["pdb_id"])
    ligand = clean_str(instance["Ligand"]).upper()
    variant = normalize_variant(instance["Variant"])
    code = clean_str(instance["RECRUITER_CODE"])
    instance_key = clean_str(instance["Instance_Key"])
    original_code = clean_str(instance["Original_Ligase_Ligand_Code"])

    audit = {
        "RECRUITER_CODE": code,
        "Ligase": ligase,
        "pdb_id": pdb_id,
        "Ligand": ligand,
        "Variant": variant,
        "PDB_File": "",
        "Status": "FAILED",
        "SMILES_Source": "",
        "SMILES_Heavy_Atoms": "",
        "PDB_Heavy_Atoms": "",
        "MCS_Atoms": "",
        "Mapped_Rows": 0,
        "Notes": "",
        "Error": "",
    }

    mapped_rows: List[Dict[str, object]] = []
    smiles_atom_rows: List[Dict[str, object]] = []
    instance_row: Optional[Dict[str, object]] = None

    try:
        pdb_path = find_pdb_file(ligase_root, ligase, pdb_id, ligand, variant)
        if pdb_path is None:
            raise FileNotFoundError(f"No PDB found for {ligase} {pdb_id} {ligand} variant {variant}")

        audit["PDB_File"] = str(pdb_path)

        pdb_block, hetatm_count = extract_ligand_pdb_block(pdb_path, ligand)
        if hetatm_count == 0:
            raise ValueError(f"No HETATM lines found for ligand residue {ligand} in {pdb_path}")

        pdb_mol, pdb_note = mol_from_pdb_block(pdb_block)
        if pdb_mol is None:
            raise ValueError(f"Could not parse PDB ligand mol: {pdb_note}")

        pdb_mol_no_h = remove_hs_keep_info(pdb_mol)
        pdb_heavy = pdb_mol_no_h.GetNumAtoms()
        audit["PDB_Heavy_Atoms"] = pdb_heavy

        lookup_key = f"{ligase}|{ligand}"
        smiles_source = "Ligase_SMILE_Codes"

        smiles = ""
        canonical_smiles = ""

        if lookup_key in smiles_lookup:
            smiles = smiles_lookup[lookup_key]["SMILES"]
            canonical_smiles = smiles_lookup[lookup_key]["Canonical_SMILES"]

        pdb_generated_smiles = generate_smiles_from_pdb_mol(pdb_mol_no_h)

        if prefer_pdb_smiles or not smiles:
            if pdb_generated_smiles:
                smiles = pdb_generated_smiles
                canonical_smiles = canonicalize_smiles(smiles)
                smiles_source = "PDB_Ligand_Mol"
            elif not smiles:
                raise ValueError("No SMILES available from table or PDB ligand mol")

        smiles_mol, smi_note = smiles_mol_from_smiles(smiles)
        if smiles_mol is None:
            raise ValueError(f"Could not parse SMILES: {smi_note}")

        smiles_heavy = smiles_mol.GetNumAtoms()
        audit["SMILES_Source"] = smiles_source
        audit["SMILES_Heavy_Atoms"] = smiles_heavy

        for atom in smiles_mol.GetAtoms():
            smiles_atom_rows.append({
                "RECRUITER_CODE": code,
                "smiles_atom_index": atom.GetIdx(),
                "smile_atom": atom.GetSymbol(),
            })

        match_smiles, match_pdb, mcs_note = run_mcs(
            smiles_mol=smiles_mol,
            pdb_mol=pdb_mol_no_h,
            timeout=timeout,
            bond_compare=bond_compare,
        )

        if match_smiles is None or match_pdb is None:
            raise ValueError(f"MCS failed: {mcs_note}")

        conf = pdb_mol_no_h.GetConformer()

        # Optional lookup from SASA atoms by atom serial.
        atom_lookup = {}
        for _, r in atoms_group.iterrows():
            try:
                atom_lookup[int(r["atom_id"])] = r
            except Exception:
                pass

        for idx2d, idx3d in zip(match_smiles, match_pdb):
            atom3d = pdb_mol_no_h.GetAtomWithIdx(int(idx3d))
            pos = conf.GetAtomPosition(int(idx3d))
            ri = atom3d.GetPDBResidueInfo()

            atom_serial = ri.GetSerialNumber() if ri is not None else atom3d.GetIdx()
            exact_atom = ri.GetName().strip() if ri is not None else atom3d.GetSymbol()
            chain = ri.GetChainId().strip() if ri is not None else ""

            # Prefer SASA table values where available.
            sasa_row = atom_lookup.get(int(atom_serial))
            if sasa_row is not None:
                chain = clean_str(sasa_row.get("Chain")) or chain
                exact_atom = clean_str(sasa_row.get("exact_atom")) or exact_atom
                atom_type = clean_str(sasa_row.get("atom_type")) or atom3d.GetSymbol()
                x = sasa_row.get("x")
                y = sasa_row.get("y")
                z = sasa_row.get("z")
            else:
                atom_type = atom3d.GetSymbol()
                x, y, z = pos.x, pos.y, pos.z

            mapped_rows.append({
                "Ligase": ligase,
                "pdb_id": pdb_id,
                "Ligand": ligand,
                "Variant": variant,
                "RECRUITER_CODE": code,
                "Chain": chain,
                "atom_id": atom_serial,
                "exact_atom": exact_atom,
                "atom_type": atom_type,
                "x": x,
                "y": y,
                "z": z,
                "smiles_atom_index": int(idx2d),
                "smile_atom": smiles_mol.GetAtomWithIdx(int(idx2d)).GetSymbol(),
            })

        mcs_atoms = len(mapped_rows)
        audit["MCS_Atoms"] = mcs_atoms
        audit["Mapped_Rows"] = mcs_atoms

        notes = [mcs_note]
        if pdb_note:
            notes.append(pdb_note)
        if smi_note:
            notes.append(smi_note)

        if mcs_atoms == smiles_heavy == pdb_heavy:
            status = "PASS_FULL_MCS"
        elif mcs_atoms == min(smiles_heavy, pdb_heavy):
            status = "PASS_SUBSET_MCS"
            notes.append("MCS covers the smaller molecule but not both heavy-atom counts")
        else:
            status = "CHECK_PARTIAL_MCS"
            notes.append("MCS does not cover full heavy-atom set")

        audit["Status"] = status
        audit["Notes"] = "; ".join([n for n in notes if n])

        instance_row = {
            "RECRUITER_CODE": code,
            "Ligase": ligase,
            "pdb_id": pdb_id,
            "Ligand": ligand,
            "Variant": variant,
            "Original_Ligase_Ligand_Code": original_code,
            "SMILES": smiles,
            "Canonical_SMILES": canonical_smiles or canonicalize_smiles(smiles),
            "SMILES_Source": smiles_source,
            "PDB_File": str(pdb_path),
            "Instance_Key": instance_key,
        }

        return mapped_rows, smiles_atom_rows, instance_row, audit

    except Exception as e:
        audit["Error"] = f"{type(e).__name__}: {e}"
        audit["Notes"] = traceback.format_exc()
        return mapped_rows, smiles_atom_rows, instance_row, audit


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build SMILES atom to 3D PDB atom mappings using RDKit MCS."
    )
    parser.add_argument("--atoms", default=str(DEFAULT_ATOMS_CSV), help="Input Ligase_Ligand_SASA_atoms.csv")
    parser.add_argument("--smiles", default=str(DEFAULT_SMILES_CSV), help="Input Ligase_SMILE_Codes.csv")
    parser.add_argument("--ligase-root", default=str(DEFAULT_LIGASE_ROOT), help="Ligases root")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory")
    parser.add_argument("--timeout", type=int, default=20, help="MCS timeout in seconds per ligand")
    parser.add_argument("--bond-compare", choices=["any", "order"], default="any", help="MCS bond comparison mode")
    parser.add_argument("--prefer-pdb-smiles", action="store_true", help="Use SMILES generated from the PDB ligand mol instead of the metadata table when possible")
    parser.add_argument("--overwrite", action="store_true", help="Backup and overwrite output files")
    parser.add_argument("--dry-run", action="store_true", help="Preview counts without writing output files")
    args = parser.parse_args()

    atoms_csv = Path(args.atoms)
    smiles_csv = Path(args.smiles)
    ligase_root = Path(args.ligase_root)
    outdir = Path(args.outdir)

    outputs = {
        "mapped": outdir / "Ligase_Ligands_Smiles_3DMapped.csv",
        "recruiter_map": outdir / "Recruiter_SMILES_Map.csv",
        "smile_codes": outdir / "Ligase_SMILE_Codes.csv",
        "smile_atoms": outdir / "Ligase_SMILE_Codes_Atoms.csv",
        "instances": outdir / "Ligand_Instance_Recruiter_Codes.csv",
        "audit": outdir / "MCS_Matching_Audit.csv",
        "log": outdir / "MCS_Matching.log",
    }

    if not args.dry_run:
        existing = [p for p in outputs.values() if p.exists()]
        if existing and not args.overwrite:
            print("❌ Output files already exist:")
            for p in existing:
                print(f"   - {p}")
            print("Run with --overwrite to back them up and regenerate.")
            return 2

        if args.overwrite:
            for p in existing:
                b = backup_if_exists(p)
                if b:
                    print(f"🧷 Backed up {p} -> {b}")

    print("=== 10 MCS Matching ===")
    print(f"Started: {now()}")
    print(f"Atoms input:  {atoms_csv}")
    print(f"SMILES input: {smiles_csv}")
    print(f"Ligase root:  {ligase_root}")
    print(f"Output dir:   {outdir}")
    print(f"MCS timeout:  {args.timeout}s")
    print(f"Bond compare: {args.bond_compare}")
    print(f"Prefer PDB SMILES: {args.prefer_pdb_smiles}")
    print(f"Dry run: {args.dry_run}")

    atoms_df = load_atoms(atoms_csv)
    smiles_df = load_smiles(smiles_csv)
    smiles_lookup = build_smiles_lookup(smiles_df)

    instances = build_instances(atoms_df)

    print(f"\n📥 Loaded heavy atom rows: {len(atoms_df)}")
    print(f"🧪 Loaded SMILES lookup entries: {len(smiles_lookup)}")
    print(f"🧬 Unique ligand instances to map: {len(instances)}")

    if args.dry_run:
        print("\n🧪 Dry run complete. No files written.")
        print("Run for real with:")
        print("  python 10_MCS_Matching.py --overwrite")
        return 0

    mapped_all: List[Dict[str, object]] = []
    smile_atom_all: List[Dict[str, object]] = []
    instance_rows: List[Dict[str, object]] = []
    audit_rows: List[Dict[str, object]] = []

    grouped_atoms = {
        (
            str(k[0]),
            str(k[1]),
            str(k[2]),
            normalize_variant(k[3]),
        ): g.copy()
        for k, g in atoms_df.groupby(["Ligase", "pdb_id", "Ligand", "Variant"], dropna=False)
    }

    for i, (_, inst) in enumerate(instances.iterrows(), 1):
        key = (
            clean_str(inst["Ligase"]),
            clean_str(inst["pdb_id"]),
            clean_str(inst["Ligand"]).upper(),
            normalize_variant(inst["Variant"]),
        )
        atom_group = grouped_atoms.get(key)
        if atom_group is None:
            atom_group = atoms_df.iloc[0:0].copy()

        mapped, smile_atoms, inst_row, audit = process_instance(
            instance=inst.to_dict(),
            atoms_group=atom_group,
            smiles_lookup=smiles_lookup,
            ligase_root=ligase_root,
            timeout=args.timeout,
            bond_compare=args.bond_compare,
            prefer_pdb_smiles=args.prefer_pdb_smiles,
        )

        mapped_all.extend(mapped)
        smile_atom_all.extend(smile_atoms)
        if inst_row is not None:
            instance_rows.append(inst_row)
        audit_rows.append(audit)

        if i % 50 == 0 or i == len(instances):
            pass_count = sum(1 for r in audit_rows if str(r["Status"]).startswith("PASS"))
            fail_count = sum(1 for r in audit_rows if r["Status"] == "FAILED")
            check_count = len(audit_rows) - pass_count - fail_count
            print(f"📊 Progress {i}/{len(instances)} | PASS={pass_count} CHECK={check_count} FAILED={fail_count}")

    # Recruiter map and SMILES codes are now instance-specific L-codes.
    recruiter_map_rows = [
        {"SMILES": r["SMILES"], "RECRUITER_CODE": r["RECRUITER_CODE"]}
        for r in instance_rows
    ]

    smile_codes_rows = [
        {"SMILES": r["SMILES"], "RECRUITER_CODE": r["RECRUITER_CODE"]}
        for r in instance_rows
    ]

    write_csv(outputs["mapped"], MAPPED_COLUMNS, mapped_all)
    write_csv(outputs["recruiter_map"], RECRUITER_MAP_COLUMNS, recruiter_map_rows)
    write_csv(outputs["smile_codes"], SMILE_CODES_COLUMNS, smile_codes_rows)
    write_csv(outputs["smile_atoms"], SMILE_ATOMS_COLUMNS, smile_atom_all)
    write_csv(outputs["instances"], INSTANCE_COLUMNS, instance_rows)
    write_csv(outputs["audit"], AUDIT_COLUMNS, audit_rows)

    status_counts = pd.Series([r["Status"] for r in audit_rows]).value_counts().to_dict()
    log_payload = {
        "started_finished": now(),
        "atoms_input": str(atoms_csv),
        "smiles_input": str(smiles_csv),
        "ligase_root": str(ligase_root),
        "instances": len(instances),
        "mapped_rows": len(mapped_all),
        "smiles_atom_rows": len(smile_atom_all),
        "instance_rows": len(instance_rows),
        "status_counts": status_counts,
        "outputs": {k: str(v) for k, v in outputs.items()},
    }

    outputs["log"].write_text(json.dumps(log_payload, indent=2) + "\n")

    print("\n=== MCS Matching Summary ===")
    print(f"Ligand instances processed: {len(instances)}")
    print(f"Instance rows written:      {len(instance_rows)}")
    print(f"3D mapped rows written:     {len(mapped_all)}")
    print(f"SMILES atom rows written:   {len(smile_atom_all)}")
    print("Status counts:")
    for k, v in status_counts.items():
        print(f"  {k}: {v}")

    print(f"\n✅ Written {outputs['mapped']}")
    print(f"✅ Written {outputs['recruiter_map']}")
    print(f"✅ Written {outputs['smile_codes']}")
    print(f"✅ Written {outputs['smile_atoms']}")
    print(f"✅ Written {outputs['instances']}")
    print(f"✅ Written {outputs['audit']}")
    print(f"🧾 Written {outputs['log']}")

    if any(r["Status"] == "FAILED" for r in audit_rows):
        print("\n⚠️ Some instances failed. Review:")
        print(f"  {outputs['audit']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
