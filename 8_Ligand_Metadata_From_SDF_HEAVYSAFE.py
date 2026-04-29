#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
8_Ligand_Metadata_From_SDF.py
===============================================================================
Author: Joseph-Michael Schulz workflow patch

Purpose:
    Build ligand/recruiter metadata from the ACTUAL local SDF structures, while
    keeping only SDFs that still correspond to ligands present in the current
    Ligand_SASA_summary.csv.

Core behavior:
    1. Recursively scans Ligases/<LIGASE>/SDF/*.sdf.
    2. Reads Ligand_SASA_summary.csv and builds the active Ligase+Ligand set.
    3. Any SDF whose Ligase+Ligand is NOT present in Ligand_SASA_summary.csv
       is moved into Retired_SDFs/.
    4. Generates SMILES directly from remaining active SDF files using RDKit.
    5. Optionally falls back to Open Babel if available.
    6. Computes RDKit metadata/descriptors.
    7. Validates:
         - SDF heavy atoms vs generated-SMILES heavy atoms.
         - SDF/SMILES heavy atoms vs SASA Heavy_atoms when present.
           Falls back to Total_atoms only if Heavy_atoms is unavailable.
    8. Writes:
         - Ligand_Metadata.csv
         - Ligase_Table/Ligase_SMILE_Codes.csv
         - Ligand_Metadata_Audit.csv
         - Ligand_Metadata_From_SDF.log
         - Retired_SDFs/Retired_SDF_Manifest.csv

Important:
    PASS_NO_PDB_SUMMARY_MATCH rows are no longer kept in the active metadata.
    They are retired because those SDFs do not correspond to current SASA/PDB data.

Usage:
    python 8_Ligand_Metadata_From_SDF.py --overwrite

Optional:
    python 8_Ligand_Metadata_From_SDF.py --overwrite --dry-run
    python 8_Ligand_Metadata_From_SDF.py --overwrite --retire-mismatches
    python 8_Ligand_Metadata_From_SDF.py --overwrite --no-openbabel
===============================================================================
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from tqdm import tqdm

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


# =============================================================================
# Output schemas
# =============================================================================

METADATA_COLUMNS = [
    "Ligand", "Name", "Formula", "Type",
    "SMILES", "Canonical_SMILES", "InChI", "InChIKey",
    "Formal_Charge", "Atom_Count", "Chiral_Atom_Count",
    "Bond_Count", "Aromatic_Bond_Count",

    "Ligase", "RECRUITER_CODE", "Source_SDF",
    "SDF_Title", "SMILES_Source",
    "SDF_Total_Atom_Count", "SDF_Heavy_Atom_Count",
    "SMILES_Atom_Count", "SMILES_Heavy_Atom_Count",
    "SDF_vs_SMILES_HeavyAtom_Match",
    "SASA_Atom_Count_Source",
    "SASA_Atom_Counts_Unique",
    "SASA_Atom_Count_Min", "SASA_Atom_Count_Max",
    "Atom_Count_Validation", "Validation_Notes",
]

SMILE_CODES_COLUMNS = [
    "RECRUITER_CODE", "Ligase", "Ligand",
    "SMILES", "Canonical_SMILES", "InChIKey",
    "Formula", "MW", "Source_SDF",
    "SDF_Heavy_Atom_Count", "SMILES_Heavy_Atom_Count",
    "SASA_Atom_Count_Source", "SASA_Atom_Counts_Unique",
    "Atom_Count_Validation", "Validation_Notes",
]

AUDIT_COLUMNS = [
    "RECRUITER_CODE", "Ligase", "Ligand", "Source_SDF",
    "Status", "SMILES_Source", "SMILES",
    "SDF_Total_Atom_Count", "SDF_Heavy_Atom_Count",
    "SMILES_Heavy_Atom_Count",
    "SASA_Atom_Count_Source", "SASA_Atom_Counts_Unique",
    "Atom_Count_Validation", "Validation_Notes", "Error",
]

RETIRED_COLUMNS = [
    "RECRUITER_CODE", "Ligase", "Ligand",
    "Original_SDF", "Retired_SDF",
    "Reason", "Timestamp",
]


# =============================================================================
# Utility helpers
# =============================================================================

def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def clean_str(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and math.isnan(x):
        return ""
    return str(x).strip()


def backup_if_exists(path: Path) -> Optional[Path]:
    if not path.exists():
        return None
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = path.with_name(f"{path.stem}.backup_{stamp}{path.suffix}")
    shutil.copy2(path, backup)
    return backup


def find_obabel() -> Optional[str]:
    return shutil.which("obabel") or shutil.which("babel")


def ensure_unique_destination(dest: Path) -> Path:
    """
    Avoid overwriting retired files if the script is run more than once.
    """
    if not dest.exists():
        return dest

    stem = dest.stem
    suffix = dest.suffix
    parent = dest.parent

    i = 1
    while True:
        candidate = parent / f"{stem}.retired_{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1


# =============================================================================
# Discovery / path parsing
# =============================================================================

def discover_sdfs(ligase_root: Path, retired_root: Path) -> List[Path]:
    """
    Find active SDFs under Ligases/<Ligase>/SDF/*.sdf.
    Explicitly ignores anything already under Retired_SDFs.
    """
    if not ligase_root.exists():
        raise FileNotFoundError(f"Ligase root not found: {ligase_root}")

    sdfs = []
    for p in ligase_root.rglob("SDF/*.sdf"):
        try:
            p.relative_to(retired_root)
            continue
        except Exception:
            pass
        sdfs.append(p)

    return sorted(sdfs)


def parse_ligase_and_code(sdf_path: Path) -> Tuple[str, str, str]:
    """
    Infer ligase and legacy ligand tag from:
        Ligases/CRBN/SDF/CRBN_A1A.sdf

    Returns:
        ligase='CRBN', ligand='A1A', recruiter_code='CRBN_A1A'
    """
    ligase = sdf_path.parent.parent.name.strip()
    stem = sdf_path.stem.strip()
    prefix = f"{ligase}_"

    if stem.startswith(prefix):
        ligand = stem[len(prefix):]
    else:
        ligand = stem.split("_")[-1]

    ligand = ligand.upper().strip()
    recruiter_code = f"{ligase}_{ligand}"
    return ligase, ligand, recruiter_code


# =============================================================================
# SASA summary map
# =============================================================================

def build_summary_atom_map(summary_path: Path) -> Tuple[Dict[Tuple[str, str], List[int]], str]:
    """
    Map (Ligase, Ligand) -> sorted unique atom-count values from SASA summary.

    Prefer Heavy_atoms when present.
    Fall back to Total_atoms only if Heavy_atoms is unavailable.
    """
    atom_map: Dict[Tuple[str, str], List[int]] = {}

    if not summary_path.exists():
        raise FileNotFoundError(f"Missing SASA summary file: {summary_path}")

    df = pd.read_csv(summary_path)

    required = {"Ligase", "Ligand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{summary_path} is missing required columns: {missing}")

    if "Heavy_atoms" in df.columns:
        count_col = "Heavy_atoms"
    elif "Total_atoms" in df.columns:
        count_col = "Total_atoms"
    else:
        raise ValueError(
            f"{summary_path} must contain either Heavy_atoms or Total_atoms for validation."
        )

    tmp = df.copy()
    tmp["Ligase"] = tmp["Ligase"].astype(str).str.strip()
    tmp["Ligand"] = tmp["Ligand"].astype(str).str.strip().str.upper()
    tmp["_validation_atom_count"] = pd.to_numeric(tmp[count_col], errors="coerce")

    tmp = tmp.dropna(subset=["_validation_atom_count"])

    for (ligase, ligand), g in tmp.groupby(["Ligase", "Ligand"]):
        vals = sorted({int(v) for v in g["_validation_atom_count"].tolist()})
        atom_map[(ligase, ligand)] = vals

    return atom_map, count_col


# =============================================================================
# SDF parsing / atom counts
# =============================================================================

def count_sdf_atoms_from_text(sdf_path: Path) -> Tuple[Optional[int], Optional[int]]:
    """
    Lightweight SDF count fallback.
    Returns:
        total_atom_count, heavy_atom_count
    """
    try:
        lines = sdf_path.read_text(errors="replace").splitlines()
    except Exception:
        return None, None

    for i, line in enumerate(lines[:80]):
        if "V2000" in line:
            total = None
            try:
                total = int(line[0:3])
            except Exception:
                parts = line.split()
                if parts and parts[0].isdigit():
                    total = int(parts[0])

            heavy = None
            if total is not None:
                atom_lines = lines[i + 1:i + 1 + total]
                heavy_count = 0
                for a in atom_lines:
                    elem = a[31:34].strip() if len(a) >= 34 else ""
                    if not elem:
                        parts = a.split()
                        elem = parts[3] if len(parts) >= 4 else ""
                    elem = re.sub(r"[^A-Za-z]", "", elem)
                    if elem and elem.upper() != "H":
                        heavy_count += 1
                heavy = heavy_count

            return total, heavy

        if "V3000" in line:
            for j in range(i, min(i + 150, len(lines))):
                if "COUNTS" in lines[j]:
                    parts = lines[j].split()
                    try:
                        idx = parts.index("COUNTS")
                        total = int(parts[idx + 1])
                        return total, None
                    except Exception:
                        return None, None

    return None, None


def load_sdf_mol(sdf_path: Path) -> Tuple[Optional[Chem.Mol], str, str]:
    """
    Load first molecule from SDF using tolerant RDKit parsing.
    """
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
        mol = suppl[0] if suppl and len(suppl) else None
        if mol is None:
            return None, "", "RDKit returned no molecule from SDF"

        title = mol.GetProp("_Name") if mol.HasProp("_Name") else sdf_path.stem

        try:
            mol2 = Chem.Mol(mol)
            Chem.SanitizeMol(mol2)
            mol = mol2
        except Exception:
            try:
                mol2 = Chem.Mol(mol)
                Chem.SanitizeMol(
                    mol2,
                    sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                )
                mol = mol2
            except Exception:
                pass

        return mol, title, ""

    except Exception as e:
        return None, "", f"RDKit SDF load failed: {e}"


# =============================================================================
# SMILES generation
# =============================================================================

def smiles_from_rdkit_mol(mol: Chem.Mol) -> Tuple[Optional[str], Optional[str], str]:
    try:
        try:
            mol_no_h = Chem.RemoveHs(mol, sanitize=False)
        except Exception:
            mol_no_h = mol

        smi = Chem.MolToSmiles(mol_no_h, canonical=True, isomericSmiles=True)
        if not smi:
            return None, None, "RDKit produced empty SMILES"

        parsed = Chem.MolFromSmiles(smi)
        if parsed is not None:
            can = Chem.MolToSmiles(parsed, canonical=True, isomericSmiles=True)
        else:
            can = smi

        return smi, can, ""

    except Exception as e:
        return None, None, f"RDKit SMILES generation failed: {e}"


def smiles_from_openbabel(sdf_path: Path, obabel_path: str) -> Tuple[Optional[str], str]:
    try:
        cmd = [obabel_path, str(sdf_path), "-ocan"]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        output = (p.stdout or "").strip()

        if p.returncode != 0 or not output:
            err = (p.stderr or "").strip()
            return None, f"OpenBabel failed: {err}"

        first = output.splitlines()[0].strip()
        smi = first.split()[0].strip()
        return smi, ""

    except Exception as e:
        return None, f"OpenBabel exception: {e}"


def mol_from_smiles_tolerant(smiles: str) -> Tuple[Optional[Chem.Mol], str]:
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is not None:
            return mol, ""
    except Exception as e:
        first_error = str(e)
    else:
        first_error = "MolFromSmiles returned None"

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is not None:
            return mol, "parsed sanitize=False"
    except Exception as e:
        return None, f"SMILES parse failed: {first_error}; fallback failed: {e}"

    return None, f"SMILES parse failed: {first_error}"


# =============================================================================
# Metadata / validation
# =============================================================================

def safe_inchi(mol: Chem.Mol) -> Tuple[str, str]:
    try:
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi) if inchi else ""
        return inchi or "", inchikey or ""
    except Exception:
        return "", ""


def compute_metadata(smiles: str) -> Tuple[Dict[str, object], Optional[Chem.Mol], str]:
    mol, parse_note = mol_from_smiles_tolerant(smiles)
    if mol is None:
        return {}, None, parse_note

    try:
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        canonical = smiles

    inchi, inchikey = safe_inchi(mol)

    try:
        formula = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        formula = ""

    try:
        chiral_count = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    except Exception:
        chiral_count = ""

    try:
        aromatic_bonds = sum(1 for b in mol.GetBonds() if b.GetIsAromatic())
    except Exception:
        aromatic_bonds = ""

    try:
        mw = Descriptors.MolWt(mol)
    except Exception:
        mw = ""

    meta = {
        "Name": "",
        "Formula": formula,
        "Type": "SDF-derived recruiter",
        "SMILES": smiles,
        "Canonical_SMILES": canonical,
        "InChI": inchi,
        "InChIKey": inchikey,
        "Formal_Charge": Chem.GetFormalCharge(mol),
        "Atom_Count": mol.GetNumAtoms(),
        "Chiral_Atom_Count": chiral_count,
        "Bond_Count": mol.GetNumBonds(),
        "Aromatic_Bond_Count": aromatic_bonds,
        "MW": mw,
    }

    return meta, mol, parse_note


def validate_counts(
    sdf_heavy: Optional[int],
    smiles_heavy: Optional[int],
    sasa_counts: List[int],
) -> Tuple[str, str, bool]:
    notes = []

    if sdf_heavy is not None and smiles_heavy is not None:
        sdf_smiles_match = int(sdf_heavy) == int(smiles_heavy)
    else:
        sdf_smiles_match = False
        notes.append("Could not compare SDF heavy atoms to SMILES heavy atoms")

    if not sdf_smiles_match:
        notes.append(f"SDF heavy atoms ({sdf_heavy}) != SMILES heavy atoms ({smiles_heavy})")
        status = "FAIL_SDF_SMILES_MISMATCH"
    else:
        status = "PASS"

    if sasa_counts:
        if len(sasa_counts) > 1:
            notes.append(f"Multiple SASA atom-count values observed: {sasa_counts}")

        comparable = []
        if sdf_heavy is not None:
            comparable.append(int(sdf_heavy))
        if smiles_heavy is not None:
            comparable.append(int(smiles_heavy))

        if comparable and any(v in sasa_counts for v in comparable):
            if status == "PASS":
                status = "PASS"
            else:
                status = "WARN_SASA_MATCH_BUT_SDF_SMILES_MISMATCH"
        else:
            notes.append(
                f"No SDF/SMILES heavy atom count matched SASA validation atom counts {sasa_counts}"
            )
            if status == "PASS":
                status = "FAIL_PDB_ATOM_COUNT_MISMATCH"
            else:
                status = f"{status};FAIL_PDB_ATOM_COUNT_MISMATCH"
    else:
        # In this replacement script, this should never appear in active rows,
        # because no-SASA-match SDFs are retired before processing.
        notes.append("No Ligase+Ligand match found in Ligand_SASA_summary.csv")
        if status == "PASS":
            status = "PASS_NO_PDB_SUMMARY_MATCH"

    return status, "; ".join(notes), sdf_smiles_match


# =============================================================================
# Retirement logic
# =============================================================================

def retire_sdf(
    sdf_path: Path,
    retired_root: Path,
    ligase: str,
    recruiter_code: str,
    reason: str,
    dry_run: bool = False,
) -> Path:
    """
    Move an SDF to:
        Retired_SDFs/<Ligase>/SDF/<filename>
    """
    dest_dir = retired_root / ligase / "SDF"
    dest_dir.mkdir(parents=True, exist_ok=True)

    dest = ensure_unique_destination(dest_dir / sdf_path.name)

    if not dry_run:
        shutil.move(str(sdf_path), str(dest))

    return dest


# =============================================================================
# Per-SDF processing
# =============================================================================

def process_sdf(
    sdf_path: Path,
    atom_map: Dict[Tuple[str, str], List[int]],
    atom_count_source: str,
    use_openbabel: bool,
    obabel_path: Optional[str],
) -> Tuple[Optional[Dict[str, object]], Dict[str, object]]:
    ligase, ligand, recruiter_code = parse_ligase_and_code(sdf_path)
    source_sdf = str(sdf_path)

    audit = {
        "RECRUITER_CODE": recruiter_code,
        "Ligase": ligase,
        "Ligand": ligand,
        "Source_SDF": source_sdf,
        "Status": "FAILED",
        "SMILES_Source": "",
        "SMILES": "",
        "SDF_Total_Atom_Count": "",
        "SDF_Heavy_Atom_Count": "",
        "SMILES_Heavy_Atom_Count": "",
        "SASA_Atom_Count_Source": atom_count_source,
        "SASA_Atom_Counts_Unique": "",
        "Atom_Count_Validation": "",
        "Validation_Notes": "",
        "Error": "",
    }

    sdf_total_text, sdf_heavy_text = count_sdf_atoms_from_text(sdf_path)
    sdf_mol, sdf_title, sdf_error = load_sdf_mol(sdf_path)

    sdf_total = sdf_total_text
    sdf_heavy = sdf_heavy_text

    if sdf_mol is not None:
        try:
            sdf_total = sdf_mol.GetNumAtoms()
            sdf_heavy = sum(1 for a in sdf_mol.GetAtoms() if a.GetAtomicNum() > 1)
        except Exception:
            pass

    audit["SDF_Total_Atom_Count"] = sdf_total if sdf_total is not None else ""
    audit["SDF_Heavy_Atom_Count"] = sdf_heavy if sdf_heavy is not None else ""

    errors = []
    if sdf_error:
        errors.append(sdf_error)

    smiles = None
    canonical = None
    smiles_source = ""

    if sdf_mol is not None:
        smiles, canonical, err = smiles_from_rdkit_mol(sdf_mol)
        if smiles:
            smiles_source = "RDKit_SDF"
        elif err:
            errors.append(err)

    if not smiles and use_openbabel and obabel_path:
        smiles, err = smiles_from_openbabel(sdf_path, obabel_path)
        if smiles:
            smiles_source = "OpenBabel_SDF"
            canonical = smiles
        elif err:
            errors.append(err)

    if not smiles:
        audit["Error"] = "; ".join(errors) if errors else "Could not generate SMILES"
        return None, audit

    meta, smiles_mol, parse_note = compute_metadata(smiles)
    if not meta or smiles_mol is None:
        audit["SMILES"] = smiles
        audit["SMILES_Source"] = smiles_source
        audit["Error"] = parse_note or "Could not parse generated SMILES"
        return None, audit

    smiles_heavy = sum(1 for a in smiles_mol.GetAtoms() if a.GetAtomicNum() > 1)
    sasa_counts = atom_map.get((ligase, ligand), [])

    status, notes, sdf_smiles_match = validate_counts(
        sdf_heavy=sdf_heavy,
        smiles_heavy=smiles_heavy,
        sasa_counts=sasa_counts,
    )

    row = {col: "" for col in METADATA_COLUMNS}
    row.update({k: v for k, v in meta.items() if k in row})

    row.update({
        "Ligand": ligand,
        "Ligase": ligase,
        "RECRUITER_CODE": recruiter_code,
        "Source_SDF": source_sdf,
        "SDF_Title": sdf_title,
        "SMILES_Source": smiles_source,
        "SDF_Total_Atom_Count": sdf_total if sdf_total is not None else "",
        "SDF_Heavy_Atom_Count": sdf_heavy if sdf_heavy is not None else "",
        "SMILES_Atom_Count": smiles_mol.GetNumAtoms(),
        "SMILES_Heavy_Atom_Count": smiles_heavy,
        "SDF_vs_SMILES_HeavyAtom_Match": sdf_smiles_match,
        "SASA_Atom_Count_Source": atom_count_source,
        "SASA_Atom_Counts_Unique": ";".join(map(str, sasa_counts)),
        "SASA_Atom_Count_Min": min(sasa_counts) if sasa_counts else "",
        "SASA_Atom_Count_Max": max(sasa_counts) if sasa_counts else "",
        "Atom_Count_Validation": status,
        "Validation_Notes": notes,
    })

    audit.update({
        "Status": "OK" if status.startswith("PASS") else "CHECK",
        "SMILES_Source": smiles_source,
        "SMILES": smiles,
        "SMILES_Heavy_Atom_Count": smiles_heavy,
        "SASA_Atom_Counts_Unique": ";".join(map(str, sasa_counts)),
        "Atom_Count_Validation": status,
        "Validation_Notes": notes,
        "Error": "; ".join(errors),
    })

    return row, audit


def make_smile_codes(metadata_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in metadata_df.iterrows():
        rows.append({
            "RECRUITER_CODE": r.get("RECRUITER_CODE", ""),
            "Ligase": r.get("Ligase", ""),
            "Ligand": r.get("Ligand", ""),
            "SMILES": r.get("SMILES", ""),
            "Canonical_SMILES": r.get("Canonical_SMILES", ""),
            "InChIKey": r.get("InChIKey", ""),
            "Formula": r.get("Formula", ""),
            "MW": r.get("MW", ""),
            "Source_SDF": r.get("Source_SDF", ""),
            "SDF_Heavy_Atom_Count": r.get("SDF_Heavy_Atom_Count", ""),
            "SMILES_Heavy_Atom_Count": r.get("SMILES_Heavy_Atom_Count", ""),
            "SASA_Atom_Count_Source": r.get("SASA_Atom_Count_Source", ""),
            "SASA_Atom_Counts_Unique": r.get("SASA_Atom_Counts_Unique", ""),
            "Atom_Count_Validation": r.get("Atom_Count_Validation", ""),
            "Validation_Notes": r.get("Validation_Notes", ""),
        })

    return pd.DataFrame(rows, columns=SMILE_CODES_COLUMNS)


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build active ligand metadata from local SDFs and retire SDFs with no current SASA/PDB match."
    )

    parser.add_argument("--ligase-root", default="Ligases")
    parser.add_argument("--summary", default="Ligand_SASA_summary.csv")
    parser.add_argument("--metadata-out", default="Ligand_Metadata.csv")
    parser.add_argument("--smile-codes-out", default="Ligase_Table/Ligase_SMILE_Codes.csv")
    parser.add_argument("--audit-out", default="Ligand_Metadata_Audit.csv")
    parser.add_argument("--log", default="Ligand_Metadata_From_SDF.log")
    parser.add_argument("--retired-root", default="Retired_SDFs")
    parser.add_argument("--retired-manifest", default="Retired_SDFs/Retired_SDF_Manifest.csv")

    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-openbabel", action="store_true")
    parser.add_argument("--strict", action="store_true")

    parser.add_argument(
        "--retire-mismatches",
        action="store_true",
        help="Also move atom-count mismatch SDFs into Retired_SDFs. Default is false.",
    )

    args = parser.parse_args()

    ligase_root = Path(args.ligase_root)
    summary_path = Path(args.summary)
    metadata_out = Path(args.metadata_out)
    smile_codes_out = Path(args.smile_codes_out)
    audit_out = Path(args.audit_out)
    log_path = Path(args.log)
    retired_root = Path(args.retired_root)
    retired_manifest = Path(args.retired_manifest)

    use_openbabel = not args.no_openbabel
    obabel_path = find_obabel() if use_openbabel else None

    if use_openbabel and not obabel_path:
        print("⚠️ Open Babel fallback requested, but obabel/babel was not found in PATH. RDKit-only mode will be used.")

    atom_map, atom_count_source = build_summary_atom_map(summary_path)
    active_keys = set(atom_map.keys())

    all_sdfs = discover_sdfs(ligase_root, retired_root)

    if not all_sdfs:
        print(f"❌ No SDF files found under {ligase_root}/<Ligase>/SDF/*.sdf")
        return 2

    if not args.overwrite:
        blocked = [p for p in [metadata_out, smile_codes_out, audit_out, log_path] if p.exists()]
        if blocked:
            print("❌ Output file(s) already exist:")
            for p in blocked:
                print(f"   - {p}")
            print("Run with --overwrite to back them up and regenerate cleanly.")
            return 2

    for p in [metadata_out, smile_codes_out, audit_out, log_path, retired_manifest]:
        p.parent.mkdir(parents=True, exist_ok=True)
        if args.overwrite and p.exists() and not args.dry_run:
            backup = backup_if_exists(p)
            if backup:
                print(f"🧷 Backed up {p} -> {backup}")

    print(f"📂 Found {len(all_sdfs)} SDF files under {ligase_root}")
    print(f"📊 Loaded SASA validation map for {len(atom_map)} Ligase+Ligand pairs from {summary_path}")
    print(f"🧮 SASA atom-count validation column: {atom_count_source}")
    print(f"🧪 SMILES source priority: RDKit SDF" + (" -> OpenBabel SDF fallback" if obabel_path else ""))
    print(f"🗃️ Retired SDF root: {retired_root}")
    if args.dry_run:
        print("🧪 DRY RUN: no files will be moved or written.")

    # -------------------------------------------------------------------------
    # Phase 1: retire no-SASA-match SDFs before metadata generation
    # -------------------------------------------------------------------------
    active_sdfs = []
    retired_rows = []

    for sdf_path in all_sdfs:
        ligase, ligand, recruiter_code = parse_ligase_and_code(sdf_path)
        key = (ligase, ligand)

        if key not in active_keys:
            reason = "NO_LIGASE_LIGAND_MATCH_IN_CURRENT_SASA_SUMMARY"
            dest = retire_sdf(
                sdf_path=sdf_path,
                retired_root=retired_root,
                ligase=ligase,
                recruiter_code=recruiter_code,
                reason=reason,
                dry_run=args.dry_run,
            )

            retired_rows.append({
                "RECRUITER_CODE": recruiter_code,
                "Ligase": ligase,
                "Ligand": ligand,
                "Original_SDF": str(sdf_path),
                "Retired_SDF": str(dest),
                "Reason": reason,
                "Timestamp": now(),
            })
        else:
            active_sdfs.append(sdf_path)

    print(f"\n🧹 Retired no-SASA-match SDFs: {len(retired_rows)}")
    print(f"✅ Active SDFs remaining for metadata: {len(active_sdfs)}")

    # -------------------------------------------------------------------------
    # Phase 2: process active SDFs
    # -------------------------------------------------------------------------
    rows: List[Dict[str, object]] = []
    audit_rows: List[Dict[str, object]] = []

    for sdf_path in tqdm(active_sdfs, desc="Processing active SDF ligands"):
        row, audit = process_sdf(
            sdf_path=sdf_path,
            atom_map=atom_map,
            atom_count_source=atom_count_source,
            use_openbabel=use_openbabel,
            obabel_path=obabel_path,
        )

        # Optional: retire mismatches too.
        if row and args.retire_mismatches:
            validation = clean_str(row.get("Atom_Count_Validation"))
            if not validation.startswith("PASS"):
                ligase = clean_str(row.get("Ligase"))
                ligand = clean_str(row.get("Ligand"))
                recruiter_code = clean_str(row.get("RECRUITER_CODE"))
                reason = f"ATOM_COUNT_VALIDATION_{validation}"

                dest = retire_sdf(
                    sdf_path=sdf_path,
                    retired_root=retired_root,
                    ligase=ligase,
                    recruiter_code=recruiter_code,
                    reason=reason,
                    dry_run=args.dry_run,
                )

                retired_rows.append({
                    "RECRUITER_CODE": recruiter_code,
                    "Ligase": ligase,
                    "Ligand": ligand,
                    "Original_SDF": str(sdf_path),
                    "Retired_SDF": str(dest),
                    "Reason": reason,
                    "Timestamp": now(),
                })

                # Do not include retired mismatches in active metadata.
                audit["Status"] = "RETIRED"
                audit["Validation_Notes"] = (
                    clean_str(audit.get("Validation_Notes"))
                    + "; Retired due to --retire-mismatches"
                ).strip("; ")

                audit_rows.append(audit)
                continue

        audit_rows.append(audit)
        if row:
            rows.append(row)

    metadata_df = pd.DataFrame(rows, columns=METADATA_COLUMNS)
    audit_df = pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS)
    smile_codes_df = make_smile_codes(metadata_df)
    retired_df = pd.DataFrame(retired_rows, columns=RETIRED_COLUMNS)

    if not args.dry_run:
        metadata_df.to_csv(metadata_out, index=False)
        audit_df.to_csv(audit_out, index=False)
        smile_codes_df.to_csv(smile_codes_out, index=False)
        retired_df.to_csv(retired_manifest, index=False)

    validation_counts = (
        metadata_df["Atom_Count_Validation"].value_counts(dropna=False).to_dict()
        if not metadata_df.empty else {}
    )

    status_counts = (
        audit_df["Status"].value_counts(dropna=False).to_dict()
        if not audit_df.empty else {}
    )

    failed_validation = (
        metadata_df[
            ~metadata_df["Atom_Count_Validation"].astype(str).str.startswith("PASS", na=False)
        ]
        if not metadata_df.empty else metadata_df
    )

    log_lines = []
    log_lines.append("=== Ligand metadata from active SDF run ===")
    log_lines.append(f"Started/finished: {now()}")
    log_lines.append(f"Ligase root: {ligase_root}")
    log_lines.append(f"SASA summary: {summary_path}")
    log_lines.append(f"SASA atom-count source column: {atom_count_source}")
    log_lines.append(f"SDF files discovered: {len(all_sdfs)}")
    log_lines.append(f"No-SASA-match SDFs retired: {len(retired_rows)}")
    log_lines.append(f"Active SDFs processed: {len(active_sdfs)}")
    log_lines.append(f"Metadata rows written: {len(metadata_df)}")
    log_lines.append(f"SMILES table rows written: {len(smile_codes_df)}")
    log_lines.append(f"Audit rows written: {len(audit_df)}")
    log_lines.append(f"OpenBabel fallback enabled: {bool(obabel_path)} ({obabel_path or 'not found/disabled'})")
    log_lines.append(f"Validation counts: {json.dumps(validation_counts, indent=2)}")
    log_lines.append(f"Audit status counts: {json.dumps(status_counts, indent=2)}")
    log_lines.append(f"Retired manifest: {retired_manifest}")

    if len(failed_validation) > 0:
        log_lines.append("")
        log_lines.append("=== Active rows requiring review ===")
        for _, r in failed_validation.iterrows():
            log_lines.append(
                f"{r.get('RECRUITER_CODE')}\t{r.get('Atom_Count_Validation')}\t"
                f"SDF_H={r.get('SDF_Heavy_Atom_Count')}\t"
                f"SMI_H={r.get('SMILES_Heavy_Atom_Count')}\t"
                f"SASA={r.get('SASA_Atom_Counts_Unique')}\t"
                f"{r.get('Validation_Notes')}"
            )

    if not args.dry_run:
        with open(log_path, "w") as log:
            log.write("\n".join(log_lines) + "\n")

    print(f"\n✅ Written {metadata_out} ({len(metadata_df)} rows)" if not args.dry_run else f"\n🧪 Would write {metadata_out} ({len(metadata_df)} rows)")
    print(f"✅ Written {smile_codes_out} ({len(smile_codes_df)} rows)" if not args.dry_run else f"🧪 Would write {smile_codes_out} ({len(smile_codes_df)} rows)")
    print(f"✅ Written {audit_out} ({len(audit_df)} rows)" if not args.dry_run else f"🧪 Would write {audit_out} ({len(audit_df)} rows)")
    print(f"✅ Written {retired_manifest} ({len(retired_df)} rows)" if not args.dry_run else f"🧪 Would write {retired_manifest} ({len(retired_df)} rows)")
    print(f"🧾 Written {log_path}" if not args.dry_run else f"🧪 Would write {log_path}")

    print("\n📊 Atom-count validation summary:")
    for k, v in validation_counts.items():
        print(f"  {k}: {v}")

    print("\n📦 Retirement summary:")
    print(f"  Retired no-SASA-match SDFs: {len(retired_rows)}")
    print(f"  Active metadata rows: {len(metadata_df)}")

    if len(failed_validation) > 0:
        print(f"\n⚠️ {len(failed_validation)} active metadata rows still require review.")
        print("   These have SASA/PDB matches, so they were NOT retired by default.")
        print("   Use --retire-mismatches only if you want those moved too.")
        if args.strict:
            return 1
    else:
        print("\n🎯 All active metadata rows passed atom-count validation.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())