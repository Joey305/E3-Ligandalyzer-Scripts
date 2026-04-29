#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
12_LigaseChemicalDescriptors.py
===============================================================================
Purpose:
    Generate chemical descriptor tables for LR00001-style recruiter codes.

Input:
    Ligase_Table/Ligase_SMILE_Codes.csv
        Required columns:
            SMILES
            RECRUITER_CODE

Existing optional input:
    Ligase_Table/Ligase_SMILE_Codes_Atoms.csv
        Created earlier by 10_MCS_Matching.py / 10B_CLean.py.
        This script validates it by default but does NOT overwrite it unless
        --write-atoms is supplied.

Outputs:
    Ligase_Table/Ligase_Chemical_Descriptors.csv
    Ligase_Table/Ligase_Chemical_Descriptors_Audit.csv
    Ligase_Table/Ligase_SMILE_Codes_Atoms_Validation.csv
    Ligase_Table/Ligase_Chemical_Descriptors.log

Optional output with --write-atoms:
    Ligase_Table/Ligase_SMILE_Codes_Atoms.csv

Descriptors include:
    - MW / ExactMW / Formula
    - LogP, TPSA, HBA, HBD
    - Rotatable bonds, rings, aromatic rings, CSP3
    - Chiral atoms, formal charge
    - QED
    - Synthetic Accessibility score when sascorer is available
    - PAINS / BRENK alerts
    - Reactive SMARTS alerts
    - Lipinski / Veber / Egan / Ghose / Muegge filters

Usage:
    python 12_LigaseChemicalDescriptors.py --dry-run
    python 12_LigaseChemicalDescriptors.py --overwrite
    python 12_LigaseChemicalDescriptors.py --overwrite --write-atoms
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

RDLogger.DisableLog("rdApp.*")


# =============================================================================
# Defaults
# =============================================================================

DEFAULT_INPUT = Path("Ligase_Table/Ligase_SMILE_Codes.csv")
DEFAULT_ATOMS_OUT = Path("Ligase_Table/Ligase_SMILE_Codes_Atoms.csv")
DEFAULT_DESC_OUT = Path("Ligase_Table/Ligase_Chemical_Descriptors.csv")
DEFAULT_AUDIT_OUT = Path("Ligase_Table/Ligase_Chemical_Descriptors_Audit.csv")
DEFAULT_ATOM_VALIDATION_OUT = Path("Ligase_Table/Ligase_SMILE_Codes_Atoms_Validation.csv")
DEFAULT_LOG_FILE = Path("Ligase_Table/Ligase_Chemical_Descriptors.log")

LR_PATTERN = re.compile(r"^LR\d{5}$")


# =============================================================================
# Optional SA scorer
# =============================================================================

try:
    import sascorer
    HAS_SA = True
except Exception:
    HAS_SA = False


# =============================================================================
# Reactive SMARTS patterns
# =============================================================================

REACTIVE_SMARTS = {
    "Acid chloride": "[CX3](=O)Cl",
    "Sulfonyl chloride": "S(=O)(=O)Cl",
    "Isocyanate": "N=C=O",
    "Isothiocyanate": "N=C=S",
    "Aldehyde": "[CX3H1](=O)[#6]",
    "Epoxide": "[OX2r3]1[CX4r3][CX4r3]1",
    "Anhydride": "[CX3](=O)O[CX3](=O)",
    "Peroxide": "[OX2][OX2]",
    "Thiol": "[SX2H]",
    "Hydrazine": "[NX3][NX3]",
    "Alkyl halide": "[CX4][Cl,Br,I]",
    "Michael acceptor": "[C;H1,H2]=[C;H0,H1][CX3](=O)",
    "Diazo": "[N-]=[N+]=[*]",
}

COMPILED_REACTIVE_SMARTS = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in REACTIVE_SMARTS.items()
    if Chem.MolFromSmarts(smarts) is not None
}


# =============================================================================
# Output schemas
# =============================================================================

ATOM_COLUMNS = [
    "RECRUITER_CODE",
    "smiles_atom_index",
    "smile_atom",
]

DESCRIPTOR_COLUMNS = [
    "RECRUITER_CODE",
    "SMILES",
    "Canonical_SMILES",
    "Formula",
    "MW",
    "Exact_MW",
    "LogP",
    "TPSA",
    "HBA",
    "HBD",
    "Rotatable_Bonds",
    "Ring_Count",
    "Aromatic_Rings",
    "Fraction_CSP3",
    "Heavy_Atom_Count",
    "Atom_Count",
    "Chiral_Atoms",
    "Formal_Charge",
    "QED",
    "SA_Score",
    "BertzCT",
    "HallKierAlpha",
    "Kappa1",
    "Kappa2",
    "Kappa3",
    "NumSpiroAtoms",
    "NumBridgeheadAtoms",
    "NumAliphaticRings",
    "NumAromaticRings",
    "NumSaturatedRings",
    "NumHeteroAtoms",
    "MolMR",
    "Lipinski_Violations",
    "Lipinski_Pass",
    "Veber_Pass",
    "Egan_Pass",
    "Ghose_Pass",
    "Muegge_Pass",
    "PAINS_Hits",
    "Brenk_Hits",
    "Reactive_SMARTS_Hits",
]

AUDIT_COLUMNS = [
    "RECRUITER_CODE",
    "SMILES",
    "Status",
    "Notes",
    "Error",
]

ATOM_VALIDATION_COLUMNS = [
    "RECRUITER_CODE",
    "Expected_SMILES_Atom_Count",
    "Observed_Atom_Table_Count",
    "Min_Atom_Index",
    "Max_Atom_Index",
    "Duplicate_Atom_Indices",
    "Status",
    "Notes",
]


# =============================================================================
# Helpers
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


def backup_if_exists(path: Path) -> Optional[Path]:
    if not path.exists():
        return None
    backup = path.with_name(f"{path.stem}.backup_{stamp()}{path.suffix}")
    shutil.copy2(path, backup)
    return backup


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")


def write_csv(path: Path, columns: List[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(path, index=False)


# =============================================================================
# RDKit molecule parsing
# =============================================================================

def safe_mol(smiles: str) -> Tuple[Optional[Chem.Mol], str]:
    """
    Return molecule with tolerant fallback parsing.
    """
    smi = clean_str(smiles).replace('"', "").replace("'", "")
    if not smi:
        return None, "Empty SMILES"

    try:
        mol = Chem.MolFromSmiles(smi, sanitize=True)
        if mol is not None:
            return mol, ""
    except Exception as e:
        first_error = f"sanitize=True failed: {type(e).__name__}: {e}"
    else:
        first_error = "MolFromSmiles returned None"

    try:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is not None:
            try:
                mol.UpdatePropertyCache(strict=False)
                Chem.GetSymmSSSR(mol)
            except Exception:
                pass
            return mol, first_error + "; parsed sanitize=False"
    except Exception as e:
        return None, first_error + f"; sanitize=False failed: {type(e).__name__}: {e}"

    return None, first_error


def canonicalize_smiles(mol: Chem.Mol, fallback: str) -> str:
    try:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return fallback


def safe_formula(mol: Chem.Mol) -> str:
    try:
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return ""


def safe_exact_mw(mol: Chem.Mol) -> float:
    try:
        return float(rdMolDescriptors.CalcExactMolWt(mol))
    except Exception:
        return math.nan


def safe_bertz(mol: Chem.Mol) -> float:
    try:
        return float(rdMolDescriptors.CalcBertzCT(mol))
    except Exception:
        try:
            return float(Descriptors.BertzCT(mol))
        except Exception:
            return math.nan


def safe_descriptor(fn, mol: Chem.Mol, default=math.nan):
    try:
        return fn(mol)
    except Exception:
        return default


# =============================================================================
# Alerts
# =============================================================================

def build_filter_catalog(catalogs: List[FilterCatalogParams.FilterCatalogs]):
    params = FilterCatalogParams()
    for cat in catalogs:
        params.AddCatalog(cat)
    return FilterCatalog.FilterCatalog(params)


def check_catalog(mol: Chem.Mol, catalog) -> str:
    try:
        entries = catalog.GetMatches(mol)
        return "; ".join(sorted({e.GetDescription() for e in entries})) if entries else "None"
    except Exception:
        return "CATALOG_CHECK_FAILED"


def check_reactive_smarts(mol: Chem.Mol) -> str:
    hits = []
    for name, patt in COMPILED_REACTIVE_SMARTS.items():
        try:
            if mol.HasSubstructMatch(patt):
                hits.append(name)
        except Exception:
            continue
    return "; ".join(sorted(hits)) if hits else "None"


# =============================================================================
# Drug-likeness rules
# =============================================================================

def druglikeness_rules(d: Dict[str, object]) -> Dict[str, object]:
    mw = float(d.get("MW", math.nan))
    logp = float(d.get("LogP", math.nan))
    tpsa = float(d.get("TPSA", math.nan))
    hbd = int(d.get("HBD", 0))
    hba = int(d.get("HBA", 0))
    rot = int(d.get("Rotatable_Bonds", 0))
    heavy = int(d.get("Heavy_Atom_Count", 0))

    lipinski_violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10,
    ])

    lipinski_pass = lipinski_violations <= 1
    veber_pass = (tpsa <= 140) and (rot <= 10)
    egan_pass = (tpsa <= 131) and (logp <= 5.88)
    ghose_pass = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and (20 <= heavy <= 70)
    muegge_pass = (200 <= mw <= 600) and (-2 <= logp <= 5) and (hba <= 10) and (hbd <= 5) and (tpsa <= 150)

    return {
        "Lipinski_Violations": lipinski_violations,
        "Lipinski_Pass": lipinski_pass,
        "Veber_Pass": veber_pass,
        "Egan_Pass": egan_pass,
        "Ghose_Pass": ghose_pass,
        "Muegge_Pass": muegge_pass,
    }


# =============================================================================
# Descriptor calculation
# =============================================================================

def compute_descriptors(
    code: str,
    smiles: str,
    mol: Chem.Mol,
    pains_catalog,
    brenk_catalog,
) -> Dict[str, object]:

    canonical = canonicalize_smiles(mol, smiles)

    d = {
        "RECRUITER_CODE": code,
        "SMILES": smiles,
        "Canonical_SMILES": canonical,
        "Formula": safe_formula(mol),
        "MW": safe_descriptor(Descriptors.MolWt, mol),
        "Exact_MW": safe_exact_mw(mol),
        "LogP": safe_descriptor(Crippen.MolLogP, mol),
        "TPSA": safe_descriptor(rdMolDescriptors.CalcTPSA, mol),
        "HBA": safe_descriptor(rdMolDescriptors.CalcNumHBA, mol, 0),
        "HBD": safe_descriptor(rdMolDescriptors.CalcNumHBD, mol, 0),
        "Rotatable_Bonds": safe_descriptor(rdMolDescriptors.CalcNumRotatableBonds, mol, 0),
        "Ring_Count": safe_descriptor(rdMolDescriptors.CalcNumRings, mol, 0),
        "Aromatic_Rings": safe_descriptor(rdMolDescriptors.CalcNumAromaticRings, mol, 0),
        "Fraction_CSP3": safe_descriptor(rdMolDescriptors.CalcFractionCSP3, mol),
        "Heavy_Atom_Count": mol.GetNumHeavyAtoms(),
        "Atom_Count": mol.GetNumAtoms(),
        "Chiral_Atoms": safe_descriptor(rdMolDescriptors.CalcNumAtomStereoCenters, mol, 0),
        "Formal_Charge": Chem.GetFormalCharge(mol),
        "QED": safe_descriptor(QED.qed, mol),
        "BertzCT": safe_bertz(mol),
        "HallKierAlpha": safe_descriptor(Descriptors.HallKierAlpha, mol),
        "Kappa1": safe_descriptor(Descriptors.Kappa1, mol),
        "Kappa2": safe_descriptor(Descriptors.Kappa2, mol),
        "Kappa3": safe_descriptor(Descriptors.Kappa3, mol),
        "NumSpiroAtoms": safe_descriptor(rdMolDescriptors.CalcNumSpiroAtoms, mol, 0),
        "NumBridgeheadAtoms": safe_descriptor(rdMolDescriptors.CalcNumBridgeheadAtoms, mol, 0),
        "NumAliphaticRings": safe_descriptor(rdMolDescriptors.CalcNumAliphaticRings, mol, 0),
        "NumAromaticRings": safe_descriptor(rdMolDescriptors.CalcNumAromaticRings, mol, 0),
        "NumSaturatedRings": safe_descriptor(rdMolDescriptors.CalcNumSaturatedRings, mol, 0),
        "NumHeteroAtoms": safe_descriptor(rdMolDescriptors.CalcNumHeteroatoms, mol, 0),
        "MolMR": safe_descriptor(Descriptors.MolMR, mol),
    }

    if HAS_SA:
        try:
            d["SA_Score"] = sascorer.calculateScore(mol)
        except Exception:
            d["SA_Score"] = math.nan
    else:
        d["SA_Score"] = math.nan

    d.update(druglikeness_rules(d))

    d["PAINS_Hits"] = check_catalog(mol, pains_catalog)
    d["Brenk_Hits"] = check_catalog(mol, brenk_catalog)
    d["Reactive_SMARTS_Hits"] = check_reactive_smarts(mol)

    return d


def make_atom_rows(code: str, mol: Chem.Mol) -> List[Dict[str, object]]:
    rows = []
    for atom in mol.GetAtoms():
        rows.append({
            "RECRUITER_CODE": code,
            "smiles_atom_index": atom.GetIdx(),
            "smile_atom": atom.GetSymbol(),
        })
    return rows


# =============================================================================
# Atom table validation
# =============================================================================

def validate_existing_atom_table(
    atom_table_path: Path,
    expected_counts: Dict[str, int],
) -> List[Dict[str, object]]:

    rows = []

    if not atom_table_path.exists():
        for code, expected in sorted(expected_counts.items()):
            rows.append({
                "RECRUITER_CODE": code,
                "Expected_SMILES_Atom_Count": expected,
                "Observed_Atom_Table_Count": "",
                "Min_Atom_Index": "",
                "Max_Atom_Index": "",
                "Duplicate_Atom_Indices": "",
                "Status": "ATOM_TABLE_MISSING",
                "Notes": f"{atom_table_path} does not exist",
            })
        return rows

    atoms = pd.read_csv(atom_table_path)

    required = {"RECRUITER_CODE", "smiles_atom_index", "smile_atom"}
    missing = required - set(atoms.columns)

    if missing:
        for code, expected in sorted(expected_counts.items()):
            rows.append({
                "RECRUITER_CODE": code,
                "Expected_SMILES_Atom_Count": expected,
                "Observed_Atom_Table_Count": "",
                "Min_Atom_Index": "",
                "Max_Atom_Index": "",
                "Duplicate_Atom_Indices": "",
                "Status": "ATOM_TABLE_BAD_SCHEMA",
                "Notes": f"Missing columns: {sorted(missing)}",
            })
        return rows

    atoms = atoms.copy()
    atoms["RECRUITER_CODE"] = atoms["RECRUITER_CODE"].astype(str).str.strip()
    atoms["smiles_atom_index"] = pd.to_numeric(atoms["smiles_atom_index"], errors="coerce")

    grouped = {code: g.copy() for code, g in atoms.groupby("RECRUITER_CODE", dropna=False)}

    for code, expected in sorted(expected_counts.items()):
        g = grouped.get(code)

        if g is None or len(g) == 0:
            rows.append({
                "RECRUITER_CODE": code,
                "Expected_SMILES_Atom_Count": expected,
                "Observed_Atom_Table_Count": 0,
                "Min_Atom_Index": "",
                "Max_Atom_Index": "",
                "Duplicate_Atom_Indices": "",
                "Status": "MISSING_CODE_IN_ATOM_TABLE",
                "Notes": "",
            })
            continue

        observed = len(g)
        idxs = g["smiles_atom_index"].dropna().astype(int).tolist()
        duplicates = sorted({i for i in idxs if idxs.count(i) > 1})

        if observed == expected and not duplicates:
            status = "PASS"
            notes = ""
        elif observed != expected:
            status = "COUNT_MISMATCH"
            notes = f"expected={expected}, observed={observed}"
        else:
            status = "DUPLICATE_ATOM_INDICES"
            notes = f"duplicate_indices={duplicates}"

        rows.append({
            "RECRUITER_CODE": code,
            "Expected_SMILES_Atom_Count": expected,
            "Observed_Atom_Table_Count": observed,
            "Min_Atom_Index": min(idxs) if idxs else "",
            "Max_Atom_Index": max(idxs) if idxs else "",
            "Duplicate_Atom_Indices": ";".join(map(str, duplicates)) if duplicates else "",
            "Status": status,
            "Notes": notes,
        })

    extra_codes = sorted(set(grouped.keys()) - set(expected_counts.keys()))
    for code in extra_codes:
        g = grouped[code]
        idxs = g["smiles_atom_index"].dropna().astype(int).tolist()
        rows.append({
            "RECRUITER_CODE": code,
            "Expected_SMILES_Atom_Count": "",
            "Observed_Atom_Table_Count": len(g),
            "Min_Atom_Index": min(idxs) if idxs else "",
            "Max_Atom_Index": max(idxs) if idxs else "",
            "Duplicate_Atom_Indices": "",
            "Status": "EXTRA_CODE_IN_ATOM_TABLE",
            "Notes": "Code exists in atom table but not in current Ligase_SMILE_Codes.csv",
        })

    return rows


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate chemical descriptors for LR recruiter SMILES."
    )

    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input Ligase_SMILE_Codes.csv")
    parser.add_argument("--desc-out", default=str(DEFAULT_DESC_OUT), help="Output descriptor CSV")
    parser.add_argument("--audit-out", default=str(DEFAULT_AUDIT_OUT), help="Output audit CSV")
    parser.add_argument("--atoms-out", default=str(DEFAULT_ATOMS_OUT), help="SMILES atom table path")
    parser.add_argument("--atom-validation-out", default=str(DEFAULT_ATOM_VALIDATION_OUT), help="Atom table validation CSV")
    parser.add_argument("--log", default=str(DEFAULT_LOG_FILE), help="Output log file")

    parser.add_argument(
        "--write-atoms",
        action="store_true",
        help="Regenerate Ligase_SMILE_Codes_Atoms.csv from SMILES. Default only validates existing table.",
    )

    parser.add_argument(
        "--allow-non-lr-codes",
        action="store_true",
        help="Allow recruiter codes that do not match LR00001 format.",
    )

    parser.add_argument("--overwrite", action="store_true", help="Backup and overwrite outputs")
    parser.add_argument("--dry-run", action="store_true", help="Preview without writing files")

    args = parser.parse_args()

    input_path = Path(args.input)
    desc_out = Path(args.desc_out)
    audit_out = Path(args.audit_out)
    atoms_out = Path(args.atoms_out)
    atom_validation_out = Path(args.atom_validation_out)
    log_file = Path(args.log)

    require_file(input_path)

    outputs = [desc_out, audit_out, atom_validation_out, log_file]
    if args.write_atoms:
        outputs.append(atoms_out)

    if not args.dry_run:
        existing_outputs = [p for p in outputs if p.exists()]
        if existing_outputs and not args.overwrite:
            print("❌ Output file(s) already exist:")
            for p in existing_outputs:
                print(f"   - {p}")
            print("Run with --overwrite to back them up and regenerate.")
            return 2

        if args.overwrite:
            for p in existing_outputs:
                b = backup_if_exists(p)
                if b:
                    print(f"🧷 Backed up {p} -> {b}")

    print("=== Ligase Chemical Descriptor Generation ===")
    print(f"Started: {now()}")
    print(f"Input: {input_path}")
    print(f"Descriptor output: {desc_out}")
    print(f"Audit output: {audit_out}")
    print(f"Atom table: {atoms_out}")
    print(f"Atom validation output: {atom_validation_out}")
    print(f"Write atom table: {args.write_atoms}")
    print(f"SA scorer available: {HAS_SA}")
    print(f"Dry run: {args.dry_run}")

    df = pd.read_csv(input_path)

    required = {"SMILES", "RECRUITER_CODE"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{input_path} is missing required columns: {missing}")

    df = df.copy()
    df["SMILES"] = df["SMILES"].astype(str).str.strip()
    df["RECRUITER_CODE"] = df["RECRUITER_CODE"].astype(str).str.strip()
    df = df[(df["SMILES"] != "") & (df["RECRUITER_CODE"] != "")]
    df = df.drop_duplicates(subset=["RECRUITER_CODE", "SMILES"]).reset_index(drop=True)

    if not args.allow_non_lr_codes:
        bad_code_mask = ~df["RECRUITER_CODE"].str.match(LR_PATTERN)
        if bad_code_mask.any():
            examples = df.loc[bad_code_mask, "RECRUITER_CODE"].head(10).tolist()
            raise ValueError(
                "Found recruiter codes that do not match LR00001 format. "
                f"Examples: {examples}. Use --allow-non-lr-codes to bypass."
            )

    duplicated_codes = df[df.duplicated("RECRUITER_CODE", keep=False)]
    if len(duplicated_codes):
        raise ValueError(
            "Duplicate RECRUITER_CODE values found in input. "
            "Each LR code should map to one SMILES."
        )

    print(f"\n🔍 Loaded {len(df)} recruiter SMILES from {input_path.name}")
    print(f"🔹 Unique SMILES: {df['SMILES'].nunique()}")
    print(f"🔹 Unique recruiter codes: {df['RECRUITER_CODE'].nunique()}")

    pains_catalog = build_filter_catalog([
        FilterCatalogParams.FilterCatalogs.PAINS_A,
        FilterCatalogParams.FilterCatalogs.PAINS_B,
        FilterCatalogParams.FilterCatalogs.PAINS_C,
    ])

    brenk_catalog = build_filter_catalog([
        FilterCatalogParams.FilterCatalogs.BRENK,
    ])

    atom_rows: List[Dict[str, object]] = []
    desc_rows: List[Dict[str, object]] = []
    audit_rows: List[Dict[str, object]] = []
    expected_atom_counts: Dict[str, int] = {}

    for _, row in df.iterrows():
        code = clean_str(row["RECRUITER_CODE"])
        smi = clean_str(row["SMILES"])

        mol, parse_note = safe_mol(smi)

        if mol is None:
            audit_rows.append({
                "RECRUITER_CODE": code,
                "SMILES": smi,
                "Status": "FAILED",
                "Notes": "",
                "Error": parse_note,
            })
            continue

        expected_atom_counts[code] = mol.GetNumAtoms()

        atom_rows.extend(make_atom_rows(code, mol))

        try:
            desc = compute_descriptors(
                code=code,
                smiles=smi,
                mol=mol,
                pains_catalog=pains_catalog,
                brenk_catalog=brenk_catalog,
            )
            desc_rows.append(desc)

            audit_rows.append({
                "RECRUITER_CODE": code,
                "SMILES": smi,
                "Status": "OK" if not parse_note else "CHECK",
                "Notes": parse_note,
                "Error": "",
            })

        except Exception as e:
            audit_rows.append({
                "RECRUITER_CODE": code,
                "SMILES": smi,
                "Status": "FAILED_DESCRIPTOR_CALC",
                "Notes": parse_note,
                "Error": f"{type(e).__name__}: {e}",
            })

    atom_validation_rows = validate_existing_atom_table(
        atom_table_path=atoms_out,
        expected_counts=expected_atom_counts,
    )

    status_counts = pd.Series([r["Status"] for r in audit_rows]).value_counts(dropna=False).to_dict()
    atom_validation_counts = pd.Series([r["Status"] for r in atom_validation_rows]).value_counts(dropna=False).to_dict()

    print("\n=== Planned / Computed Counts ===")
    print(f"Descriptor rows: {len(desc_rows)}")
    print(f"Audit rows: {len(audit_rows)}")
    print(f"Generated atom rows in memory: {len(atom_rows)}")
    print("Descriptor audit status counts:")
    for k, v in status_counts.items():
        print(f"  {k}: {v}")
    print("Atom table validation counts:")
    for k, v in atom_validation_counts.items():
        print(f"  {k}: {v}")

    if args.dry_run:
        print("\n🧪 Dry run complete. Nothing written.")
        print("Run for real with:")
        print("  python 12_LigaseChemicalDescriptors.py --overwrite")
        return 0

    write_csv(desc_out, DESCRIPTOR_COLUMNS, desc_rows)
    write_csv(audit_out, AUDIT_COLUMNS, audit_rows)
    write_csv(atom_validation_out, ATOM_VALIDATION_COLUMNS, atom_validation_rows)

    if args.write_atoms:
        write_csv(atoms_out, ATOM_COLUMNS, atom_rows)

    log_payload = {
        "started_finished": now(),
        "input": str(input_path),
        "desc_out": str(desc_out),
        "audit_out": str(audit_out),
        "atoms_out": str(atoms_out),
        "atom_validation_out": str(atom_validation_out),
        "write_atoms": args.write_atoms,
        "sa_scorer_available": HAS_SA,
        "counts": {
            "input_rows": len(df),
            "unique_smiles": int(df["SMILES"].nunique()),
            "unique_recruiter_codes": int(df["RECRUITER_CODE"].nunique()),
            "descriptor_rows": len(desc_rows),
            "audit_rows": len(audit_rows),
            "generated_atom_rows": len(atom_rows),
        },
        "descriptor_audit_status_counts": status_counts,
        "atom_validation_status_counts": atom_validation_counts,
    }

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text(json.dumps(log_payload, indent=2) + "\n")

    print(f"\n✅ Written {desc_out}")
    print(f"✅ Written {audit_out}")
    print(f"✅ Written {atom_validation_out}")
    if args.write_atoms:
        print(f"✅ Written {atoms_out}")
    else:
        print(f"ℹ️ Did not overwrite {atoms_out}; used validation-only mode.")
    print(f"🧾 Written {log_file}")

    failed = [r for r in audit_rows if str(r["Status"]).startswith("FAILED")]
    if failed:
        print(f"\n⚠️ {len(failed)} descriptor failures. Review {audit_out}")

    atom_bad = [r for r in atom_validation_rows if r["Status"] != "PASS"]
    if atom_bad:
        print(f"\n⚠️ {len(atom_bad)} atom-table validation issues. Review {atom_validation_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())