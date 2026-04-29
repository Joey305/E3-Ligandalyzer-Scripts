#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
16_TBLCreation.py
===============================================================================
Purpose:
    Generate database-ready CSV tables that connect legacy ligand names to the
    new LR recruiter codes and their SMILES atom indices.

Creates:
    1. Ligase_Table/Ligase_Duplicate_Ligands.csv

       Columns:
         - Ligand
         - MATCHED_RECRUITERS

       Meaning:
         One row per ligand code that appears in more than one recruiter instance.
         MATCHED_RECRUITERS contains comma-separated LR recruiter codes.

    2. Ligase_Table/Ligase_Ligands_Smiles.csv

       Columns:
         - Ligand
         - RECRUITER_CODE
         - smiles_atom_index
         - smile_atom

       Meaning:
         Atom-level SMILES atom table linked back to the original ligand code and
         the new LR recruiter code.

Primary inputs:
    Ligase_Table/Ligand_Instance_Recruiter_Codes.csv
    Ligase_Table/Ligase_SMILE_Codes_Atoms.csv

Fallback inputs:
    Ligase_Table/Ligase_Ligand_Metadata.csv
    Ligase_Table/Ligase_SMILE_Codes.csv

Usage:
    python 16_TBLCreation.py --dry-run
    python 16_TBLCreation.py --overwrite
===============================================================================
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd


# =============================================================================
# Defaults
# =============================================================================

TABLE_DIR = Path("Ligase_Table")

INSTANCE_FILE = TABLE_DIR / "Ligand_Instance_Recruiter_Codes.csv"
METADATA_FILE = TABLE_DIR / "Ligase_Ligand_Metadata.csv"
SMILES_ATOMS_FILE = TABLE_DIR / "Ligase_SMILE_Codes_Atoms.csv"

OUT_DUPLICATES = TABLE_DIR / "Ligase_Duplicate_Ligands.csv"
OUT_LIGAND_SMILES = TABLE_DIR / "Ligase_Ligands_Smiles.csv"
OUT_AUDIT = TABLE_DIR / "16_TBLCreation_Audit.csv"
OUT_LOG = TABLE_DIR / "16_TBLCreation.log"
OUT_MASTER_MAP = TABLE_DIR / "Recruiter_Master_Map.csv"
LR_PATTERN = re.compile(r"^LR\d{5}$")


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


def lr_sort_key(code: str):
    """
    Sort LR00001, LR00002, ..., LR01000 naturally.
    """
    code = clean_str(code)
    m = re.match(r"^LR(\d+)$", code)
    if m:
        return int(m.group(1))
    return code


# =============================================================================
# Load source mapping
# =============================================================================

def load_recruiter_ligand_map(instance_file: Path, metadata_file: Path) -> pd.DataFrame:
    """
    Preferred source:
        Ligand_Instance_Recruiter_Codes.csv

    Fallback:
        Ligase_Ligand_Metadata.csv
    """
    if instance_file.exists():
        df = pd.read_csv(instance_file)
        source = str(instance_file)
    else:
        require_file(metadata_file)
        df = pd.read_csv(metadata_file)
        source = str(metadata_file)

    required = {"RECRUITER_CODE", "Ligand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{source} missing required columns: {missing}")

    keep_cols = [c for c in ["Ligand", "RECRUITER_CODE", "Ligase", "pdb_id", "Variant"] if c in df.columns]
    df = df[keep_cols].copy()

    df["Ligand"] = df["Ligand"].astype(str).str.strip().str.upper()
    df["RECRUITER_CODE"] = df["RECRUITER_CODE"].astype(str).str.strip()

    df = df[(df["Ligand"] != "") & (df["RECRUITER_CODE"] != "")]
    df = df.drop_duplicates(subset=["Ligand", "RECRUITER_CODE"]).reset_index(drop=True)

    return df


def load_smiles_atoms(path: Path) -> pd.DataFrame:
    require_file(path)

    df = pd.read_csv(path)

    required = {"RECRUITER_CODE", "smiles_atom_index", "smile_atom"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")

    df = df.copy()
    df["RECRUITER_CODE"] = df["RECRUITER_CODE"].astype(str).str.strip()
    df["smiles_atom_index"] = pd.to_numeric(df["smiles_atom_index"], errors="coerce").astype("Int64")
    df["smile_atom"] = df["smile_atom"].astype(str).str.strip()

    df = df[(df["RECRUITER_CODE"] != "") & df["smiles_atom_index"].notna()]
    df = df.drop_duplicates(subset=["RECRUITER_CODE", "smiles_atom_index", "smile_atom"]).reset_index(drop=True)

    return df


# =============================================================================
# Build outputs
# =============================================================================

def build_duplicate_ligands(recruiter_map: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for ligand, group in recruiter_map.groupby("Ligand", dropna=False):
        recruiters = sorted(
            set(group["RECRUITER_CODE"].dropna().astype(str)),
            key=lr_sort_key,
        )

        # This table is specifically duplicate/shared ligand names.
        if len(recruiters) <= 1:
            continue

        rows.append({
            "Ligand": ligand,
            "MATCHED_RECRUITERS": ",".join(recruiters),
        })

    out = pd.DataFrame(rows, columns=["Ligand", "MATCHED_RECRUITERS"])
    out = out.sort_values("Ligand").reset_index(drop=True)
    return out


def build_ligand_smiles_table(recruiter_map: pd.DataFrame, smiles_atoms: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Join:
        RECRUITER_CODE -> Ligand
    onto:
        RECRUITER_CODE -> smiles_atom_index, smile_atom
    """
    map_small = recruiter_map[["Ligand", "RECRUITER_CODE"]].drop_duplicates()

    merged = smiles_atoms.merge(
        map_small,
        on="RECRUITER_CODE",
        how="left",
    )

    audit_rows = []

    missing_ligand = merged[merged["Ligand"].isna()].copy()
    for code in sorted(missing_ligand["RECRUITER_CODE"].dropna().unique(), key=lr_sort_key):
        audit_rows.append({
            "Issue": "SMILES_ATOM_CODE_MISSING_LIGAND_MAP",
            "RECRUITER_CODE": code,
            "Ligand": "",
            "Rows": len(missing_ligand[missing_ligand["RECRUITER_CODE"] == code]),
            "Notes": "RECRUITER_CODE exists in Ligase_SMILE_Codes_Atoms.csv but not in ligand-instance map.",
        })

    out = merged[merged["Ligand"].notna()].copy()
    out = out[["Ligand", "RECRUITER_CODE", "smiles_atom_index", "smile_atom"]]

    out = out.sort_values(
        by=["Ligand", "RECRUITER_CODE", "smiles_atom_index"],
        key=lambda col: col.map(lr_sort_key) if col.name == "RECRUITER_CODE" else col,
    ).reset_index(drop=True)

    audit_df = pd.DataFrame(
        audit_rows,
        columns=["Issue", "RECRUITER_CODE", "Ligand", "Rows", "Notes"],
    )

    return out, audit_df

def build_recruiter_master_map(recruiter_map: pd.DataFrame) -> pd.DataFrame:
    """
    Build Recruiter_Master_Map.csv.

    One row per RECRUITER_CODE.
    For each Ligand, the first LR code becomes the MasterRecruiter.
    If a Ligand has only one recruiter, that recruiter is its own master.
    """
    rows = []

    for ligand, group in recruiter_map.groupby("Ligand", dropna=False):
        recruiters = sorted(
            set(group["RECRUITER_CODE"].dropna().astype(str)),
            key=lr_sort_key,
        )

        if not recruiters:
            continue

        master = recruiters[0]

        for code in recruiters:
            rows.append({
                "RECRUITER_CODE": code,
                "MasterRecruiter": master,
                "Ligand": ligand,
            })

    out = pd.DataFrame(
        rows,
        columns=["RECRUITER_CODE", "MasterRecruiter", "Ligand"]
    )

    out = out.sort_values(
        by=["Ligand", "RECRUITER_CODE"],
        key=lambda col: col.map(lr_sort_key) if col.name in {"RECRUITER_CODE", "MasterRecruiter"} else col,
    ).reset_index(drop=True)

    return out
# =============================================================================
# Main
# =============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create duplicate ligand and ligand-SMILES atom CSV tables."
    )

    parser.add_argument("--table-dir", default=str(TABLE_DIR))
    parser.add_argument("--instances", default=str(INSTANCE_FILE))
    parser.add_argument("--metadata", default=str(METADATA_FILE))
    parser.add_argument("--smiles-atoms", default=str(SMILES_ATOMS_FILE))
    parser.add_argument("--duplicates-out", default=str(OUT_DUPLICATES))
    parser.add_argument("--ligand-smiles-out", default=str(OUT_LIGAND_SMILES))
    parser.add_argument("--audit-out", default=str(OUT_AUDIT))
    parser.add_argument("--log", default=str(OUT_LOG))
    parser.add_argument("--allow-non-lr-codes", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")

    args = parser.parse_args()

    table_dir = Path(args.table_dir)
    instance_file = Path(args.instances)
    metadata_file = Path(args.metadata)
    smiles_atoms_file = Path(args.smiles_atoms)
    duplicates_out = Path(args.duplicates_out)
    ligand_smiles_out = Path(args.ligand_smiles_out)
    audit_out = Path(args.audit_out)
    log_file = Path(args.log)

    master_map_out = OUT_MASTER_MAP

    outputs = [duplicates_out, ligand_smiles_out, master_map_out, audit_out, log_file]

    if not args.dry_run:
        existing = [p for p in outputs if p.exists()]
        if existing and not args.overwrite:
            print("❌ Output file(s) already exist:")
            for p in existing:
                print(f"   - {p}")
            print("Run with --overwrite to back them up and regenerate.")
            return 2

        if args.overwrite:
            for p in existing:
                b = backup_if_exists(p)
                if b:
                    print(f"🧷 Backed up {p} -> {b}")

    print("=== 16 Table Creation ===")
    print(f"Started: {now()}")
    print(f"Instance map: {instance_file if instance_file.exists() else '(missing; using metadata fallback)'}")
    print(f"Metadata fallback: {metadata_file}")
    print(f"SMILES atom table: {smiles_atoms_file}")
    print(f"Dry run: {args.dry_run}")

    recruiter_map = load_recruiter_ligand_map(instance_file, metadata_file)
    smiles_atoms = load_smiles_atoms(smiles_atoms_file)

    if not args.allow_non_lr_codes:
        bad_map = recruiter_map[~recruiter_map["RECRUITER_CODE"].str.match(LR_PATTERN)]
        bad_atoms = smiles_atoms[~smiles_atoms["RECRUITER_CODE"].str.match(LR_PATTERN)]

        if len(bad_map) or len(bad_atoms):
            examples = []
            if len(bad_map):
                examples.extend(bad_map["RECRUITER_CODE"].head(5).tolist())
            if len(bad_atoms):
                examples.extend(bad_atoms["RECRUITER_CODE"].head(5).tolist())
            raise ValueError(
                "Found non-LR recruiter codes. Examples: "
                + ", ".join(sorted(set(examples)))
                + ". Use --allow-non-lr-codes to bypass."
            )

    duplicates = build_duplicate_ligands(recruiter_map)
    master_map = build_recruiter_master_map(recruiter_map)
    ligand_smiles, audit = build_ligand_smiles_table(recruiter_map, smiles_atoms)

    print("\n=== Counts ===")
    print(f"Recruiter ↔ ligand map rows: {len(recruiter_map)}")
    print(f"Unique recruiters:           {recruiter_map['RECRUITER_CODE'].nunique()}")
    print(f"Unique ligand names:         {recruiter_map['Ligand'].nunique()}")
    print(f"Duplicate ligand rows:       {len(duplicates)}")
    print(f"Master map rows:             {len(master_map)}")
    print(f"SMILES atom input rows:      {len(smiles_atoms)}")
    print(f"Ligand-SMILES output rows:   {len(ligand_smiles)}")
    print(f"Audit rows:                  {len(audit)}")

    if args.dry_run:
        print("\n🧪 Dry run complete. Nothing written.")
        print("Run for real with:")
        print("  python 16_TBLCreation.py --overwrite")
        return 0

    table_dir.mkdir(parents=True, exist_ok=True)

    duplicates.to_csv(duplicates_out, index=False)
    ligand_smiles.to_csv(ligand_smiles_out, index=False)
    master_map.to_csv(master_map_out, index=False)
    audit.to_csv(audit_out, index=False)

    log_payload = {
        "started_finished": now(),
        "inputs": {
            "instance_file": str(instance_file),
            "metadata_file": str(metadata_file),
            "smiles_atoms_file": str(smiles_atoms_file),
        },
        "outputs": {
            "ligase_duplicate_ligands": str(duplicates_out),
            "ligase_ligands_smiles": str(ligand_smiles_out),
            "audit": str(audit_out),
            "log": str(log_file),
            "recruiter_master_map": str(master_map_out),
        },
        "counts": {
            "recruiter_ligand_map_rows": len(recruiter_map),
            "unique_recruiters": int(recruiter_map["RECRUITER_CODE"].nunique()),
            "unique_ligand_names": int(recruiter_map["Ligand"].nunique()),
            "duplicate_ligand_rows": len(duplicates),
            "smiles_atom_input_rows": len(smiles_atoms),
            "ligand_smiles_output_rows": len(ligand_smiles),
            "audit_rows": len(audit),
            "master_map_rows": len(master_map),
        },
    }

    log_file.write_text(pd.Series(log_payload).to_json(indent=2))

    print(f"\n✅ Written duplicate ligand table → {duplicates_out}")
    print(f"✅ Written ligand-SMILES atom table → {ligand_smiles_out}")
    print(f"✅ Written recruiter master map → {master_map_out}")
    print(f"✅ Written audit → {audit_out}")
    print(f"🧾 Written log → {log_file}")

    if len(audit):
        print("\n⚠️ Audit rows were written. Review:")
        print(f"  {audit_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
