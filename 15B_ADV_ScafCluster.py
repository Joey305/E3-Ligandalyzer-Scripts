#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
15B_ADV_ScafCluster.py
===============================================================================
Purpose:
    Build global scaffold "superclusters" across ligases after running:

        14_ligasescaffolds.py
        15_ADV_ScaffoldData.py

Why:
    14_ligasescaffolds.py creates ligase-local scaffold IDs, e.g.:

        CRBN_SCAF_1
        VHL_SCAF_1
        MDM2_SCAF_1

    But the same chemical scaffold can appear across multiple ligases. This script
    groups scaffolds globally by Scaffold_Hash and assigns unified supercluster IDs:

        LR-SCAF00001
        LR-SCAF00002
        LR-SCAF00003

Inputs:
    Ligase_Table/Ligase_Recruiters_Scaffold.csv
        Expected columns:
            Ligase
            Scaffold_ID
            RECRUITER_CODE
            Scaffold_SMILES
            Scaffold_Hash

    Ligase_Table/Ligase_Scaffold_Data.csv
        Created by 15_ADV_ScaffoldData.py.
        Expected useful columns:
            Ligase
            Scaffold_ID
            Recruiter_Count
            Scaffold_SMILES
            Murcko_SMILES
            Scaffold_Class
            Recruiter_Density_Score
            Shannon_Diversity_Index
            Normalized_Diversity

Outputs:
    Ligase_Table/Ligase_Scaffold_Superclusters.csv
        One row per global scaffold supercluster.

    Ligase_Table/Ligase_Recruiters_Superclustered.csv
        One row per recruiter with its ligase-local scaffold and global supercluster.

    Ligase_Table/Ligase_Scaffold_Supercluster_Frequency.csv
        One row per Ligase + Unified_Scaffold_ID.

    Ligase_Table/Ligase_Scaffold_Supercluster_Matrix.csv
        Pivot table: ligases x global superclusters.

    Ligase_Table/Ligase_Scaffold_Supercluster.log

Usage:
    python 15B_ADV_ScafCluster.py --dry-run
    python 15B_ADV_ScafCluster.py --overwrite
===============================================================================
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


# =============================================================================
# Defaults
# =============================================================================

TABLE_DIR = Path("Ligase_Table")

RECRUITER_SCAFFOLD_FILE = TABLE_DIR / "Ligase_Recruiters_Scaffold.csv"
SCAFFOLD_DATA_FILE = TABLE_DIR / "Ligase_Scaffold_Data.csv"

OUT_SUPERCLUSTERS = TABLE_DIR / "Ligase_Scaffold_Superclusters.csv"
OUT_RECRUITERS = TABLE_DIR / "Ligase_Recruiters_Superclustered.csv"
OUT_FREQUENCY = TABLE_DIR / "Ligase_Scaffold_Supercluster_Frequency.csv"
OUT_MATRIX = TABLE_DIR / "Ligase_Scaffold_Supercluster_Matrix.csv"
OUT_LOG = TABLE_DIR / "Ligase_Scaffold_Supercluster.log"


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


def hash_string(text: str, n: int = 12) -> str:
    return hashlib.sha1(str(text).encode()).hexdigest()[:n]


def first_nonempty(values) -> str:
    for v in values:
        s = clean_str(v)
        if s:
            return s
    return ""


def join_unique(values, sep: str = ";") -> str:
    vals = sorted({clean_str(v) for v in values if clean_str(v)})
    return sep.join(vals)


# =============================================================================
# Loading / normalization
# =============================================================================

def load_recruiter_scaffolds(path: Path) -> pd.DataFrame:
    require_file(path)
    df = pd.read_csv(path)

    required = {"Ligase", "Scaffold_ID", "RECRUITER_CODE", "Scaffold_SMILES"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")

    df = df.copy()
    df["Ligase"] = df["Ligase"].astype(str).str.strip()
    df["Scaffold_ID"] = df["Scaffold_ID"].astype(str).str.strip()
    df["RECRUITER_CODE"] = df["RECRUITER_CODE"].astype(str).str.strip()
    df["Scaffold_SMILES"] = df["Scaffold_SMILES"].astype(str).str.strip()

    if "Scaffold_Hash" not in df.columns:
        df["Scaffold_Hash"] = df["Scaffold_SMILES"].apply(lambda x: hash_string(x, n=8) if clean_str(x) else "")
    else:
        df["Scaffold_Hash"] = df["Scaffold_Hash"].astype(str).str.strip()

    # Fallback for empty hash.
    empty_hash = df["Scaffold_Hash"].eq("") | df["Scaffold_Hash"].str.lower().eq("nan")
    df.loc[empty_hash, "Scaffold_Hash"] = df.loc[empty_hash, "Scaffold_SMILES"].apply(
        lambda x: hash_string(x, n=8) if clean_str(x) else ""
    )

    return df


def load_scaffold_data(path: Path) -> pd.DataFrame:
    require_file(path)
    df = pd.read_csv(path)

    required = {"Ligase", "Scaffold_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}")

    df = df.copy()
    df["Ligase"] = df["Ligase"].astype(str).str.strip()
    df["Scaffold_ID"] = df["Scaffold_ID"].astype(str).str.strip()

    if "Recruiter_Count" not in df.columns:
        df["Recruiter_Count"] = 1

    if "Scaffold_SMILES" not in df.columns:
        df["Scaffold_SMILES"] = ""

    if "Murcko_SMILES" not in df.columns:
        df["Murcko_SMILES"] = ""

    if "Scaffold_Class" not in df.columns:
        df["Scaffold_Class"] = "Unknown"

    return df


def merge_inputs(recruiter_df: pd.DataFrame, data_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add advanced scaffold data onto recruiter-level scaffold rows.
    """
    data_cols = [
        c for c in [
            "Ligase",
            "Scaffold_ID",
            "Recruiter_Count",
            "Scaffold_SMILES",
            "Murcko_SMILES",
            "Scaffold_Class",
            "Scaffold_Center_of_Mass_X",
            "Scaffold_Center_of_Mass_Y",
            "Scaffold_Center_of_Mass_Z",
            "Ligase_Scaffold_Connectivity",
            "Total_Recruiters",
            "Recruiter_Density_Score",
            "Shannon_Diversity_Index",
            "Normalized_Diversity",
        ]
        if c in data_df.columns
    ]

    merged = recruiter_df.merge(
        data_df[data_cols],
        on=["Ligase", "Scaffold_ID"],
        how="left",
        suffixes=("", "_ADV"),
    )

    # Prefer advanced Scaffold_SMILES only if original is blank.
    if "Scaffold_SMILES_ADV" in merged.columns:
        blank = merged["Scaffold_SMILES"].astype(str).str.strip().eq("")
        merged.loc[blank, "Scaffold_SMILES"] = merged.loc[blank, "Scaffold_SMILES_ADV"]
        merged = merged.drop(columns=["Scaffold_SMILES_ADV"])

    return merged


# =============================================================================
# Supercluster creation
# =============================================================================

def assign_supercluster_ids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign global LR-SCAF IDs based primarily on Scaffold_Hash.

    Grouping key:
        Scaffold_Hash

    Fallback if Scaffold_Hash is missing:
        hash(Murcko_SMILES or Scaffold_SMILES or Scaffold_ID)
    """
    df = df.copy()

    def make_key(row):
        scaffold_hash = clean_str(row.get("Scaffold_Hash"))
        if scaffold_hash and scaffold_hash.lower() != "nan":
            return scaffold_hash

        murcko = clean_str(row.get("Murcko_SMILES"))
        if murcko:
            return "MURCKO_" + hash_string(murcko, n=12)

        smi = clean_str(row.get("Scaffold_SMILES"))
        if smi:
            return "SMI_" + hash_string(smi, n=12)

        return "SCAFID_" + hash_string(clean_str(row.get("Scaffold_ID")), n=12)

    df["Supercluster_Key"] = df.apply(make_key, axis=1)

    unique_keys = (
        df[["Supercluster_Key"]]
        .drop_duplicates()
        .sort_values("Supercluster_Key")
        .reset_index(drop=True)
    )

    unique_keys["Unified_Scaffold_ID"] = [
        f"LR-SCAF{i:05d}" for i in range(1, len(unique_keys) + 1)
    ]

    key_map = dict(zip(unique_keys["Supercluster_Key"], unique_keys["Unified_Scaffold_ID"]))
    df["Unified_Scaffold_ID"] = df["Supercluster_Key"].map(key_map)

    return df


def build_supercluster_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for key, group in df.groupby("Supercluster_Key", dropna=False):
        unified_id = first_nonempty(group["Unified_Scaffold_ID"])
        scaffold_hash = first_nonempty(group["Scaffold_Hash"])
        scaffold_smiles = first_nonempty(group["Scaffold_SMILES"])
        murcko_smiles = first_nonempty(group["Murcko_SMILES"]) if "Murcko_SMILES" in group.columns else ""
        scaffold_class = first_nonempty(group["Scaffold_Class"]) if "Scaffold_Class" in group.columns else "Unknown"

        ligases = sorted(group["Ligase"].dropna().astype(str).unique().tolist())
        scaffold_ids = sorted(group["Scaffold_ID"].dropna().astype(str).unique().tolist())
        recruiters = sorted(group["RECRUITER_CODE"].dropna().astype(str).unique().tolist())

        # If Recruiter_Count exists at scaffold-level, do not sum it per recruiter row,
        # because recruiter-level rows would duplicate it. Count unique recruiters instead.
        recruiter_count = len(recruiters)

        rows.append({
            "Unified_Scaffold_ID": unified_id,
            "Supercluster_Key": key,
            "Scaffold_Hash": scaffold_hash,
            "Representative_Scaffold_SMILES": scaffold_smiles,
            "Representative_Murcko_SMILES": murcko_smiles,
            "Scaffold_Class": scaffold_class,
            "Ligase_Count": len(ligases),
            "Ligases": ";".join(ligases),
            "Ligase_Local_Scaffold_Count": len(scaffold_ids),
            "Original_Scaffold_ID_Group": ";".join(scaffold_ids),
            "Recruiter_Count": recruiter_count,
            "Recruiter_Code_Group": ";".join(recruiters),
            "Is_Cross_Ligase_Supercluster": len(ligases) > 1,
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["Ligase_Count", "Recruiter_Count", "Unified_Scaffold_ID"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    return out


def build_frequency(df: pd.DataFrame) -> pd.DataFrame:
    freq = (
        df.groupby(["Ligase", "Unified_Scaffold_ID"], dropna=False)
        .agg(
            Recruiter_Count=("RECRUITER_CODE", "nunique"),
            Ligase_Local_Scaffold_Count=("Scaffold_ID", "nunique"),
            Scaffold_ID_Group=("Scaffold_ID", lambda x: ";".join(sorted(set(map(str, x))))),
            Recruiter_Code_Group=("RECRUITER_CODE", lambda x: ";".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )

    totals = freq.groupby("Ligase")["Recruiter_Count"].sum().to_dict()
    freq["Total_Recruiters_For_Ligase"] = freq["Ligase"].map(totals)
    freq["Supercluster_Density_For_Ligase"] = (
        freq["Recruiter_Count"] / freq["Total_Recruiters_For_Ligase"]
    )

    freq = freq.sort_values(
        ["Ligase", "Recruiter_Count", "Unified_Scaffold_ID"],
        ascending=[True, False, True],
    ).reset_index(drop=True)

    return freq


def build_matrix(freq: pd.DataFrame) -> pd.DataFrame:
    mat = freq.pivot_table(
        index="Ligase",
        columns="Unified_Scaffold_ID",
        values="Recruiter_Count",
        aggfunc="sum",
        fill_value=0,
    )

    mat = mat.reset_index()
    mat.columns.name = None
    return mat


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create global LR-SCAF superclusters from ligase-local scaffold hashes."
    )

    parser.add_argument("--table-dir", default=str(TABLE_DIR))
    parser.add_argument("--recruiter-scaffolds", default=str(RECRUITER_SCAFFOLD_FILE))
    parser.add_argument("--scaffold-data", default=str(SCAFFOLD_DATA_FILE))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")

    args = parser.parse_args()

    table_dir = Path(args.table_dir)
    recruiter_scaffold_file = Path(args.recruiter_scaffolds)
    scaffold_data_file = Path(args.scaffold_data)

    out_superclusters = table_dir / OUT_SUPERCLUSTERS.name
    out_recruiters = table_dir / OUT_RECRUITERS.name
    out_frequency = table_dir / OUT_FREQUENCY.name
    out_matrix = table_dir / OUT_MATRIX.name
    out_log = table_dir / OUT_LOG.name

    outputs = [
        out_superclusters,
        out_recruiters,
        out_frequency,
        out_matrix,
        out_log,
    ]

    if not args.dry_run:
        existing = [p for p in outputs if p.exists()]
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

    print("=== 15B Advanced Scaffold Superclustering ===")
    print(f"Started: {now()}")
    print(f"Recruiter scaffold input: {recruiter_scaffold_file}")
    print(f"Advanced scaffold data:   {scaffold_data_file}")
    print(f"Output dir:               {table_dir}")
    print(f"Dry run:                  {args.dry_run}")

    recruiter_df = load_recruiter_scaffolds(recruiter_scaffold_file)
    data_df = load_scaffold_data(scaffold_data_file)

    merged = merge_inputs(recruiter_df, data_df)
    clustered = assign_supercluster_ids(merged)

    superclusters = build_supercluster_summary(clustered)
    frequency = build_frequency(clustered)
    matrix = build_matrix(frequency)

    print("\n=== Counts ===")
    print(f"Recruiter scaffold rows:       {len(recruiter_df)}")
    print(f"Unique recruiters:             {recruiter_df['RECRUITER_CODE'].nunique()}")
    print(f"Ligase-local scaffold IDs:     {recruiter_df['Scaffold_ID'].nunique()}")
    print(f"Unique scaffold hashes:        {clustered['Supercluster_Key'].nunique()}")
    print(f"Global superclusters:          {len(superclusters)}")
    print(f"Cross-ligase superclusters:    {int(superclusters['Is_Cross_Ligase_Supercluster'].sum())}")
    print(f"Frequency rows:                {len(frequency)}")
    print(f"Matrix shape:                  {matrix.shape[0]} ligases x {matrix.shape[1] - 1} superclusters")

    if args.dry_run:
        print("\n🧪 Dry run complete. Nothing written.")
        print("Run for real with:")
        print("  python 15B_ADV_ScafCluster.py --overwrite")
        return 0

    table_dir.mkdir(parents=True, exist_ok=True)

    superclusters.to_csv(out_superclusters, index=False)
    clustered.to_csv(out_recruiters, index=False)
    frequency.to_csv(out_frequency, index=False)
    matrix.to_csv(out_matrix, index=False)

    log_payload = {
        "started_finished": now(),
        "inputs": {
            "recruiter_scaffolds": str(recruiter_scaffold_file),
            "scaffold_data": str(scaffold_data_file),
        },
        "outputs": {
            "superclusters": str(out_superclusters),
            "recruiters_superclustered": str(out_recruiters),
            "frequency": str(out_frequency),
            "matrix": str(out_matrix),
            "log": str(out_log),
        },
        "counts": {
            "recruiter_scaffold_rows": len(recruiter_df),
            "unique_recruiters": int(recruiter_df["RECRUITER_CODE"].nunique()),
            "ligase_local_scaffold_ids": int(recruiter_df["Scaffold_ID"].nunique()),
            "global_superclusters": len(superclusters),
            "cross_ligase_superclusters": int(superclusters["Is_Cross_Ligase_Supercluster"].sum()),
            "frequency_rows": len(frequency),
            "matrix_ligases": int(matrix.shape[0]),
            "matrix_superclusters": int(matrix.shape[1] - 1),
        },
    }

    out_log.write_text(json.dumps(log_payload, indent=2) + "\n")

    print(f"\n✅ Written global superclusters → {out_superclusters}")
    print(f"✅ Written recruiter-level supercluster map → {out_recruiters}")
    print(f"✅ Written ligase-supercluster frequency table → {out_frequency}")
    print(f"✅ Written ligase x supercluster matrix → {out_matrix}")
    print(f"🧾 Written log → {out_log}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
