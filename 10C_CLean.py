#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
10B_CLean.py
===============================================================================
Purpose:
    Propagate the new LR00001-style recruiter codes from MCS mapping outputs
    back into the final Ligase_Table SASA/metadata/SMILES tables.

Why:
    Earlier files still contain old recruiter codes like:
        CHIP_A90
        CRBN_6EL
        MDM2_7HC

    But after 10_MCS_Matching.py, the database-facing recruiter code should be:
        LR00001
        LR00002
        LR00003
        ...

This script updates:
    Ligase_Table/Ligase_Ligand_SASA_summary.csv
    Ligase_Table/Ligase_Ligand_SASA_atoms.csv
    Ligase_Table/Ligase_Ligand_Metadata.csv

And rebuilds:
    Ligase_Table/Recruiter_SMILES_Map.csv
    Ligase_Table/Ligase_SMILE_Codes.csv
    Ligase_Table/Recruiter_SMILES_Wide.csv

Using:
    Ligase_Table/Ligand_Instance_Recruiter_Codes.csv

Core mapping:
    Ligase + pdb_id + Ligand + Variant  ->  LR00001

Metadata expansion:
    Ligand metadata is usually Ligase + Ligand level, while LR codes are
    Ligase + pdb_id + Ligand + Variant level. Therefore, metadata rows are
    expanded to one row per LR instance.

Outputs:
    Updated in-place under Ligase_Table/:
      - Ligase_Ligand_SASA_summary.csv
      - Ligase_Ligand_SASA_atoms.csv
      - Ligase_Ligand_Metadata.csv
      - Recruiter_SMILES_Map.csv
      - Ligase_SMILE_Codes.csv
      - Recruiter_SMILES_Wide.csv
      - Recruiter_Code_Crosswalk.csv
      - 10B_CLean_Audit.csv
      - 10B_CLean.log

Usage:
    python 10B_CLean.py --dry-run
    python 10B_CLean.py --overwrite
===============================================================================
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


# =============================================================================
# Defaults
# =============================================================================

DEFAULT_TABLE_DIR = Path("Ligase_Table")

DEFAULT_INSTANCE_FILE = DEFAULT_TABLE_DIR / "Ligand_Instance_Recruiter_Codes.csv"
DEFAULT_SUMMARY_FILE = DEFAULT_TABLE_DIR / "Ligase_Ligand_SASA_summary.csv"
DEFAULT_ATOMS_FILE = DEFAULT_TABLE_DIR / "Ligase_Ligand_SASA_atoms.csv"
DEFAULT_METADATA_FILE = DEFAULT_TABLE_DIR / "Ligase_Ligand_Metadata.csv"

OUT_RECRUITER_MAP = DEFAULT_TABLE_DIR / "Recruiter_SMILES_Map.csv"
OUT_SMILE_CODES = DEFAULT_TABLE_DIR / "Ligase_SMILE_Codes.csv"
OUT_SMILES_WIDE = DEFAULT_TABLE_DIR / "Recruiter_SMILES_Wide.csv"
OUT_CROSSWALK = DEFAULT_TABLE_DIR / "Recruiter_Code_Crosswalk.csv"
OUT_AUDIT = DEFAULT_TABLE_DIR / "10B_CLean_Audit.csv"
OUT_LOG = DEFAULT_TABLE_DIR / "10B_CLean.log"

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
    backup = path.with_name(f"{path.stem}.backup_{stamp()}{path.suffix}")
    shutil.copy2(path, backup)
    return backup


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")


def make_instance_key(ligase, pdb_id, ligand, variant) -> str:
    return (
        clean_str(ligase) + "|" +
        clean_str(pdb_id).upper() + "|" +
        clean_str(ligand).upper() + "|" +
        normalize_variant(variant)
    )


def make_ligase_ligand_key(ligase, ligand) -> str:
    return clean_str(ligase) + "|" + clean_str(ligand).upper()


def normalize_identity_cols(df: pd.DataFrame, require_pdb: bool = False) -> pd.DataFrame:
    df = df.copy()

    required = {"Ligase", "Ligand"}
    if require_pdb:
        required.add("pdb_id")

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required identity columns: {missing}")

    df["Ligase"] = df["Ligase"].astype(str).str.strip()
    df["Ligand"] = df["Ligand"].astype(str).str.strip().str.upper()

    if "pdb_id" in df.columns:
        df["pdb_id"] = df["pdb_id"].astype(str).str.strip().str.upper()

    if "Variant" not in df.columns:
        df["Variant"] = "1"
    else:
        df["Variant"] = df["Variant"].apply(normalize_variant)

    df["_LL_KEY"] = df.apply(
        lambda r: make_ligase_ligand_key(r["Ligase"], r["Ligand"]),
        axis=1,
    )

    if "pdb_id" in df.columns:
        df["_INSTANCE_KEY"] = df.apply(
            lambda r: make_instance_key(r["Ligase"], r["pdb_id"], r["Ligand"], r["Variant"]),
            axis=1,
        )

    return df


def drop_internal_cols(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(
        columns=[c for c in df.columns if c.startswith("_")],
        errors="ignore",
    )


# =============================================================================
# Load instance LR code map
# =============================================================================

def load_instances(path: Path) -> pd.DataFrame:
    require_file(path)

    df = pd.read_csv(path)
    required = {
        "RECRUITER_CODE",
        "Ligase",
        "pdb_id",
        "Ligand",
        "Variant",
        "SMILES",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")

    df = normalize_identity_cols(df, require_pdb=True)

    df["RECRUITER_CODE"] = df["RECRUITER_CODE"].astype(str).str.strip()
    bad_codes = df[~df["RECRUITER_CODE"].str.match(LR_PATTERN)]
    if len(bad_codes):
        examples = ", ".join(bad_codes["RECRUITER_CODE"].head(10).tolist())
        raise ValueError(
            f"Found RECRUITER_CODE values that do not match LR00001 format. Examples: {examples}"
        )

    if "Canonical_SMILES" not in df.columns:
        df["Canonical_SMILES"] = df["SMILES"]

    if "Original_Ligase_Ligand_Code" not in df.columns:
        df["Original_Ligase_Ligand_Code"] = df["Ligase"] + "_" + df["Ligand"]

    if "Instance_Key" not in df.columns:
        df["Instance_Key"] = df["_INSTANCE_KEY"]

    # Ensure one LR code per instance key.
    dup_instance = df[df.duplicated("_INSTANCE_KEY", keep=False)]
    if len(dup_instance):
        raise ValueError(
            "Duplicate instance keys found in Ligand_Instance_Recruiter_Codes.csv. "
            "Check 10_MCS_Matching output."
        )

    dup_codes = df[df.duplicated("RECRUITER_CODE", keep=False)]
    if len(dup_codes):
        raise ValueError(
            "Duplicate LR recruiter codes found in Ligand_Instance_Recruiter_Codes.csv. "
            "Check 10_MCS_Matching output."
        )

    return df


# =============================================================================
# Update SASA summary / atoms
# =============================================================================

def update_instance_table(
    df: pd.DataFrame,
    instances: pd.DataFrame,
    table_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Update RECRUITER_CODE in a row-level table with Ligase+pdb_id+Ligand+Variant.
    """
    df = normalize_identity_cols(df, require_pdb=True)

    map_df = instances[[
        "_INSTANCE_KEY",
        "RECRUITER_CODE",
        "SMILES",
        "Canonical_SMILES",
        "Original_Ligase_Ligand_Code",
        "Instance_Key",
    ]].copy()

    merged = df.merge(
        map_df,
        on="_INSTANCE_KEY",
        how="left",
        suffixes=("", "_NEW"),
    )

    audit_rows = []

    missing = merged[merged["RECRUITER_CODE_NEW"].isna()].copy()
    for _, r in missing.iterrows():
        audit_rows.append({
            "Table": table_name,
            "Issue": "MISSING_LR_CODE_FOR_INSTANCE",
            "Ligase": r.get("Ligase", ""),
            "pdb_id": r.get("pdb_id", ""),
            "Ligand": r.get("Ligand", ""),
            "Variant": r.get("Variant", ""),
            "Old_RECRUITER_CODE": r.get("RECRUITER_CODE", ""),
            "New_RECRUITER_CODE": "",
            "Notes": r.get("_INSTANCE_KEY", ""),
        })

    # Preserve original code for traceability.
    if "RECRUITER_CODE" in merged.columns:
        merged["Original_RECRUITER_CODE"] = merged["RECRUITER_CODE"]
    else:
        merged["Original_RECRUITER_CODE"] = merged["Ligase"] + "_" + merged["Ligand"]

    merged["RECRUITER_CODE"] = merged["RECRUITER_CODE_NEW"]

    # Drop rows that could not be mapped, because they should not enter DB.
    cleaned = merged[merged["RECRUITER_CODE"].notna()].copy()

    # Add useful instance identity columns.
    cleaned["Original_Ligase_Ligand_Code"] = cleaned["Original_Ligase_Ligand_Code"]
    cleaned["Instance_Key"] = cleaned["Instance_Key"]

    # Remove merge helper columns.
    cleaned = cleaned.drop(
        columns=[
            "RECRUITER_CODE_NEW",
            "SMILES",
            "Canonical_SMILES",
        ],
        errors="ignore",
    )

    audit_df = pd.DataFrame(audit_rows)
    return drop_internal_cols(cleaned), audit_df


# =============================================================================
# Expand metadata to LR instance-level metadata
# =============================================================================

def infer_instance_sdf_path(row: pd.Series) -> str:
    ligase = clean_str(row.get("Ligase"))
    pdb_id = clean_str(row.get("pdb_id")).upper()
    ligand = clean_str(row.get("Ligand")).upper()
    variant = normalize_variant(row.get("Variant"))

    if not ligase or not pdb_id or not ligand:
        return ""

    if variant and variant != "1":
        p = Path("Ligases") / ligase / "SDF_4Download" / f"{pdb_id}_{ligand}_{variant}.sdf"
    else:
        p = Path("Ligases") / ligase / "SDF_4Download" / f"{pdb_id}_{ligand}.sdf"

    return str(p)


def expand_metadata_to_instances(
    meta_df: pd.DataFrame,
    instances: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Metadata is normally one row per Ligase+Ligand.
    LR codes are one row per Ligase+pdb_id+Ligand+Variant.

    This duplicates metadata rows to instance-level LR rows.
    """
    meta_df = normalize_identity_cols(meta_df, require_pdb=False)

    inst_cols = [
        "RECRUITER_CODE",
        "Ligase",
        "pdb_id",
        "Ligand",
        "Variant",
        "SMILES",
        "Canonical_SMILES",
        "Original_Ligase_Ligand_Code",
        "Instance_Key",
        "PDB_File",
        "_LL_KEY",
        "_INSTANCE_KEY",
    ]

    available_inst_cols = [c for c in inst_cols if c in instances.columns]
    inst = instances[available_inst_cols].copy()

    expanded_rows: List[dict] = []
    audit_rows: List[dict] = []

    meta_by_ll = {k: g.copy() for k, g in meta_df.groupby("_LL_KEY", dropna=False)}

    for _, inst_row in inst.iterrows():
        ll_key = inst_row["_LL_KEY"]
        matches = meta_by_ll.get(ll_key)

        if matches is None or len(matches) == 0:
            audit_rows.append({
                "Table": "Ligase_Ligand_Metadata",
                "Issue": "NO_METADATA_FOR_LR_INSTANCE",
                "Ligase": inst_row.get("Ligase", ""),
                "pdb_id": inst_row.get("pdb_id", ""),
                "Ligand": inst_row.get("Ligand", ""),
                "Variant": inst_row.get("Variant", ""),
                "Old_RECRUITER_CODE": inst_row.get("Original_Ligase_Ligand_Code", ""),
                "New_RECRUITER_CODE": inst_row.get("RECRUITER_CODE", ""),
                "Notes": ll_key,
            })
            continue

        # Usually there is one metadata row for Ligase+Ligand.
        # If there are duplicates, take the first and record a warning.
        base = matches.iloc[0].copy()

        if len(matches) > 1:
            audit_rows.append({
                "Table": "Ligase_Ligand_Metadata",
                "Issue": "MULTIPLE_METADATA_ROWS_FOR_LIGASE_LIGAND_USED_FIRST",
                "Ligase": inst_row.get("Ligase", ""),
                "pdb_id": inst_row.get("pdb_id", ""),
                "Ligand": inst_row.get("Ligand", ""),
                "Variant": inst_row.get("Variant", ""),
                "Old_RECRUITER_CODE": base.get("RECRUITER_CODE", ""),
                "New_RECRUITER_CODE": inst_row.get("RECRUITER_CODE", ""),
                "Notes": f"{ll_key}; rows={len(matches)}",
            })

        old_code = base.get("RECRUITER_CODE", "")

        # Convert to dict and overwrite identity fields with instance-level values.
        out = base.to_dict()

        out["Original_RECRUITER_CODE"] = old_code
        out["RECRUITER_CODE"] = inst_row["RECRUITER_CODE"]
        out["Ligase"] = inst_row["Ligase"]
        out["pdb_id"] = inst_row["pdb_id"]
        out["Ligand"] = inst_row["Ligand"]
        out["Variant"] = normalize_variant(inst_row["Variant"])
        out["Original_Ligase_Ligand_Code"] = inst_row.get("Original_Ligase_Ligand_Code", "")
        out["Instance_Key"] = inst_row.get("Instance_Key", inst_row.get("_INSTANCE_KEY", ""))
        out["PDB_File"] = inst_row.get("PDB_File", "")

        # Force SMILES to match the LR code source from 10_MCS.
        # This is important when one ligand code existed in multiple instances.
        out["SMILES"] = inst_row.get("SMILES", out.get("SMILES", ""))
        out["Canonical_SMILES"] = inst_row.get("Canonical_SMILES", out.get("Canonical_SMILES", out.get("SMILES", "")))

        # Point Source_SDF at the per-PDB SDF if it exists.
        instance_sdf = infer_instance_sdf_path(inst_row)
        if instance_sdf:
            out["Source_SDF"] = instance_sdf

        expanded_rows.append(out)

    expanded = pd.DataFrame(expanded_rows)

    if len(expanded):
        expanded = drop_internal_cols(expanded)

    audit_df = pd.DataFrame(audit_rows)
    return expanded, audit_df


# =============================================================================
# Rebuild SMILES tables
# =============================================================================

def rebuild_smiles_tables(instances: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build:
      Recruiter_SMILES_Map.csv: SMILES, RECRUITER_CODE
      Ligase_SMILE_Codes.csv:  SMILES, RECRUITER_CODE
      Recruiter_SMILES_Wide.csv: SMILES, Recruiter_1...
    """
    base = instances[["SMILES", "RECRUITER_CODE"]].copy()
    base = base.dropna(subset=["SMILES", "RECRUITER_CODE"])
    base["SMILES"] = base["SMILES"].astype(str).str.strip()
    base["RECRUITER_CODE"] = base["RECRUITER_CODE"].astype(str).str.strip()
    base = base[(base["SMILES"] != "") & (base["RECRUITER_CODE"] != "")]
    base = base.drop_duplicates().sort_values(["RECRUITER_CODE", "SMILES"]).reset_index(drop=True)

    recruiter_map = base[["SMILES", "RECRUITER_CODE"]].copy()
    smile_codes = base[["SMILES", "RECRUITER_CODE"]].copy()

    grouped = (
        base.groupby("SMILES")["RECRUITER_CODE"]
        .apply(lambda x: sorted(set(x)))
        .reset_index(name="Recruiter_List")
    )

    max_len = grouped["Recruiter_List"].apply(len).max() if len(grouped) else 0

    wide_data = {"SMILES": grouped["SMILES"]}
    for i in range(max_len):
        wide_data[f"Recruiter_{i + 1}"] = grouped["Recruiter_List"].apply(
            lambda lst, i=i: lst[i] if i < len(lst) else ""
        )

    wide = pd.DataFrame(wide_data)
    return recruiter_map, smile_codes, wide


def build_crosswalk(instances: pd.DataFrame) -> pd.DataFrame:
    cols = [
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
    cols = [c for c in cols if c in instances.columns]
    return drop_internal_cols(instances[cols].copy())


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Propagate LR00001 recruiter codes into SASA/metadata/SMILES tables."
    )

    parser.add_argument("--table-dir", default=str(DEFAULT_TABLE_DIR))
    parser.add_argument("--instances", default=str(DEFAULT_INSTANCE_FILE))
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY_FILE))
    parser.add_argument("--atoms", default=str(DEFAULT_ATOMS_FILE))
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA_FILE))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")

    args = parser.parse_args()

    table_dir = Path(args.table_dir)
    instances_file = Path(args.instances)
    summary_file = Path(args.summary)
    atoms_file = Path(args.atoms)
    metadata_file = Path(args.metadata)

    output_files = [
        summary_file,
        atoms_file,
        metadata_file,
        OUT_RECRUITER_MAP,
        OUT_SMILE_CODES,
        OUT_SMILES_WIDE,
        OUT_CROSSWALK,
        OUT_AUDIT,
        OUT_LOG,
    ]

    if not args.dry_run and not args.overwrite:
        print("❌ Refusing to modify files without --overwrite.")
        print("Preview first:")
        print("  python 10B_CLean.py --dry-run")
        print("Then run:")
        print("  python 10B_CLean.py --overwrite")
        return 2

    for p in [instances_file, summary_file, atoms_file, metadata_file]:
        require_file(p)

    print("=== 10B CLean: Propagate LR Recruiter Codes ===")
    print(f"Started: {now()}")
    print(f"Instances: {instances_file}")
    print(f"Summary:   {summary_file}")
    print(f"Atoms:     {atoms_file}")
    print(f"Metadata:  {metadata_file}")
    print(f"Dry run:   {args.dry_run}")

    instances = load_instances(instances_file)

    summary_df = pd.read_csv(summary_file)
    atoms_df = pd.read_csv(atoms_file)
    meta_df = pd.read_csv(metadata_file)

    print(f"\n📥 Loaded LR instances: {len(instances)}")
    print(f"📥 Loaded summary rows: {len(summary_df)}")
    print(f"📥 Loaded atom rows: {len(atoms_df)}")
    print(f"📥 Loaded metadata rows: {len(meta_df)}")

    summary_clean, audit_summary = update_instance_table(
        summary_df,
        instances,
        table_name="Ligase_Ligand_SASA_summary",
    )

    atoms_clean, audit_atoms = update_instance_table(
        atoms_df,
        instances,
        table_name="Ligase_Ligand_SASA_atoms",
    )

    meta_clean, audit_meta = expand_metadata_to_instances(
        meta_df,
        instances,
    )

    recruiter_map, smile_codes, smiles_wide = rebuild_smiles_tables(instances)
    crosswalk = build_crosswalk(instances)

    audit_parts = [audit_summary, audit_atoms, audit_meta]
    audit_df = pd.concat(
        [x for x in audit_parts if x is not None and len(x)],
        ignore_index=True,
    ) if any(len(x) for x in audit_parts) else pd.DataFrame(
        columns=[
            "Table", "Issue", "Ligase", "pdb_id", "Ligand", "Variant",
            "Old_RECRUITER_CODE", "New_RECRUITER_CODE", "Notes"
        ]
    )

    # Check LR-code consistency.
    summary_codes = set(summary_clean["RECRUITER_CODE"].dropna().astype(str))
    atom_codes = set(atoms_clean["RECRUITER_CODE"].dropna().astype(str))
    meta_codes = set(meta_clean["RECRUITER_CODE"].dropna().astype(str))
    inst_codes = set(instances["RECRUITER_CODE"].dropna().astype(str))

    consistency = {
        "instances_codes": len(inst_codes),
        "summary_codes": len(summary_codes),
        "atom_codes": len(atom_codes),
        "metadata_codes": len(meta_codes),
        "summary_missing_from_instances": sorted(summary_codes - inst_codes),
        "atoms_missing_from_instances": sorted(atom_codes - inst_codes),
        "metadata_missing_from_instances": sorted(meta_codes - inst_codes),
        "instances_missing_from_summary": sorted(inst_codes - summary_codes),
        "instances_missing_from_atoms": sorted(inst_codes - atom_codes),
        "instances_missing_from_metadata": sorted(inst_codes - meta_codes),
    }

    print("\n=== Planned Output Counts ===")
    print(f"Summary rows:      {len(summary_clean)}")
    print(f"Atom rows:         {len(atoms_clean)}")
    print(f"Metadata rows:     {len(meta_clean)}")
    print(f"Recruiter map rows:{len(recruiter_map)}")
    print(f"SMILE code rows:   {len(smile_codes)}")
    print(f"Wide SMILES rows:  {len(smiles_wide)}")
    print(f"Crosswalk rows:    {len(crosswalk)}")
    print(f"Audit rows:        {len(audit_df)}")

    print("\n=== Unique LR Code Counts ===")
    print(f"Instances: {len(inst_codes)}")
    print(f"Summary:   {len(summary_codes)}")
    print(f"Atoms:     {len(atom_codes)}")
    print(f"Metadata:  {len(meta_codes)}")

    if args.dry_run:
        print("\n🧪 Dry run complete. Nothing written.")
        print("Run for real with:")
        print("  python 10B_CLean.py --overwrite")
        return 0

    # Backup old outputs.
    table_dir.mkdir(parents=True, exist_ok=True)
    for p in output_files:
        if p.exists():
            backup = backup_if_exists(p)
            if backup:
                print(f"🧷 Backed up {p} -> {backup}")

    # Write cleaned outputs.
    summary_clean.to_csv(summary_file, index=False)
    atoms_clean.to_csv(atoms_file, index=False)
    meta_clean.to_csv(metadata_file, index=False)

    recruiter_map.to_csv(OUT_RECRUITER_MAP, index=False)
    smile_codes.to_csv(OUT_SMILE_CODES, index=False)
    smiles_wide.to_csv(OUT_SMILES_WIDE, index=False)
    crosswalk.to_csv(OUT_CROSSWALK, index=False)
    audit_df.to_csv(OUT_AUDIT, index=False)

    log_payload = {
        "started_finished": now(),
        "instances_file": str(instances_file),
        "summary_file": str(summary_file),
        "atoms_file": str(atoms_file),
        "metadata_file": str(metadata_file),
        "counts": {
            "instances": len(instances),
            "summary_rows": len(summary_clean),
            "atom_rows": len(atoms_clean),
            "metadata_rows": len(meta_clean),
            "recruiter_map_rows": len(recruiter_map),
            "smile_code_rows": len(smile_codes),
            "wide_smiles_rows": len(smiles_wide),
            "crosswalk_rows": len(crosswalk),
            "audit_rows": len(audit_df),
        },
        "consistency": consistency,
        "outputs": {
            "summary": str(summary_file),
            "atoms": str(atoms_file),
            "metadata": str(metadata_file),
            "recruiter_smiles_map": str(OUT_RECRUITER_MAP),
            "ligase_smile_codes": str(OUT_SMILE_CODES),
            "recruiter_smiles_wide": str(OUT_SMILES_WIDE),
            "crosswalk": str(OUT_CROSSWALK),
            "audit": str(OUT_AUDIT),
        },
    }

    OUT_LOG.write_text(json.dumps(log_payload, indent=2) + "\n")

    print("\n✅ Cleaned LR recruiter codes propagated.")
    print(f"✅ Updated {summary_file}")
    print(f"✅ Updated {atoms_file}")
    print(f"✅ Updated {metadata_file}")
    print(f"✅ Written {OUT_RECRUITER_MAP}")
    print(f"✅ Written {OUT_SMILE_CODES}")
    print(f"✅ Written {OUT_SMILES_WIDE}")
    print(f"✅ Written {OUT_CROSSWALK}")
    print(f"✅ Written {OUT_AUDIT}")
    print(f"🧾 Written {OUT_LOG}")

    if len(audit_df):
        print("\n⚠️ Audit rows were written. Review:")
        print(f"  {OUT_AUDIT}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
