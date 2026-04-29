#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Build_Ligase_Table_Final.py
===============================================================================
Purpose:
    Finalize clean ligase tables after:
      - PDB cleanup
      - duplicate ligand renaming
      - SASA rerun
      - SDF rebuild
      - metadata regeneration

This safer version:
    1. Removes BRD Recruiter rows from SASA summary using exact row identity.
    2. Filters atom rows to match the remaining SASA summary rows.
    3. Filters metadata using RECRUITER_CODE = Ligase_Ligand.
    4. Ensures final SASA tables only contain ligands with validated metadata.
    5. Writes a filtered Ligase_Table/Ligase_SMILE_Codes.csv for downstream scripts.

Outputs:
    Ligase_Table/
      ├── Ligase_Ligand_SASA_atoms.csv
      ├── Ligase_Ligand_SASA_summary.csv
      ├── Ligase_Ligand_Metadata.csv
      ├── Ligase_SMILE_Codes.csv
      ├── Final_Table_Audit.csv

    Ligase_Filter.log
===============================================================================
"""

import os
import pandas as pd
from pathlib import Path

# ====================== Setup ==========================

summary_file = Path("Ligand_SASA_summary.csv")
atoms_file   = Path("Ligand_SASA_atoms.csv")
meta_file    = Path("Ligand_Metadata.csv")
smiles_file  = Path("Ligase_Table/Ligase_SMILE_Codes.csv")

output_dir   = Path("Ligase_Table")
log_file     = Path("Ligase_Filter.log")

output_dir.mkdir(exist_ok=True)


# ====================== Helpers ==========================

def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")


def normalize_common_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    for col in ["Ligase", "Ligand", "pdb_id", "Variant"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    if "Ligand" in df.columns:
        df["Ligand"] = df["Ligand"].str.upper()

    if "Variant" not in df.columns:
        df["Variant"] = "1"
    else:
        df["Variant"] = df["Variant"].replace({"nan": "1", "NA": "1", "": "1"}).fillna("1").astype(str)

    if {"Ligase", "Ligand"}.issubset(df.columns):
        df["RECRUITER_CODE"] = df["Ligase"].astype(str).str.strip() + "_" + df["Ligand"].astype(str).str.strip()

    if {"Ligase", "pdb_id", "Ligand", "Variant"}.issubset(df.columns):
        df["_ROW_KEY"] = (
            df["Ligase"].astype(str).str.strip() + "|" +
            df["pdb_id"].astype(str).str.strip() + "|" +
            df["Ligand"].astype(str).str.strip() + "|" +
            df["Variant"].astype(str).str.strip()
        )

    return df


def write_count_block(log, title, series):
    log.write(f"\n=== {title} ===\n")
    if series is None or len(series) == 0:
        log.write("None\n")
        return
    for key, value in series.items():
        log.write(f"{str(key):<20} : {value}\n")


# ====================== Load Files =====================

for f in [summary_file, atoms_file, meta_file]:
    require_file(f)

summary_df = normalize_common_cols(pd.read_csv(summary_file))
atoms_df   = normalize_common_cols(pd.read_csv(atoms_file))
meta_df    = normalize_common_cols(pd.read_csv(meta_file))

print(
    f"📥 Loaded {len(summary_df)} summary rows, "
    f"{len(atoms_df)} atom rows, and {len(meta_df)} metadata rows."
)

if "RECRUITER_CODE" not in meta_df.columns:
    if {"Ligase", "Ligand"}.issubset(meta_df.columns):
        meta_df["RECRUITER_CODE"] = meta_df["Ligase"] + "_" + meta_df["Ligand"]
    else:
        raise ValueError("Metadata file must contain either RECRUITER_CODE or Ligase + Ligand columns.")

# ====================== Remove BRD Recruiter Rows =====================

if "Recruiter_Class" in summary_df.columns:
    brd_mask = summary_df["Recruiter_Class"].astype(str).str.contains("BRD Recruiter", case=False, na=False)
    brd_rows = summary_df[brd_mask].copy()
else:
    brd_mask = pd.Series(False, index=summary_df.index)
    brd_rows = summary_df.iloc[0:0].copy()
    print("⚠️ No Recruiter_Class column found in summary — skipping BRD filter.")

brd_row_keys = set(brd_rows["_ROW_KEY"]) if "_ROW_KEY" in brd_rows.columns else set()

summary_no_brd = summary_df[~summary_df["_ROW_KEY"].isin(brd_row_keys)].copy()

print(f"🚫 Removed {len(brd_rows)} exact BRD summary rows.")

# ====================== Keep Only Rows With Validated Metadata =====================

metadata_codes = set(meta_df["RECRUITER_CODE"].dropna().astype(str))
summary_codes_after_brd = set(summary_no_brd["RECRUITER_CODE"].dropna().astype(str))

codes_missing_metadata = sorted(summary_codes_after_brd - metadata_codes)
metadata_without_summary = sorted(metadata_codes - summary_codes_after_brd)

final_codes = summary_codes_after_brd & metadata_codes

summary_filtered = summary_no_brd[summary_no_brd["RECRUITER_CODE"].isin(final_codes)].copy()

# Filter atoms by exact summary rows remaining.
final_row_keys = set(summary_filtered["_ROW_KEY"].dropna().astype(str))
atoms_filtered = atoms_df[atoms_df["_ROW_KEY"].isin(final_row_keys)].copy()

# Filter metadata to final active recruiter codes.
meta_filtered = meta_df[meta_df["RECRUITER_CODE"].isin(final_codes)].copy()

print(f"✅ Final summary rows:  {len(summary_filtered)}")
print(f"✅ Final atom rows:     {len(atoms_filtered)}")
print(f"✅ Final metadata rows: {len(meta_filtered)}")
print(f"⚠️ Summary codes missing metadata and dropped: {len(codes_missing_metadata)}")
print(f"ℹ️ Metadata codes without remaining summary rows: {len(metadata_without_summary)}")

# ====================== Filter / Rebuild SMILES Table =====================

smiles_filtered = None

if smiles_file.exists():
    smiles_df = pd.read_csv(smiles_file)
    if "RECRUITER_CODE" not in smiles_df.columns:
        if {"Ligase", "Ligand"}.issubset(smiles_df.columns):
            smiles_df["RECRUITER_CODE"] = smiles_df["Ligase"].astype(str) + "_" + smiles_df["Ligand"].astype(str)
        else:
            raise ValueError(f"{smiles_file} lacks RECRUITER_CODE and Ligase/Ligand columns.")

    smiles_filtered = smiles_df[smiles_df["RECRUITER_CODE"].isin(final_codes)].copy()
else:
    # Rebuild a minimal compatible smiles table from metadata.
    needed = [
        "RECRUITER_CODE", "Ligase", "Ligand", "SMILES", "Canonical_SMILES",
        "InChIKey", "Formula", "MW", "Source_SDF",
        "SDF_Heavy_Atom_Count", "SMILES_Heavy_Atom_Count",
        "SASA_Atom_Count_Source", "SASA_Atom_Counts_Unique",
        "Atom_Count_Validation", "Validation_Notes",
    ]
    smiles_filtered = meta_filtered[[c for c in needed if c in meta_filtered.columns]].copy()

# ====================== Save Outputs =====================

drop_internal = ["_ROW_KEY"]

def clean_for_output(df):
    return df.drop(columns=[c for c in drop_internal if c in df.columns], errors="ignore")

summary_out = output_dir / "Ligase_Ligand_SASA_summary.csv"
atoms_out   = output_dir / "Ligase_Ligand_SASA_atoms.csv"
meta_out    = output_dir / "Ligase_Ligand_Metadata.csv"
smiles_out  = output_dir / "Ligase_SMILE_Codes.csv"
audit_out   = output_dir / "Final_Table_Audit.csv"

clean_for_output(summary_filtered).to_csv(summary_out, index=False)
clean_for_output(atoms_filtered).to_csv(atoms_out, index=False)
clean_for_output(meta_filtered).to_csv(meta_out, index=False)
clean_for_output(smiles_filtered).to_csv(smiles_out, index=False)

audit_rows = []

for code in codes_missing_metadata:
    audit_rows.append({
        "Issue": "SUMMARY_CODE_DROPPED_NO_METADATA",
        "RECRUITER_CODE": code,
        "Notes": "Present in non-BRD SASA summary but absent from validated metadata/SMILES table.",
    })

for code in metadata_without_summary:
    audit_rows.append({
        "Issue": "METADATA_CODE_UNUSED",
        "RECRUITER_CODE": code,
        "Notes": "Present in metadata but absent from non-BRD SASA summary after filtering.",
    })

for _, r in brd_rows.iterrows():
    audit_rows.append({
        "Issue": "BRD_ROW_REMOVED",
        "RECRUITER_CODE": r.get("RECRUITER_CODE", ""),
        "Notes": f"{r.get('Ligase', '')} {r.get('pdb_id', '')} {r.get('Ligand', '')} variant {r.get('Variant', '')}",
    })

pd.DataFrame(audit_rows).to_csv(audit_out, index=False)

print("✅ Cleaned datasets saved in Ligase_Table/")
print(f"✅ Written {summary_out}")
print(f"✅ Written {atoms_out}")
print(f"✅ Written {meta_out}")
print(f"✅ Written {smiles_out}")
print(f"✅ Written {audit_out}")

# ====================== Logging =====================

ligase_counts = summary_filtered["Ligase"].value_counts().sort_index() if "Ligase" in summary_filtered.columns else pd.Series(dtype=int)
class_counts = summary_filtered["Recruiter_Class"].value_counts() if "Recruiter_Class" in summary_filtered.columns else pd.Series(dtype=int)

with open(log_file, "w") as log:
    log.write("=== Ligase Table Finalization Log ===\n\n")

    log.write(f"Input summary rows:  {len(summary_df)}\n")
    log.write(f"Input atom rows:     {len(atoms_df)}\n")
    log.write(f"Input metadata rows: {len(meta_df)}\n\n")

    log.write(f"Exact BRD summary rows removed: {len(brd_rows)}\n")
    log.write(f"Summary codes missing metadata and dropped: {len(codes_missing_metadata)}\n")
    log.write(f"Metadata codes without remaining summary rows: {len(metadata_without_summary)}\n\n")

    log.write(f"Final summary rows:  {len(summary_filtered)}\n")
    log.write(f"Final atom rows:     {len(atoms_filtered)}\n")
    log.write(f"Final metadata rows: {len(meta_filtered)}\n")
    log.write(f"Final SMILES rows:   {len(smiles_filtered)}\n\n")

    if codes_missing_metadata:
        log.write("=== Summary codes dropped because metadata was missing ===\n")
        for code in codes_missing_metadata:
            log.write(f" - {code}\n")
        log.write("\n")

    if metadata_without_summary:
        log.write("=== Metadata codes not used after summary filtering ===\n")
        for code in metadata_without_summary:
            log.write(f" - {code}\n")
        log.write("\n")

    write_count_block(log, "Remaining Summary Rows per Ligase", ligase_counts)
    write_count_block(log, "Remaining Recruiter_Class Counts", class_counts)

print(f"🧾 Log written to {log_file}")