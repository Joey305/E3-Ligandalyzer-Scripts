#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Remove_From_All_Ligase_Tables_and_Structures.py
===============================================================================
Author: Joseph-Michael Schulz (UM BMB)
Purpose:
    Removes all rows that contain a given ligand or pdb_id (or both)
    from all three Ligase tables inside ./Ligase_Table/,
    AND deletes any matching structure files (.pdb, .sdf)
    from all relevant subdirectories.

Input:
    Ligase_Table/Ligase_Ligand_SASA_summary.csv
    Ligase_Table/Ligase_Ligand_SASA_atoms.csv
    Ligase_Table/Ligase_Ligand_Metadata.csv

Usage examples:
    python Remove_From_All_Ligase_Tables_and_Structures.py --ligand DTR
    python Remove_From_All_Ligase_Tables_and_Structures.py --pdb 4VHK
    python Remove_From_All_Ligase_Tables_and_Structures.py --ligand DTR --pdb 4VHK
    python Remove_From_All_Ligase_Tables_and_Structures.py      ← interactive mode

Output:
    Overwrites filtered CSVs and deletes structure files.
    Logs all actions in: Ligase_Table/Remove_Filter.log
===============================================================================
"""

import os
import pandas as pd
import argparse
import re

# -------------------- CLI setup --------------------
parser = argparse.ArgumentParser(description="Remove entries by ligand or pdb_id from all Ligase tables and structure folders.")
parser.add_argument("--ligand", help="Ligand 3-letter code(s), comma-separated")
parser.add_argument("--pdb", help="PDB ID(s), comma-separated")
args = parser.parse_args()

ligands = []
pdbs = []

# Interactive prompt if no args given
if not args.ligand and not args.pdb:
    print("⚙️ No ligand or pdb_id flag provided. Enter manually:")
    mode = input("Remove by (L)igand, (P)DB, or (B)oth? ").strip().lower()
    if mode in ["l", "ligand"]:
        args.ligand = input("Enter ligand code(s), comma-separated: ").strip()
    elif mode in ["p", "pdb"]:
        args.pdb = input("Enter PDB ID(s), comma-separated: ").strip()
    elif mode in ["b", "both"]:
        args.ligand = input("Enter ligand code(s), comma-separated: ").strip()
        args.pdb = input("Enter PDB ID(s), comma-separated: ").strip()
    else:
        print("❌ Invalid option. Exiting.")
        exit()

if args.ligand:
    ligands = [x.strip().upper() for x in args.ligand.split(",") if x.strip()]
if args.pdb:
    pdbs = [x.strip().upper() for x in args.pdb.split(",") if x.strip()]

if not ligands and not pdbs:
    print("⚠️ Nothing to remove. Exiting.")
    exit()

# -------------------- Setup --------------------
table_dir = "Ligase_Table"
log_file = os.path.join(table_dir, "Remove_Filter.log")

files = {
    "summary": os.path.join(table_dir, "Ligase_Ligand_SASA_summary.csv"),
    "atoms":   os.path.join(table_dir, "Ligase_Ligand_SASA_atoms.csv"),
    "meta":    os.path.join(table_dir, "Ligase_Ligand_Metadata.csv")
}

# -------------------- Helpers --------------------
def remove_entries(path, ligands, pdbs):
    """Remove matching rows from CSV if it exists"""
    if not os.path.exists(path):
        print(f"⚠️ File not found: {path}")
        return 0, 0

    df = pd.read_csv(path)
    if df.empty:
        print(f"⚠️ Skipped empty file: {path}")
        return 0, 0

    before = len(df)
    cond = pd.Series(False, index=df.index)

    if ligands and "Ligand" in df.columns:
        cond |= df["Ligand"].astype(str).str.upper().isin(ligands)
    if pdbs and "pdb_id" in df.columns:
        cond |= df["pdb_id"].astype(str).str.upper().isin(pdbs)

    removed = df[cond]
    kept = df[~cond]
    kept.to_csv(path, index=False)
    return before, len(removed)


def match_pdb_filename(filename, ligands, pdbs):
    """Check if a file matches ligand or pdb ID pattern."""
    name = filename.upper()

    # Matches like 4VHK_LIG.pdb or 4VHK_LIG_1.pdb
    if pdbs:
        for pdb in pdbs:
            if re.match(rf"^{pdb}_[A-Z0-9]+(_\d+)?\.PDB$", name):
                return True

    # Matches like ####_LIG.pdb or ####_LIG_1.pdb
    if ligands:
        for lig in ligands:
            if re.match(rf"^[A-Z0-9]+_{lig}(_\d+)?\.PDB$", name):
                return True

    return False


def match_sdf_filename(filename, ligands):
    """Check if a .sdf file matches pattern SOMETHING(LIG)_LIG.sdf"""
    name = filename.upper()
    if not name.endswith(".SDF"):
        return False
    for lig in ligands:
        if f"({lig})_LIG.SDF" in name:
            return True
    return False


def delete_structure_files(base_dir, ligands, pdbs):
    """Recursively delete .pdb and .sdf files that match ligand/pdb rules"""
    deleted_files = []

    for root, dirs, files in os.walk(base_dir):
        for f in files:
            f_upper = f.upper()
            full_path = os.path.join(root, f)

            # Handle .pdb
            if f_upper.endswith(".PDB") and match_pdb_filename(f_upper, ligands, pdbs):
                os.remove(full_path)
                deleted_files.append(full_path)
                continue

            # Handle .sdf
            if ligands and f_upper.endswith(".SDF") and match_sdf_filename(f_upper, ligands):
                os.remove(full_path)
                deleted_files.append(full_path)

    return deleted_files

# -------------------- Run filtering --------------------
summary_total, summary_removed = remove_entries(files["summary"], ligands, pdbs)
atoms_total, atoms_removed     = remove_entries(files["atoms"], ligands, pdbs)
meta_total, meta_removed       = remove_entries(files["meta"], ligands, pdbs)

deleted_structures = delete_structure_files(".", ligands, pdbs)

# -------------------- Logging --------------------
with open(log_file, "a") as log:
    log.write("\n=== Removal Run ===\n")
    if ligands:
        log.write(f"Removed Ligands: {', '.join(ligands)}\n")
    if pdbs:
        log.write(f"Removed PDB IDs: {', '.join(pdbs)}\n")
    log.write(f"Summary: {summary_removed}/{summary_total} rows removed\n")
    log.write(f"Atoms:   {atoms_removed}/{atoms_total} rows removed\n")
    log.write(f"Metadata:{meta_removed}/{meta_total} rows removed\n")
    if deleted_structures:
        log.write(f"\nDeleted structure files ({len(deleted_structures)}):\n")
        for p in deleted_structures:
            log.write(f" - {p}\n")
    else:
        log.write("\nNo matching PDB or SDF files found.\n")
    log.write("\n")

print(f"✅ Done! Updated Ligase_Table and removed {len(deleted_structures)} structure files.")
print(f"🧾 Logged changes in {log_file}")
