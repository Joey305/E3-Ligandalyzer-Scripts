#!/usr/bin/env python3
"""
CleanSASA.py

Deletes all files whose names match:
   *_CAF.pdb, *_CAF_*.pdb, *_CAF.sdf

Searches recursively under current directory.
"""

from pathlib import Path
import re
import os

# Define filename patterns to remove
patterns = [r"_CAF\.pdb$", r"_CAF_\d+\.pdb$", r"_CAF\.sdf$"]
regexes = [re.compile(p, re.IGNORECASE) for p in patterns]  # ✅ fixed typo

removed = []

for root, _, files in os.walk("."):
    for f in files:
        for rx in regexes:
            if rx.search(f):
                full = Path(root) / f
                removed.append(full)
                try:
                    os.remove(full)
                    print(f"🗑️ Removed {full}")
                except Exception as e:
                    print(f"⚠️ Could not delete {full}: {e}")
                break  # don’t check other regexes once matched

print(f"\n🧹 Removed {len(removed)} _CAF-related files.\n")

if removed:
    with open("Removed_CAF_Files.txt", "w") as out:
        for path in removed:
            out.write(str(path) + "\n")
    print("🧾 Written Removed_CAF_Files.txt for audit.")
else:
    print("✅ No _CAF files found.")



#!/usr/bin/env python3
"""
reclassify_by_mw.py
Re-labels recruiter class purely based on MW thresholds.

Classes:
  < 200 → Fragment Recruiter
  200–500 → Drug-like Recruiter
  > 500 → Peptide-like Recruiter

Overwrites Ligand_SASA_summary.csv in place.
"""

import pandas as pd

FILE = "Ligand_SASA_summary.csv"

def classify_mw(mw):
    try:
        mw = float(mw)
    except Exception:
        return "Unknown"
    if mw < 200:
        return "Fragment Recruiter"
    elif mw <= 500:
        return "Drug-like Recruiter"
    else:
        return "Peptide-like Recruiter"

def main():
    print(f"📄 Loading {FILE} ...")
    df = pd.read_csv(FILE)
    print(f"✅ Loaded {len(df)} rows.")

    df["Recruiter_Class"] = df["MW"].apply(classify_mw)

    # Show counts
    counts = df["Recruiter_Class"].value_counts().to_dict()
    print("\n📊 Updated Recruiter_Class counts:")
    for k, v in counts.items():
        print(f"  {k}: {v}")

    df.to_csv(FILE, index=False)
    print(f"\n💾 Overwritten file: {FILE}")
    print("🎯 Classification complete.\n")

if __name__ == "__main__":
    main()

"""
update_brd_from_summary.py

Post-processing for Ligand_SASA_summary.csv:
- If a single pdb_id appears in 4+ rows → classify as 'BRD Recruiter'
- Prints summary to console and writes BRD_PDB_Report.txt
- Overwrites the input CSV with updated Recruiter_Class values

Usage:
    python update_brd_from_summary.py --input Ligand_SASA_summary.csv
"""

import pandas as pd
import argparse
from collections import Counter, defaultdict

def main():
    parser = argparse.ArgumentParser(description="Reclassify BRD recruiters based on PDB repetition.")
    parser.add_argument("--input", default="Ligand_SASA_summary.csv", help="Path to SASA summary CSV file.")
    args = parser.parse_args()

    # Load CSV
    df = pd.read_csv(args.input)
    print(f"📄 Loaded {len(df)} rows from {args.input}\n")

    # Count number of entries per PDB
    pdb_counts = Counter(df["pdb_id"])

    # Identify BRD candidates (≥4 occurrences)
    brd_candidates = [pdb for pdb, count in pdb_counts.items() if count >= 4]
    print(f"🔍 Found {len(brd_candidates)} PDBs with ≥4 ligands/fragments.\n")

    # Build ligase grouping for report
    ligase_groups = defaultdict(list)
    for pdb in brd_candidates:
        ligases = df.loc[df["pdb_id"] == pdb, "Ligase"].unique().tolist()
        for ligase in ligases:
            ligase_groups[ligase].append(pdb)

    # Update recruiter class
    before = df["Recruiter_Class"].value_counts().to_dict()
    df.loc[df["pdb_id"].isin(brd_candidates), "Recruiter_Class"] = "BRD Recruiter"
    after = df["Recruiter_Class"].value_counts().to_dict()

    # Overwrite CSV
    df.to_csv(args.input, index=False)
    print(f"✅ Updated Recruiter_Class for {len(brd_candidates)} PDB structures.")
    print(f"💾 Overwritten file: {args.input}\n")

    # Write BRD report
    report_path = "BRD_PDB_Report.txt"
    with open(report_path, "w") as f:
        f.write("BRD Recruiter Identification Report\n")
        f.write("=================================\n\n")
        f.write(f"Total PDBs classified as BRD Recruiters: {len(brd_candidates)}\n\n")

        for ligase, pdb_list in sorted(ligase_groups.items()):
            f.write(f"{ligase} ({len(pdb_list)} structures):\n")
            for pdb in sorted(pdb_list):
                count = pdb_counts[pdb]
                f.write(f"  - {pdb}  ({count} ligands)\n")
            f.write("\n")

    print("🧾 BRD PDB Report Written → BRD_PDB_Report.txt\n")

    # Print summary preview to console
    print("📊 Summary (Recruiter_Class counts):")
    print("Before update:", before)
    print("After update:", after)
    print("\n📋 BRD Candidates by Ligase:")
    for ligase, pdb_list in sorted(ligase_groups.items()):
        print(f"  {ligase}: {', '.join(sorted(pdb_list))}")

if __name__ == "__main__":
    main()