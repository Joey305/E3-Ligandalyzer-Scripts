#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map_Recruiter_SMILES.py
----------------------------------------
Builds a cross-reference table linking each unique SMILES string
to all recruiter IDs (RECRUITER_CODE) that use it.

Input:  Ligase_SMILE_Codes.csv
Output: 
  1️⃣ Recruiter_SMILES_Map.csv   (SMILES + comma-separated RECRUITER_CODEs)
  2️⃣ Recruiter_SMILES_Wide.csv  (wide matrix: one column per recruiter)
"""

import pandas as pd
from pathlib import Path

# ============================================================
# ⚙️ Configuration
# ============================================================
INPUT  = Path("Ligase_Table/Ligase_SMILE_Codes.csv")
OUT1   = Path("Ligase_Table/Recruiter_SMILES_Map.csv")
OUT2   = Path("Ligase_Table/Recruiter_SMILES_Wide.csv")

# ============================================================
# 🚀 Load Data
# ============================================================
df = pd.read_csv(INPUT)
df = df.dropna(subset=["SMILES", "RECRUITER_CODE"])

print(f"🔍 Loaded {len(df)} entries from {INPUT.name}")
print(f"🔹 Found {df['SMILES'].nunique()} unique SMILES")

# ============================================================
# 🧩 Group by SMILES
# ============================================================
grouped = (
    df.groupby("SMILES")["RECRUITER_CODE"]
    .apply(lambda x: sorted(set(x)))  # remove duplicates, sort
    .reset_index(name="Recruiter_List")
)

# Join as comma-separated string
grouped["Recruiters_CommaSeparated"] = grouped["Recruiter_List"].apply(lambda x: ", ".join(x))

# ============================================================
# 💾 Save Format 1 (Compact)
# ============================================================
grouped[["SMILES", "Recruiters_CommaSeparated"]].to_csv(OUT1, index=False)
print(f"✅ Written {OUT1}")

# ============================================================
# 💾 Save Format 2 (Wide Matrix)
# ============================================================
# Pad lists to equal length for DataFrame expansion
max_len = grouped["Recruiter_List"].apply(len).max()
wide = pd.DataFrame({
    "SMILES": grouped["SMILES"],
    **{
        f"Recruiter_{i+1}": grouped["Recruiter_List"].apply(
            lambda lst, i=i: lst[i] if i < len(lst) else ""
        )
        for i in range(max_len)
    }
})

wide.to_csv(OUT2, index=False)
print(f"✅ Written {OUT2}")

print(f"🧪 {len(grouped)} unique SMILES processed.")
