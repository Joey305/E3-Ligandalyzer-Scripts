#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligase_Recruiter_Scaffolds_FinalPlus.py
----------------------------------------
Generates scaffold groupings per ligase using recruiter–SMILES relationships.
All recruiters are retained — malformed SMILES simply receive unique hash scaffolds.

Outputs:
  1️⃣ Ligase_Table/Ligase_Recruiters_Scaffold.csv
  2️⃣ Ligase_Table/Ligase_Scaffold_Summary.csv
  3️⃣ Ligase_Table/Ligase_Scaffold_Frequency.csv
  4️⃣ Ligase_Table/Ligase_Recruiters_Scaffold.log

Diversity metrics:
  - Diversity_Score = unique scaffolds / total recruiters
  - Shannon_Index   = -Σ(pᵢ * ln(pᵢ)) across scaffolds per ligase
"""

import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import hashlib, math
from collections import Counter

# ============================================================
# ⚙️ Configuration
# ============================================================
TABLE_DIR = Path("Ligase_Table")
TABLE_DIR.mkdir(exist_ok=True)

MAP3D = TABLE_DIR / "Ligase_Ligands_Smiles_3DMapped.csv"
SMILESCODES = TABLE_DIR / "Ligase_SMILE_Codes.csv"

OUT1 = TABLE_DIR / "Ligase_Recruiters_Scaffold.csv"
OUT2 = TABLE_DIR / "Ligase_Scaffold_Summary.csv"
OUT3 = TABLE_DIR / "Ligase_Scaffold_Frequency.csv"
LOG  = TABLE_DIR / "Ligase_Recruiters_Scaffold.log"

# ============================================================
# 🧩 Helper Functions
# ============================================================
def safe_mol(smi):
    try:
        return Chem.MolFromSmiles(smi)
    except Exception:
        return None

def scaffold_smiles(smi):
    mol = safe_mol(smi)
    if not mol:
        return None
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaf, canonical=True) if scaf else None
    except Exception:
        return None

def hash_string(text):
    """Deterministic short hash."""
    return hashlib.sha1(str(text).encode()).hexdigest()[:8]

def shannon_diversity(counts):
    """Compute Shannon index for diversity."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((n / total) * math.log(n / total) for n in counts.values() if n > 0)

# ============================================================
# 🚀 Main Routine
# ============================================================
def main():
    # Load source data
    df_3d = pd.read_csv(MAP3D)
    df_smiles = pd.read_csv(SMILESCODES)
    print(f"🔍 Loaded {len(df_smiles)} SMILES and {len(df_3d)} 3D mapping rows.")

    recruiter_map = df_3d[["Ligase", "RECRUITER_CODE"]].drop_duplicates()
    df = pd.merge(recruiter_map, df_smiles, on="RECRUITER_CODE", how="left")
    df = df.dropna(subset=["SMILES", "Ligase", "RECRUITER_CODE"])

    # Compute scaffolds
    df["Scaffold_SMILES"] = df["SMILES"].apply(scaffold_smiles)
    df["Scaffold_Hash"] = df["Scaffold_SMILES"].apply(lambda x: hash_string(x) if x else None)

    out_rows, summary_rows, freq_rows = [], [], []

    for ligase, sub in df.groupby("Ligase"):
        valid_scaffolds = sub["Scaffold_Hash"].dropna().unique().tolist()
        scaffold_map = {sc: f"{ligase}_SCAF_{i+1}" for i, sc in enumerate(valid_scaffolds)}

        for _, row in sub.iterrows():
            scaf_hash = row["Scaffold_Hash"]
            if pd.isna(scaf_hash):
                # Unique hashed scaffold fallback
                hash_val = hash_string(row["SMILES"])
                scaf_id = f"{ligase}_SCAF_{hash_val}"
                scaf_smi = row["SMILES"]
                scaf_hash = hash_val
            else:
                scaf_id = scaffold_map[scaf_hash]
                scaf_smi = row["Scaffold_SMILES"]

            out_rows.append({
                "Ligase": ligase,
                "Scaffold_ID": scaf_id,
                "RECRUITER_CODE": row["RECRUITER_CODE"],
                "Scaffold_SMILES": scaf_smi,
                "Scaffold_Hash": scaf_hash
            })

        # Frequency table
        sub_freq = (
            pd.DataFrame(out_rows)
            .query("Ligase == @ligase")
            .groupby(["Ligase", "Scaffold_ID"])
            .size()
            .reset_index(name="Recruiter_Count")
        )
        freq_rows.extend(sub_freq.to_dict("records"))

        # Diversity summary
        counts = Counter(sub_freq["Recruiter_Count"])
        total = sub_freq["Recruiter_Count"].sum()
        unique = len(sub_freq)
        diversity = unique / total if total > 0 else 0
        shannon = shannon_diversity(Counter(sub_freq["Recruiter_Count"]))

        summary_rows.append({
            "Ligase": ligase,
            "Unique_Scaffolds": unique,
            "Total_Recruiters": total,
            "Diversity_Score": round(diversity, 3),
            "Shannon_Index": round(shannon, 3)
        })

    # Save outputs
    pd.DataFrame(out_rows).to_csv(OUT1, index=False)
    pd.DataFrame(summary_rows).to_csv(OUT2, index=False)
    pd.DataFrame(freq_rows).to_csv(OUT3, index=False)

    with open(LOG, "w") as f:
        f.write(f"✅ {len(out_rows)} recruiter–scaffold mappings written.\n")
        f.write(f"🔹 {len(summary_rows)} ligase diversity summaries.\n")
        f.write(f"📊 {len(freq_rows)} scaffold frequency rows.\n")

    print(f"✅ Written recruiter-scaffold mapping → {OUT1}")
    print(f"✅ Written scaffold summary → {OUT2}")
    print(f"✅ Written scaffold frequency table → {OUT3}")
    print(f"🧩 All recruiters assigned deterministic scaffold IDs — no data lost.")

if __name__ == "__main__":
    main()
