#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute_Advanced_Scaffold_Data.py
---------------------------------
Enhances Ligase_Scaffold_Frequency.csv with:
  ✅ Bemis–Murcko scaffold classification
  ✅ Canonical Murcko_SMILES (core framework)
  ✅ Scaffold_Center_of_Mass_X/Y/Z
  ✅ Ligase_Scaffold_Connectivity
  ✅ Recruiter_Density_Score
  ✅ Shannon_Diversity_Index + Normalized_Diversity

Outputs:
  Ligase_Table/Ligase_Scaffold_Data.csv
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
import numpy as np
import networkx as nx
from math import log
from pathlib import Path

# ============================================================
# ⚙️ Configuration
# ============================================================
INPUT = Path("Ligase_Table/Ligase_Scaffold_Frequency.csv")
OUTPUT = Path("Ligase_Table/Ligase_Scaffold_Data.csv")

# ============================================================
# 🧩 Generate Bemis–Murcko Scaffold + Classify
# ============================================================
def generate_murcko(smiles):
    """Return (Murcko_SMILES, Scaffold_Class) from input SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "", "Unknown"
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        if core.GetNumAtoms() == 0:
            return "", "Acyclic"
        murcko_smi = Chem.MolToSmiles(core, canonical=True)
        ring_info = core.GetRingInfo()
        ring_count = ring_info.NumRings()
        heteroatoms = sum(1 for a in core.GetAtoms() if a.GetAtomicNum() not in (6, 1))

        if ring_count == 0:
            cls = "Acyclic"
        elif ring_count == 1:
            cls = "Monocyclic_Hetero" if heteroatoms > 0 else "Monocyclic"
        elif ring_count == 2:
            cls = "Bicyclic_Hetero" if heteroatoms > 0 else "Bicyclic"
        else:
            cls = "Polycyclic_Hetero" if heteroatoms > 0 else "Polycyclic"

        return murcko_smi, cls
    except Exception:
        return "", "Unknown"

# ============================================================
# ⚛️ Compute Scaffold Center of Mass
# ============================================================
def scaffold_center_of_mass(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan, np.nan, np.nan
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        conf = mol.GetConformer()
        coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        center = coords.mean(axis=0)
        return tuple(center)
    except Exception:
        return np.nan, np.nan, np.nan

# ============================================================
# 🔗 Ligase_Scaffold_Connectivity
# ============================================================
def compute_connectivity(df):
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row["Ligase"], row["Scaffold_ID"])
    connectivity = {node: G.degree(node) for node in G.nodes() if "_SCAF_" in node}
    df["Ligase_Scaffold_Connectivity"] = df["Scaffold_ID"].map(connectivity)
    return df

# ============================================================
# 📈 Recruiter_Density_Score
# ============================================================
def compute_density(df):
    ligase_totals = df.groupby("Ligase")["Recruiter_Count"].sum()
    df["Total_Recruiters"] = df["Ligase"].map(ligase_totals)
    df["Recruiter_Density_Score"] = df["Recruiter_Count"] / df["Total_Recruiters"]
    return df

# ============================================================
# 🔢 Shannon Diversity Index + Normalized Diversity
# ============================================================
def compute_shannon_index(df):
    ligase_groups = df.groupby("Ligase")
    shannon = {}
    normalized = {}
    for ligase, group in ligase_groups:
        counts = group["Recruiter_Count"].values
        total = counts.sum()
        p = counts / total
        H = -sum(p_i * log(p_i) for p_i in p if p_i > 0)
        Hmax = log(len(p)) if len(p) > 1 else 0
        shannon[ligase] = H
        normalized[ligase] = H / Hmax if Hmax > 0 else 0
    df["Shannon_Diversity_Index"] = df["Ligase"].map(shannon)
    df["Normalized_Diversity"] = df["Ligase"].map(normalized)
    return df

# ============================================================
# 🚀 Main Routine
# ============================================================
def main():
    df = pd.read_csv(INPUT)
    print(f"🔍 Loaded {len(df)} scaffolds from {INPUT.name}")

    smiles_file = Path("Ligase_Table/Ligase_Recruiters_Scaffold.csv")
    smiles_map = {}
    if smiles_file.exists():
        sdf = pd.read_csv(smiles_file)
        if "Scaffold_SMILES" in sdf.columns:
            smiles_map = dict(zip(sdf["Scaffold_ID"], sdf["Scaffold_SMILES"]))

    df["Scaffold_SMILES"] = df["Scaffold_ID"].map(smiles_map).fillna("")

    # Generate Murcko scaffold + classification
    murcko_smiles, scaffold_class = [], []
    com_x, com_y, com_z = [], [], []
    for smi in df["Scaffold_SMILES"]:
        murcko, cls = generate_murcko(smi)
        murcko_smiles.append(murcko)
        scaffold_class.append(cls)
        x, y, z = scaffold_center_of_mass(smi)
        com_x.append(x)
        com_y.append(y)
        com_z.append(z)

    df["Murcko_SMILES"] = murcko_smiles
    df["Scaffold_Class"] = scaffold_class
    df["Scaffold_Center_of_Mass_X"] = com_x
    df["Scaffold_Center_of_Mass_Y"] = com_y
    df["Scaffold_Center_of_Mass_Z"] = com_z

    df = compute_connectivity(df)
    df = compute_density(df)
    df = compute_shannon_index(df)

    df.to_csv(OUTPUT, index=False)
    print(f"✅ Saved enhanced scaffold data → {OUTPUT}")
    print("🧪 Added: Murcko_SMILES, Scaffold_Class, Center_of_Mass_(X/Y/Z), "
          "Connectivity, Recruiter_Density_Score, Shannon_Diversity_Index, Normalized_Diversity")

if __name__ == "__main__":
    main()
