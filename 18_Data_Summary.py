#!/usr/bin/env python3
"""
18_Data_Summary.py

Generate a manuscript-ready text summary for the E3 Recruiter Ligandalyzer dataset.

IMPORTANT DATA MODEL USED BY THIS VERSION
----------------------------------------
This version fixes the key counting issue:

    Unique ligand = the original PDB/component ligand name, for example EF2, VH5, LEN, etc.
    Recruiter code = the LR##### tag assigned to a specific recruiter entry/conformation/instance.

Therefore:

    - 602 LR codes does NOT mean 602 unique ligands.
    - A single ligand can have multiple LR codes when the same ligand appears in multiple
      structures, poses, crystal conformations, chains, or curated recruiter instances.
    - Manuscript-facing unique ligand counts should be based on the ligand/component name,
      not the LR recruiter code.

Primary source for this mapping:

    Ligase_Table/Recruiter_Master_Map.csv

Expected useful columns in Recruiter_Master_Map.csv may include names like:

    RECRUITER_CODE, Recruiter_Code, Recruiter, LR_Code, LR_ID
    LIGAND, Ligand, Ligand_Name, Ligand_Code, HET, Component_ID, Original_Ligand_Code

The script is defensive and tries to infer columns, but it strongly prefers explicit
RECRUITER_CODE and LIGAND columns when they exist.

Run from the E3Ligandalyzer project root:

    python 18_Data_Summary.py

Optional:

    python 18_Data_Summary.py --table-dir Ligase_Table --ligase-dir Ligases --out E3Recruiter_Data_Summary.txt --top-n 25

Outputs:

    E3Recruiter_Data_Summary.txt
    Ligase_Table/Derived_Ligand_Conformation_Summary.csv
    Ligase_Table/Derived_RecruiterCode_To_Ligand_Map.csv

Purpose:

    Summarizes:
      - CSV table inventory
      - curated number of E3 ligases represented
      - unique ligand/component counts based on ligand name, not LR code
      - LR recruiter code / conformation / instance counts
      - number of conformations per ligand
      - unique PDB/structure counts
      - scaffold and supercluster counts
      - per-ligase ligand/recruiter/scaffold summaries
      - physicochemical descriptor summaries
      - SASA summary and atom-level exposure data
      - data completeness / sanity checks
"""

from __future__ import annotations

import argparse
import re
import sys
import textwrap
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# Column name synonym registry
# =============================================================================

RECRUITER_CODE_CANDIDATES = [
    "RECRUITER_CODE", "Recruiter_Code", "RecruiterCode", "Recruiter_Code_ID",
    "Recruiter_ID", "Recruiter", "Recruiter_Tag", "LR_Code", "LR_ID", "LR",
    "Unified_Recruiter_Code", "recruiter_code", "recruiter", "recruiter_id",
]

LIGAND_NAME_CANDIDATES = [
    "LIGAND", "Ligand", "ligand",
    "Ligand_Name", "LIGAND_NAME", "ligand_name",
    "Ligand_Code", "LIGAND_CODE", "ligand_code",
    "Original_Ligand", "ORIGINAL_LIGAND", "original_ligand",
    "Original_Ligand_Code", "ORIGINAL_LIGAND_CODE", "original_ligand_code",
    "PDB_Ligand", "PDB_LIGAND", "pdb_ligand",
    "HET", "HET_CODE", "HetCode", "Het_Code", "het_code",
    "Component", "Component_ID", "COMPONENT_ID", "component_id",
    "Chemical_Component", "Chemical_Component_ID",
]

COLUMN_SYNONYMS = {
    "ligase": [
        "Ligase", "E3", "E3_Ligase", "E3 Ligase", "Ligase_Name", "Protein",
        "Target", "Target_Ligase", "Recruiter_Ligase", "ligase", "e3_ligase"
    ],
    "recruiter_code": RECRUITER_CODE_CANDIDATES,
    "ligand_name": LIGAND_NAME_CANDIDATES,
    "instance": [
        "Instance_ID", "Ligand_Instance", "Recruiter_Instance", "PDB_Ligand_Instance",
        "Ligand_Instance_ID", "Instance", "Pose_ID", "Pose", "pose", "Conformation",
        "Conformer", "Conformer_ID", "Conformation_ID"
    ],
    "pdb": [
        "PDB", "PDB_ID", "PDBID", "PDB_Code", "Structure", "Structure_ID",
        "StructureID", "pdb", "pdb_id"
    ],
    "smiles": [
        "SMILES", "smiles", "Canonical_SMILES", "canonical_smiles",
        "Ligand_SMILES", "Recruiter_SMILES", "Mapped_SMILES", "Isomeric_SMILES",
        "Smiles", "LigandSmiles", "Scaffold_SMILES", "ScaffoldSmiles"
    ],
    "scaffold": [
        "Scaffold_ID", "Scaffold", "Murcko_Scaffold", "Murcko_Scaffold_ID",
        "Unified_Scaffold_ID", "Bemis_Murcko", "Bemis_Murcko_Scaffold",
        "scaffold_id"
    ],
    "supercluster": [
        "Supercluster_ID", "Supercluster", "Scaffold_Supercluster_ID",
        "Supercluster_Label", "Supercluster_Name", "Supercluster_Key"
    ],
    "atom": [
        "Atom", "Atom_Name", "Atom_ID", "Atom_Index", "atom", "atom_id",
        "Ligand_Atom", "LigandAtom"
    ],
    "sasa": [
        "SASA", "sasa", "Atom_SASA", "Ligand_SASA", "Total_SASA",
        "SASA_A2", "SASA_angstrom2", "SASA_Å2", "SASA_Value"
    ],
}

DESCRIPTOR_SYNONYMS = {
    "Molecular Weight": ["MW", "MolWt", "MolecularWeight", "Molecular_Weight", "ExactMolWt"],
    "LogP": ["LogP", "MolLogP", "Crippen_LogP", "cLogP"],
    "TPSA": ["TPSA", "TopoPSA", "Polar_Surface_Area"],
    "HBD": ["HBD", "NumHDonors", "H_Donors", "Hydrogen_Bond_Donors"],
    "HBA": ["HBA", "NumHAcceptors", "H_Acceptors", "Hydrogen_Bond_Acceptors"],
    "Rotatable Bonds": ["RotatableBonds", "NumRotatableBonds", "RotB", "Num_Rotatable_Bonds"],
    "Ring Count": ["RingCount", "NumRings", "Num_Rings"],
    "Fraction Csp3": ["FractionCSP3", "FracCSP3", "Fsp3", "Fraction_Csp3"],
    "QED": ["QED", "qed"],
    "Bertz Complexity": ["BertzCT", "Bertz", "Bertz_Index", "Molecular_Complexity", "Complexity"],
    "MR": ["MolMR", "Molecular_Refractivity", "MR"],
    "Lipinski Pass": ["Lipinski_Pass", "Lipinski", "RuleOfFive", "RO5_Pass", "Pass_Lipinski"],
}

PREFERRED_MASTER_MAP = "Recruiter_Master_Map.csv"

PREFERRED_LIGASE_SUMMARY_TABLES = [
    "Ligand_Instance_Recruiter_Codes.csv",
    "Ligase_Ligand_Metadata.csv",
    "Recruiter_Code_Crosswalk.csv",
    "Ligase_Recruiters_Scaffold.csv",
    "Recruiter_Master_Map.csv",
]

PREFERRED_SCAFFOLD_TABLES = [
    "Ligase_Recruiters_Scaffold.csv",
    "Ligase_Scaffold_Data.csv",
    "Ligase_Scaffold_Frequency.csv",
    "Ligase_Scaffold_Summary.csv",
    "Scaffold_Unified_Map.csv",
]

PREFERRED_DESCRIPTOR_TABLES = [
    "Ligase_Chemical_Descriptors.csv",
]

PREFERRED_SASA_SUMMARY_TABLES = [
    "Ligase_Ligand_SASA_summary.csv",
]

PREFERRED_SASA_ATOM_TABLES = [
    "Ligase_Ligand_SASA_atoms.csv",
]


# =============================================================================
# Formatting helpers
# =============================================================================

def line(char: str = "=", width: int = 100) -> str:
    return char * width


def section(title: str) -> str:
    return f"\n{line('=')}\n{title}\n{line('=')}\n"


def subsection(title: str) -> str:
    return f"\n{line('-')}\n{title}\n{line('-')}\n"


def clean_value(x):
    if pd.isna(x):
        return np.nan
    if isinstance(x, str):
        x = x.strip()
        if x == "" or x.lower() in {"nan", "none", "null", "na", "n/a"}:
            return np.nan
    return x


def normalize_text_series(s: pd.Series) -> pd.Series:
    return s.map(clean_value).astype("string").str.strip()


def normalize_code_series(s: pd.Series, uppercase: bool = True) -> pd.Series:
    out = normalize_text_series(s)
    out = out.str.replace(r"\s+", "", regex=True)
    if uppercase:
        out = out.str.upper()
    return out.replace("", pd.NA)


def normalize_ligand_series(s: pd.Series) -> pd.Series:
    """Normalize original PDB/component ligand names, not LR codes."""
    out = normalize_code_series(s, uppercase=True)
    return out.replace("", pd.NA)


def normalize_recruiter_code_series(s: pd.Series) -> pd.Series:
    """Normalize LR recruiter tags."""
    out = normalize_code_series(s, uppercase=True)
    return out.replace("", pd.NA)


def normalize_smiles_series(s: pd.Series) -> pd.Series:
    """Simple string normalization only; not chemistry-level canonicalization."""
    return (
        s.map(clean_value)
        .astype("string")
        .str.strip()
        .str.replace(r"\s+", "", regex=True)
        .replace("", pd.NA)
    )


def count_unique_nonempty(series: pd.Series) -> int:
    if series is None or len(series) == 0:
        return 0
    return int(series.dropna().astype("string").str.strip().replace("", pd.NA).dropna().nunique())


def safe_float(x, digits: int = 3) -> str:
    try:
        if pd.isna(x):
            return "NA"
        return f"{float(x):,.{digits}f}"
    except Exception:
        return "NA"


def format_df(df: pd.DataFrame, max_rows: int = 25, max_colwidth: int = 48) -> str:
    if df is None or df.empty:
        return "No data available.\n"

    out = df.copy()

    for col in out.columns:
        if out[col].dtype == object or str(out[col].dtype).startswith("string"):
            out[col] = out[col].astype(str).map(
                lambda x: x if len(x) <= max_colwidth else x[: max_colwidth - 3] + "..."
            )

    if len(out) > max_rows:
        out = out.head(max_rows)

    return out.to_string(index=False) + "\n"


def wrap_paragraph(text: str, width: int = 100) -> str:
    return "\n".join(textwrap.wrap(text, width=width)) + "\n"


# =============================================================================
# Data loading and column detection
# =============================================================================

def read_csv_safely(path: Path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path)
    except UnicodeDecodeError:
        try:
            return pd.read_csv(path, encoding="latin-1")
        except Exception as e:
            print(f"[WARN] Could not read {path}: {e}", file=sys.stderr)
            return None
    except Exception as e:
        print(f"[WARN] Could not read {path}: {e}", file=sys.stderr)
        return None


def load_tables(table_dir: Path) -> Dict[str, pd.DataFrame]:
    tables = {}
    if not table_dir.exists():
        raise FileNotFoundError(f"Table directory not found: {table_dir}")

    for path in sorted(table_dir.glob("*.csv")):
        df = read_csv_safely(path)
        if df is not None:
            tables[path.name] = df

    return tables


def find_column_exact_first(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    if df is None or df.empty:
        return None

    lower_map = {c.lower().strip(): c for c in df.columns}

    for cand in candidates:
        key = cand.lower().strip()
        if key in lower_map:
            return lower_map[key]

    return None


def find_column_relaxed(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    if df is None or df.empty:
        return None

    exact = find_column_exact_first(df, candidates)
    if exact:
        return exact

    # conservative relaxed matching: split on separators and require a meaningful token match
    for col in df.columns:
        col_l = col.lower().strip()
        tokens = set(re.split(r"[^a-z0-9]+", col_l))
        for cand in candidates:
            cand_l = cand.lower().strip()
            cand_tokens = set(re.split(r"[^a-z0-9]+", cand_l))
            if not cand_tokens:
                continue
            if cand_l in col_l or cand_tokens.issubset(tokens):
                return col

    return None


def column_lr_score(s: pd.Series) -> float:
    vals = normalize_recruiter_code_series(s).dropna().astype(str)
    if vals.empty:
        return 0.0
    sample = vals.head(500)
    return float(sample.str.match(r"^LR\d{3,}$", na=False).mean())


def column_ligand_code_score(s: pd.Series) -> float:
    vals = normalize_ligand_series(s).dropna().astype(str)
    if vals.empty:
        return 0.0
    sample = vals.head(500)
    # PDB chemical component IDs are commonly 1-5 alphanumeric chars.
    # Exclude LR recruiter codes and obvious SMILES-like strings.
    is_component_like = sample.str.match(r"^[A-Z0-9]{1,6}$", na=False)
    is_lr = sample.str.match(r"^LR\d{3,}$", na=False)
    has_smiles_chars = sample.str.contains(r"[=\[\]\(\)#@/\\]", regex=True, na=False)
    return float((is_component_like & ~is_lr & ~has_smiles_chars).mean())


def infer_recruiter_code_column(df: pd.DataFrame) -> Optional[str]:
    col = find_column_exact_first(df, RECRUITER_CODE_CANDIDATES)
    if col:
        return col

    best_col = None
    best_score = 0.0
    for c in df.columns:
        score = column_lr_score(df[c])
        if score > best_score:
            best_score = score
            best_col = c

    if best_score >= 0.70:
        return best_col

    return find_column_relaxed(df, RECRUITER_CODE_CANDIDATES)


def infer_ligand_name_column(df: pd.DataFrame, recruiter_col: Optional[str] = None) -> Optional[str]:
    # Strongly prefer exact ligand/component columns.
    col = find_column_exact_first(df, LIGAND_NAME_CANDIDATES)
    if col and col != recruiter_col:
        return col

    # Then score columns by whether they look like PDB component names and are not LR codes.
    best_col = None
    best_score = 0.0
    for c in df.columns:
        if c == recruiter_col:
            continue
        c_l = c.lower()
        if "smiles" in c_l or "atom" in c_l or "sasa" in c_l:
            continue
        score = column_ligand_code_score(df[c])
        if score > best_score:
            best_score = score
            best_col = c

    if best_score >= 0.70:
        return best_col

    relaxed = find_column_relaxed(df, LIGAND_NAME_CANDIDATES)
    if relaxed and relaxed != recruiter_col:
        return relaxed

    return None


def find_column(df: pd.DataFrame, kind: str) -> Optional[str]:
    if kind == "recruiter_code":
        return infer_recruiter_code_column(df)
    if kind == "ligand_name":
        return infer_ligand_name_column(df, recruiter_col=infer_recruiter_code_column(df))
    return find_column_relaxed(df, COLUMN_SYNONYMS.get(kind, []))


def find_descriptor_column(df: pd.DataFrame, descriptor_name: str) -> Optional[str]:
    return find_column_relaxed(df, DESCRIPTOR_SYNONYMS.get(descriptor_name, []))


def table_profile(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    rec = find_column(df, "recruiter_code")
    lig = infer_ligand_name_column(df, recruiter_col=rec)
    return {
        "ligase": find_column(df, "ligase"),
        "recruiter_code": rec,
        "ligand_name": lig,
        "instance": find_column(df, "instance"),
        "pdb": find_column(df, "pdb"),
        "smiles": find_column(df, "smiles"),
        "scaffold": find_column(df, "scaffold"),
        "supercluster": find_column(df, "supercluster"),
        "atom": find_column(df, "atom"),
        "sasa": find_column(df, "sasa"),
    }


def choose_table(
    tables: Dict[str, pd.DataFrame],
    preferred_names: List[str],
    required_kinds: List[str],
) -> Tuple[Optional[str], Optional[pd.DataFrame], Dict[str, Optional[str]]]:
    def ok(df):
        prof = table_profile(df)
        return all(prof.get(kind) is not None for kind in required_kinds), prof

    for name in preferred_names:
        if name in tables:
            good, prof = ok(tables[name])
            if good:
                return name, tables[name], prof

    for name, df in tables.items():
        good, prof = ok(df)
        if good:
            return name, df, prof

    return None, None, {}


def numeric_summary(series: pd.Series) -> Dict[str, float]:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return {}
    return {
        "n": int(s.count()),
        "mean": float(s.mean()),
        "median": float(s.median()),
        "std": float(s.std(ddof=1)) if len(s) > 1 else 0.0,
        "min": float(s.min()),
        "max": float(s.max()),
    }


# =============================================================================
# Master map handling
# =============================================================================

def build_master_mapping(tables: Dict[str, pd.DataFrame]) -> Tuple[pd.DataFrame, Dict[str, str], List[str]]:
    """
    Build the authoritative mapping:

        LR recruiter code -> original ligand/component name

    Returns:
        master_clean, columns_used, warnings
    """
    warnings: List[str] = []

    if PREFERRED_MASTER_MAP not in tables:
        raise FileNotFoundError(
            f"Required {PREFERRED_MASTER_MAP} not found in Ligase_Table. "
            "This script needs it to map LR recruiter codes back to original ligand names."
        )

    master = tables[PREFERRED_MASTER_MAP].copy()
    rec_col = infer_recruiter_code_column(master)
    lig_col = infer_ligand_name_column(master, recruiter_col=rec_col)

    if rec_col is None:
        raise ValueError(
            f"Could not identify recruiter-code column in {PREFERRED_MASTER_MAP}.\n"
            f"Columns found: {list(master.columns)}\n"
            "Expected something like RECRUITER_CODE, Recruiter_Code, Recruiter, LR_Code, or LR_ID."
        )

    if lig_col is None:
        raise ValueError(
            f"Could not identify original ligand/component column in {PREFERRED_MASTER_MAP}.\n"
            f"Columns found: {list(master.columns)}\n"
            "Expected something like LIGAND, Ligand, Ligand_Name, Ligand_Code, HET, or Component_ID.\n"
            "This is important because unique ligand counts must be based on the original ligand name, not LR codes."
        )

    clean = pd.DataFrame({
        "Recruiter_Code": normalize_recruiter_code_series(master[rec_col]),
        "Ligand": normalize_ligand_series(master[lig_col]),
    })

    # Carry through useful optional columns if present.
    for kind, out_name in [
        ("ligase", "Ligase"),
        ("pdb", "PDB"),
        ("smiles", "SMILES"),
        ("scaffold", "Scaffold_ID"),
        ("supercluster", "Supercluster_ID"),
    ]:
        col = find_column(master, kind)
        if col and col not in {rec_col, lig_col}:
            if kind == "ligase":
                clean[out_name] = normalize_text_series(master[col])
            elif kind == "pdb":
                clean[out_name] = normalize_code_series(master[col], uppercase=True)
            elif kind == "smiles":
                clean[out_name] = normalize_smiles_series(master[col])
            else:
                clean[out_name] = normalize_text_series(master[col])

    clean = clean.dropna(subset=["Recruiter_Code", "Ligand"]).drop_duplicates()

    # If a recruiter code maps to multiple ligands, flag it and keep the first for downstream joins.
    multi = clean.groupby("Recruiter_Code")["Ligand"].nunique().reset_index(name="Ligand_Count")
    bad = multi[multi["Ligand_Count"] > 1]
    if not bad.empty:
        warnings.append(
            f"WARNING: {len(bad):,} recruiter codes map to more than one ligand in {PREFERRED_MASTER_MAP}. "
            "This should be inspected. The report still uses all rows for ligand-level summaries."
        )

    columns_used = {
        "recruiter_code_column": rec_col,
        "ligand_name_column": lig_col,
    }

    return clean, columns_used, warnings


def add_ligand_names_from_master(df: pd.DataFrame, master_map: pd.DataFrame) -> pd.DataFrame:
    """Attach original Ligand names to any table that has Recruiter_Code/LR tags."""
    out = df.copy()
    prof = table_profile(out)

    if "Ligand" in out.columns:
        out["_Ligand"] = normalize_ligand_series(out["Ligand"])
    else:
        lig_col = prof.get("ligand_name")
        if lig_col:
            out["_Ligand"] = normalize_ligand_series(out[lig_col])
        else:
            out["_Ligand"] = pd.NA

    rec_col = prof.get("recruiter_code")
    if rec_col:
        out["_Recruiter_Code"] = normalize_recruiter_code_series(out[rec_col])
    else:
        out["_Recruiter_Code"] = pd.NA

    mapper = (
        master_map[["Recruiter_Code", "Ligand"]]
        .dropna()
        .drop_duplicates(subset=["Recruiter_Code"], keep="first")
    )

    out = out.merge(
        mapper.rename(columns={"Recruiter_Code": "_Recruiter_Code", "Ligand": "_Ligand_From_Master"}),
        how="left",
        on="_Recruiter_Code",
    )

    out["_Ligand"] = out["_Ligand"].fillna(out["_Ligand_From_Master"])
    out = out.drop(columns=["_Ligand_From_Master"], errors="ignore")

    return out


# =============================================================================
# Core summaries
# =============================================================================

def summarize_table_inventory(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for name, df in sorted(tables.items()):
        prof = table_profile(df)
        detected = [k for k, v in prof.items() if v is not None]
        rows.append({
            "Table": name,
            "Rows": len(df),
            "Columns": len(df.columns),
            "Detected key fields": ", ".join(detected) if detected else "none detected",
        })
    return pd.DataFrame(rows)


def summarize_folder_ligases(ligase_dir: Path) -> List[str]:
    if not ligase_dir.exists():
        return []
    dirs = [
        p.name for p in ligase_dir.iterdir()
        if p.is_dir() and not p.name.startswith(".") and p.name.lower() not in {"__pycache__", "csvcache"}
    ]
    return sorted(dirs)


def clean_ligase_values(series: pd.Series) -> pd.Series:
    vals = normalize_text_series(series).dropna().astype(str).str.strip()
    # Remove bogus numeric tokens like 1, 2, 3, 6 that can appear from matrix/table parsing.
    vals = vals[~vals.str.match(r"^\d+$", na=False)]
    vals = vals[~vals.str.lower().isin({"nan", "none", "null", "ligase", "e3"})]
    return vals.drop_duplicates()


def get_unique_values_from_all_tables(tables: Dict[str, pd.DataFrame], kind: str) -> pd.Series:
    values = []
    for _, df in tables.items():
        col = find_column(df, kind)
        if col is not None:
            if kind == "ligase":
                values.append(clean_ligase_values(df[col]))
            elif kind == "pdb":
                values.append(normalize_code_series(df[col], uppercase=True))
            elif kind == "ligand_name":
                values.append(normalize_ligand_series(df[col]))
            elif kind == "recruiter_code":
                values.append(normalize_recruiter_code_series(df[col]))
            else:
                values.append(normalize_text_series(df[col]))

    if not values:
        return pd.Series(dtype="string")

    return pd.concat(values, ignore_index=True).dropna().drop_duplicates()


def summarize_ligand_conformations(master_map: pd.DataFrame, tables: Dict[str, pd.DataFrame]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate ligand-level conformation summary.

    Conformation count = number of unique LR recruiter codes associated with that ligand.
    This intentionally preserves multiple crystal/conformation instances of the same ligand.
    """
    base = master_map.copy()

    # Add ligase/PDB/scaffold information from richer tables when possible.
    enrichment_sources = [
        "Ligand_Instance_Recruiter_Codes.csv",
        "Recruiter_Code_Crosswalk.csv",
        "Ligase_Ligand_Metadata.csv",
        "Ligase_Recruiters_Scaffold.csv",
        "Ligase_Recruiters_Superclustered.csv",
        "Ligase_Ligands_Smiles_3DMapped.csv",
    ]

    enriched_parts = [base]
    for name in enrichment_sources:
        if name not in tables:
            continue
        df = add_ligand_names_from_master(tables[name], base)
        prof = table_profile(df)
        rec_col = prof.get("recruiter_code")
        if rec_col is None and "_Recruiter_Code" not in df.columns:
            continue

        part = pd.DataFrame()
        part["Recruiter_Code"] = df["_Recruiter_Code"] if "_Recruiter_Code" in df.columns else normalize_recruiter_code_series(df[rec_col])
        part["Ligand"] = df["_Ligand"]

        ligase_col = prof.get("ligase")
        pdb_col = prof.get("pdb")
        smiles_col = prof.get("smiles")
        scaffold_col = prof.get("scaffold")
        supercluster_col = prof.get("supercluster")

        if ligase_col:
            part["Ligase"] = normalize_text_series(df[ligase_col])
        if pdb_col:
            part["PDB"] = normalize_code_series(df[pdb_col], uppercase=True)
        if smiles_col:
            part["SMILES"] = normalize_smiles_series(df[smiles_col])
        if scaffold_col:
            part["Scaffold_ID"] = normalize_text_series(df[scaffold_col])
        if supercluster_col:
            part["Supercluster_ID"] = normalize_text_series(df[supercluster_col])

        enriched_parts.append(part.dropna(subset=["Recruiter_Code", "Ligand"]))

    full = pd.concat(enriched_parts, ignore_index=True, sort=False).drop_duplicates()

    # Clean any bogus ligase values.
    if "Ligase" in full.columns:
        full["Ligase"] = normalize_text_series(full["Ligase"])
        full.loc[full["Ligase"].astype("string").str.match(r"^\d+$", na=False), "Ligase"] = pd.NA

    def join_examples(x, max_items=12):
        vals = [str(v) for v in pd.Series(x).dropna().drop_duplicates().tolist() if str(v).strip()]
        vals = sorted(vals)
        if len(vals) > max_items:
            return ", ".join(vals[:max_items]) + f", ... (+{len(vals) - max_items} more)"
        return ", ".join(vals)

    grouped = (
        full.dropna(subset=["Ligand"])
        .groupby("Ligand", dropna=True)
        .agg(
            Conformations_LR_Codes=("Recruiter_Code", lambda x: count_unique_nonempty(x)),
            LR_Code_List=("Recruiter_Code", join_examples),
            Unique_Ligases=("Ligase", lambda x: count_unique_nonempty(x) if "Ligase" in full.columns else 0),
            Ligase_List=("Ligase", join_examples if "Ligase" in full.columns else lambda x: ""),
            Unique_PDBs=("PDB", lambda x: count_unique_nonempty(x) if "PDB" in full.columns else 0),
            PDB_List=("PDB", join_examples if "PDB" in full.columns else lambda x: ""),
            Unique_SMILES=("SMILES", lambda x: count_unique_nonempty(x) if "SMILES" in full.columns else 0),
            Unique_Scaffolds=("Scaffold_ID", lambda x: count_unique_nonempty(x) if "Scaffold_ID" in full.columns else 0),
            Scaffold_List=("Scaffold_ID", join_examples if "Scaffold_ID" in full.columns else lambda x: ""),
        )
        .reset_index()
        .sort_values(["Conformations_LR_Codes", "Unique_PDBs", "Unique_Ligases", "Ligand"], ascending=[False, False, False, True])
    )

    return grouped, full


def summarize_dataset_overview(
    tables: Dict[str, pd.DataFrame],
    ligase_dir: Path,
    master_map: pd.DataFrame,
    ligand_conf_summary: pd.DataFrame,
) -> Dict[str, object]:
    csv_ligases = get_unique_values_from_all_tables(tables, "ligase")
    ligase_folders = summarize_folder_ligases(ligase_dir)

    all_pdbs = get_unique_values_from_all_tables(tables, "pdb")
    all_scaffolds = get_unique_values_from_all_tables(tables, "scaffold")
    all_superclusters = get_unique_values_from_all_tables(tables, "supercluster")

    unique_ligands = count_unique_nonempty(master_map["Ligand"])
    unique_lr_codes = count_unique_nonempty(master_map["Recruiter_Code"])

    multi_conf_ligands = int((ligand_conf_summary["Conformations_LR_Codes"] > 1).sum()) if not ligand_conf_summary.empty else 0
    max_conf = int(ligand_conf_summary["Conformations_LR_Codes"].max()) if not ligand_conf_summary.empty else 0

    return {
        "csv_ligase_count_cleaned": count_unique_nonempty(csv_ligases),
        "csv_ligases_cleaned": sorted(csv_ligases.dropna().astype(str).tolist()),
        "ligase_folder_count": len(ligase_folders),
        "ligase_folders": ligase_folders,
        "unique_ligand_count": unique_ligands,
        "unique_ligand_basis": f"original ligand/component names from {PREFERRED_MASTER_MAP}:Ligand",
        "lr_recruiter_code_count": unique_lr_codes,
        "lr_recruiter_code_basis": f"unique LR recruiter codes from {PREFERRED_MASTER_MAP}:Recruiter_Code",
        "ligand_conformation_count": unique_lr_codes,
        "ligands_with_multiple_conformations": multi_conf_ligands,
        "max_conformations_for_single_ligand": max_conf,
        "unique_pdb_count": count_unique_nonempty(all_pdbs),
        "unique_scaffold_count": count_unique_nonempty(all_scaffolds),
        "unique_supercluster_count": count_unique_nonempty(all_superclusters),
    }


def summarize_ligase_level(tables: Dict[str, pd.DataFrame], master_map: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
    source_name, df, prof = choose_table(
        tables,
        PREFERRED_LIGASE_SUMMARY_TABLES,
        required_kinds=["ligase"],
    )

    if df is None:
        return pd.DataFrame(), "No ligase-level source table detected."

    temp = add_ligand_names_from_master(df, master_map)
    prof = table_profile(temp)

    ligase_col = prof.get("ligase")
    rec_col = prof.get("recruiter_code")
    pdb_col = prof.get("pdb")
    scaffold_col = prof.get("scaffold")
    supercluster_col = prof.get("supercluster")

    temp["_Ligase"] = normalize_text_series(temp[ligase_col])
    temp.loc[temp["_Ligase"].astype("string").str.match(r"^\d+$", na=False), "_Ligase"] = pd.NA

    if rec_col:
        temp["_Recruiter_Code"] = normalize_recruiter_code_series(temp[rec_col])

    temp["_PDB"] = normalize_code_series(temp[pdb_col], uppercase=True) if pdb_col else pd.NA
    temp["_Scaffold"] = normalize_text_series(temp[scaffold_col]) if scaffold_col else pd.NA
    temp["_Supercluster"] = normalize_text_series(temp[supercluster_col]) if supercluster_col else pd.NA

    grouped = (
        temp.dropna(subset=["_Ligase"])
        .groupby("_Ligase", dropna=True)
        .agg(
            Observations=("_Ligase", "size"),
            Unique_Ligands=("_Ligand", lambda x: count_unique_nonempty(x)),
            LR_Conformation_Codes=("_Recruiter_Code", lambda x: count_unique_nonempty(x)),
            Unique_PDBs=("_PDB", lambda x: count_unique_nonempty(x)),
            Unique_Scaffolds=("_Scaffold", lambda x: count_unique_nonempty(x)),
            Unique_Superclusters=("_Supercluster", lambda x: count_unique_nonempty(x)),
        )
        .reset_index()
        .rename(columns={"_Ligase": "Ligase"})
    )

    grouped = grouped.sort_values(
        ["Unique_Ligands", "LR_Conformation_Codes", "Unique_Scaffolds", "Observations"],
        ascending=[False, False, False, False],
    )

    note = (
        f"Ligase-level summary source: {source_name}. "
        f"Unique ligands are counted from original ligand/component names mapped through {PREFERRED_MASTER_MAP}; "
        "LR recruiter codes are counted separately as conformation/instance tags."
    )

    return grouped, note


def summarize_top_ligands(ligand_conf_summary: pd.DataFrame, top_n: int) -> Tuple[pd.DataFrame, str]:
    if ligand_conf_summary is None or ligand_conf_summary.empty:
        return pd.DataFrame(), "No ligand conformation summary available."

    cols = [
        "Ligand", "Conformations_LR_Codes", "Unique_Ligases", "Unique_PDBs",
        "Unique_SMILES", "Unique_Scaffolds", "LR_Code_List", "Ligase_List", "PDB_List"
    ]
    cols = [c for c in cols if c in ligand_conf_summary.columns]

    out = ligand_conf_summary[cols].head(top_n).copy()
    note = (
        f"Top ligand table uses original ligand/component names from {PREFERRED_MASTER_MAP}. "
        "Conformations_LR_Codes is the number of LR recruiter codes associated with each ligand."
    )
    return out, note


def summarize_scaffolds(tables: Dict[str, pd.DataFrame], master_map: pd.DataFrame, top_n: int) -> Tuple[pd.DataFrame, Dict[str, object], str]:
    source_name, df, prof = choose_table(
        tables,
        PREFERRED_SCAFFOLD_TABLES,
        required_kinds=["scaffold"],
    )

    if df is None:
        return pd.DataFrame(), {}, "No scaffold source table detected."

    temp = add_ligand_names_from_master(df, master_map)
    prof = table_profile(temp)

    scaffold_col = prof.get("scaffold")
    ligase_col = prof.get("ligase")
    rec_col = prof.get("recruiter_code")
    pdb_col = prof.get("pdb")
    supercluster_col = prof.get("supercluster")

    temp["_Scaffold"] = normalize_text_series(temp[scaffold_col])
    temp["_Ligase"] = normalize_text_series(temp[ligase_col]) if ligase_col else pd.NA
    if ligase_col:
        temp.loc[temp["_Ligase"].astype("string").str.match(r"^\d+$", na=False), "_Ligase"] = pd.NA
    temp["_Recruiter_Code"] = normalize_recruiter_code_series(temp[rec_col]) if rec_col else temp.get("_Recruiter_Code", pd.NA)
    temp["_PDB"] = normalize_code_series(temp[pdb_col], uppercase=True) if pdb_col else pd.NA
    temp["_Supercluster"] = normalize_text_series(temp[supercluster_col]) if supercluster_col else pd.NA

    grouped = (
        temp.dropna(subset=["_Scaffold"])
        .groupby("_Scaffold", dropna=True)
        .agg(
            Observations=("_Scaffold", "size"),
            Unique_Ligands=("_Ligand", lambda x: count_unique_nonempty(x)),
            LR_Conformation_Codes=("_Recruiter_Code", lambda x: count_unique_nonempty(x)),
            Unique_Ligases=("_Ligase", lambda x: count_unique_nonempty(x)),
            Unique_PDBs=("_PDB", lambda x: count_unique_nonempty(x)),
            Unique_Superclusters=("_Supercluster", lambda x: count_unique_nonempty(x)),
        )
        .reset_index()
        .rename(columns={"_Scaffold": "Scaffold"})
        .sort_values(["Unique_Ligands", "LR_Conformation_Codes", "Observations", "Unique_Ligases"], ascending=False)
        .head(top_n)
    )

    stats = {
        "source_table": source_name,
        "unique_scaffolds_in_source": count_unique_nonempty(temp["_Scaffold"]),
        "scaffold_column": scaffold_col,
    }

    note = (
        f"Scaffold summary source: {source_name}; scaffold column: {scaffold_col}. "
        "Unique_Ligands is based on original ligand/component names, not LR recruiter codes."
    )

    return grouped, stats, note


def summarize_superclusters(tables: Dict[str, pd.DataFrame], master_map: pd.DataFrame, top_n: int) -> Tuple[pd.DataFrame, str]:
    source_name, df, prof = choose_table(
        tables,
        [
            "Ligase_Scaffold_Superclusters.csv",
            "Ligase_Recruiters_Superclustered.csv",
            "Ligase_Scaffold_Supercluster_Frequency.csv",
            "Ligase_Scaffold_Supercluster_Matrix.csv",
        ],
        required_kinds=["supercluster"],
    )

    if df is None:
        return pd.DataFrame(), "No supercluster source table detected."

    temp = add_ligand_names_from_master(df, master_map)
    prof = table_profile(temp)

    supercluster_col = prof.get("supercluster")
    scaffold_col = prof.get("scaffold")
    ligase_col = prof.get("ligase")
    rec_col = prof.get("recruiter_code")

    temp["_Supercluster"] = normalize_text_series(temp[supercluster_col])
    temp["_Scaffold"] = normalize_text_series(temp[scaffold_col]) if scaffold_col else pd.NA
    temp["_Ligase"] = normalize_text_series(temp[ligase_col]) if ligase_col else pd.NA
    if ligase_col:
        temp.loc[temp["_Ligase"].astype("string").str.match(r"^\d+$", na=False), "_Ligase"] = pd.NA
    temp["_Recruiter_Code"] = normalize_recruiter_code_series(temp[rec_col]) if rec_col else temp.get("_Recruiter_Code", pd.NA)

    grouped = (
        temp.dropna(subset=["_Supercluster"])
        .groupby("_Supercluster", dropna=True)
        .agg(
            Observations=("_Supercluster", "size"),
            Unique_Ligands=("_Ligand", lambda x: count_unique_nonempty(x)),
            LR_Conformation_Codes=("_Recruiter_Code", lambda x: count_unique_nonempty(x)),
            Unique_Scaffolds=("_Scaffold", lambda x: count_unique_nonempty(x)),
            Unique_Ligases=("_Ligase", lambda x: count_unique_nonempty(x)),
        )
        .reset_index()
        .rename(columns={"_Supercluster": "Supercluster"})
        .sort_values(["Unique_Ligands", "Unique_Scaffolds", "Observations"], ascending=False)
        .head(top_n)
    )

    note = f"Supercluster summary source: {source_name}; supercluster column: {supercluster_col}."
    return grouped, note


def summarize_descriptors(tables: Dict[str, pd.DataFrame], master_map: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    source_name, df, prof = choose_table(
        tables,
        PREFERRED_DESCRIPTOR_TABLES,
        required_kinds=[],
    )

    if df is None:
        return pd.DataFrame(), pd.DataFrame(), "No descriptor table detected."

    temp = add_ligand_names_from_master(df, master_map)

    rows = []
    for descriptor in DESCRIPTOR_SYNONYMS.keys():
        col = find_descriptor_column(temp, descriptor)
        if col is None:
            continue

        if descriptor == "Lipinski Pass":
            s = temp[col].dropna()
            if s.empty:
                continue

            normalized = s.astype(str).str.strip().str.lower()
            pass_values = {"true", "yes", "y", "pass", "passes", "1", "1.0"}
            fail_values = {"false", "no", "n", "fail", "fails", "0", "0.0"}
            passes = int(normalized.isin(pass_values).sum())
            fails = int(normalized.isin(fail_values).sum())
            total = passes + fails

            if total > 0:
                rows.append({
                    "Descriptor": descriptor,
                    "Column": col,
                    "N": total,
                    "Mean": np.nan,
                    "Median": np.nan,
                    "Std": np.nan,
                    "Min": np.nan,
                    "Max": np.nan,
                    "Extra": f"{passes:,}/{total:,} pass ({passes / total * 100:.1f}%)",
                })
            continue

        stats = numeric_summary(temp[col])
        if not stats:
            continue

        rows.append({
            "Descriptor": descriptor,
            "Column": col,
            "N": stats["n"],
            "Mean": stats["mean"],
            "Median": stats["median"],
            "Std": stats["std"],
            "Min": stats["min"],
            "Max": stats["max"],
            "Extra": "",
        })

    summary = pd.DataFrame(rows)
    if not summary.empty:
        display_summary = summary.copy()
        for c in ["Mean", "Median", "Std", "Min", "Max"]:
            display_summary[c] = display_summary[c].map(lambda x: safe_float(x, 3))
    else:
        display_summary = pd.DataFrame()

    prof = table_profile(temp)
    ligase_col = prof.get("ligase")
    by_ligase = pd.DataFrame()

    if ligase_col:
        temp["_Ligase"] = normalize_text_series(temp[ligase_col])
        temp.loc[temp["_Ligase"].astype("string").str.match(r"^\d+$", na=False), "_Ligase"] = pd.NA

        selected_cols = {}
        for descriptor in ["Molecular Weight", "LogP", "TPSA", "QED", "Bertz Complexity"]:
            col = find_descriptor_column(temp, descriptor)
            if col:
                selected_cols[descriptor] = col

        rows2 = []
        for ligase, sub in temp.dropna(subset=["_Ligase"]).groupby("_Ligase"):
            entry = {
                "Ligase": ligase,
                "Rows": len(sub),
                "Unique_Ligands": count_unique_nonempty(sub["_Ligand"]),
                "LR_Conformation_Codes": count_unique_nonempty(sub["_Recruiter_Code"]),
            }
            for desc, col in selected_cols.items():
                vals = pd.to_numeric(sub[col], errors="coerce").dropna()
                if not vals.empty:
                    entry[f"{desc} median"] = vals.median()
                    entry[f"{desc} mean"] = vals.mean()
            rows2.append(entry)

        by_ligase = pd.DataFrame(rows2)
        if not by_ligase.empty:
            by_ligase = by_ligase.sort_values("Unique_Ligands", ascending=False)
            for col in by_ligase.columns:
                if col not in {"Ligase", "Rows", "Unique_Ligands", "LR_Conformation_Codes"}:
                    by_ligase[col] = by_ligase[col].map(lambda x: safe_float(x, 3))

    note = (
        f"Descriptor summary source: {source_name}. Descriptor rows are usually LR recruiter-code/conformation rows. "
        "The by-ligase table additionally reports unique original ligand names."
    )

    return display_summary, by_ligase, note


def summarize_sasa(tables: Dict[str, pd.DataFrame], master_map: pd.DataFrame, top_n: int) -> Tuple[str, pd.DataFrame, pd.DataFrame]:
    summary_text = []

    summary_name, summary_df, summary_prof = choose_table(
        tables,
        PREFERRED_SASA_SUMMARY_TABLES,
        required_kinds=[],
    )

    atom_name, atom_df, atom_prof = choose_table(
        tables,
        PREFERRED_SASA_ATOM_TABLES,
        required_kinds=[],
    )

    ligand_sasa_summary = pd.DataFrame()
    atom_sasa_top = pd.DataFrame()

    if summary_df is not None:
        temp = add_ligand_names_from_master(summary_df, master_map)
        prof = table_profile(temp)
        sasa_col = prof.get("sasa")
        ligase_col = prof.get("ligase")
        pdb_col = prof.get("pdb")

        summary_text.append(f"Ligand-level SASA source: {summary_name}")

        if sasa_col:
            stats = numeric_summary(temp[sasa_col])
            if stats:
                summary_text.append(
                    "Ligand-level SASA rows/conformations: "
                    f"N={stats['n']:,}, mean={stats['mean']:.3f}, median={stats['median']:.3f}, "
                    f"min={stats['min']:.3f}, max={stats['max']:.3f}"
                )

        if ligase_col:
            temp["_Ligase"] = normalize_text_series(temp[ligase_col])
            temp.loc[temp["_Ligase"].astype("string").str.match(r"^\d+$", na=False), "_Ligase"] = pd.NA
        else:
            temp["_Ligase"] = pd.NA

        temp["_PDB"] = normalize_code_series(temp[pdb_col], uppercase=True) if pdb_col else pd.NA

        if sasa_col:
            temp["_SASA"] = pd.to_numeric(temp[sasa_col], errors="coerce")
            ligand_sasa_summary = (
                temp.dropna(subset=["_Ligase"])
                .groupby("_Ligase")
                .agg(
                    Rows=("_Ligase", "size"),
                    Unique_Ligands=("_Ligand", lambda x: count_unique_nonempty(x)),
                    LR_Conformation_Codes=("_Recruiter_Code", lambda x: count_unique_nonempty(x)),
                    Unique_PDBs=("_PDB", lambda x: count_unique_nonempty(x)),
                    Mean_SASA=("_SASA", "mean"),
                    Median_SASA=("_SASA", "median"),
                    Max_SASA=("_SASA", "max"),
                )
                .reset_index()
                .rename(columns={"_Ligase": "Ligase"})
                .sort_values("Rows", ascending=False)
            )

            for col in ["Mean_SASA", "Median_SASA", "Max_SASA"]:
                ligand_sasa_summary[col] = ligand_sasa_summary[col].map(lambda x: safe_float(x, 3))
    else:
        summary_text.append("Ligand-level SASA source: not detected")

    if atom_df is not None:
        temp = add_ligand_names_from_master(atom_df, master_map)
        prof = table_profile(temp)
        sasa_col = prof.get("sasa")
        ligase_col = prof.get("ligase")
        pdb_col = prof.get("pdb")
        atom_col = prof.get("atom")

        summary_text.append(f"Atom-level SASA source: {atom_name}")

        if sasa_col:
            stats = numeric_summary(temp[sasa_col])
            if stats:
                summary_text.append(
                    "Atom-level SASA rows: "
                    f"N={stats['n']:,}, mean={stats['mean']:.3f}, median={stats['median']:.3f}, "
                    f"min={stats['min']:.3f}, max={stats['max']:.3f}"
                )

            temp["_SASA"] = pd.to_numeric(temp[sasa_col], errors="coerce")
            temp["_Ligase"] = normalize_text_series(temp[ligase_col]) if ligase_col else pd.NA
            if ligase_col:
                temp.loc[temp["_Ligase"].astype("string").str.match(r"^\d+$", na=False), "_Ligase"] = pd.NA
            temp["_PDB"] = normalize_code_series(temp[pdb_col], uppercase=True) if pdb_col else pd.NA
            temp["_Atom"] = normalize_text_series(temp[atom_col]) if atom_col else pd.NA

            atom_sasa_top = (
                temp.dropna(subset=["_SASA"])
                .sort_values("_SASA", ascending=False)
                .head(top_n)
                [["_Ligase", "_Ligand", "_Recruiter_Code", "_PDB", "_Atom", "_SASA"]]
                .rename(columns={
                    "_Ligase": "Ligase",
                    "_Ligand": "Ligand",
                    "_Recruiter_Code": "LR_Code",
                    "_PDB": "PDB",
                    "_Atom": "Atom",
                    "_SASA": "SASA",
                })
            )
            atom_sasa_top["SASA"] = atom_sasa_top["SASA"].map(lambda x: safe_float(x, 3))
    else:
        summary_text.append("Atom-level SASA source: not detected")

    return "\n".join(summary_text), ligand_sasa_summary, atom_sasa_top


def data_completeness_checks(tables: Dict[str, pd.DataFrame], master_map: pd.DataFrame) -> str:
    lines = []

    checks = [
        ("Recruiter master map", "Recruiter_Master_Map.csv"),
        ("Recruiter SMILES map", "Recruiter_SMILES_Map.csv"),
        ("Ligand/recruiter instance map", "Ligand_Instance_Recruiter_Codes.csv"),
        ("Ligase/recruiter/scaffold map", "Ligase_Recruiters_Scaffold.csv"),
        ("Chemical descriptors", "Ligase_Chemical_Descriptors.csv"),
        ("SASA ligand summary", "Ligase_Ligand_SASA_summary.csv"),
        ("SASA atom table", "Ligase_Ligand_SASA_atoms.csv"),
        ("Duplicate ligand table", "Ligase_Duplicate_Ligands.csv"),
        ("Unified scaffold map", "Scaffold_Unified_Map.csv"),
    ]

    for label, name in checks:
        if name in tables:
            df = tables[name]
            lines.append(f"[OK] {label}: {name} present with {len(df):,} rows and {len(df.columns):,} columns.")
        else:
            lines.append(f"[MISSING] {label}: {name} not found.")

    lines.append("")

    missing_lig = int(master_map["Ligand"].isna().sum())
    missing_rec = int(master_map["Recruiter_Code"].isna().sum())
    lines.append(f"[CHECK] Master map rows after cleaning: {len(master_map):,}")
    lines.append(f"[CHECK] Missing ligand names in cleaned master map: {missing_lig:,}")
    lines.append(f"[CHECK] Missing LR recruiter codes in cleaned master map: {missing_rec:,}")

    duplicated_lr = master_map[master_map.duplicated("Recruiter_Code", keep=False)].sort_values("Recruiter_Code")
    if not duplicated_lr.empty:
        lines.append(f"[CHECK] Recruiter codes appearing in multiple rows: {duplicated_lr['Recruiter_Code'].nunique():,}")
    else:
        lines.append("[OK] Each LR recruiter code appears once in the cleaned master map.")

    lines.append("")

    for name, df in tables.items():
        prof = table_profile(df)
        missing_notes = []

        for kind in ["ligase", "recruiter_code", "ligand_name", "smiles", "pdb", "scaffold"]:
            col = prof.get(kind)
            if col is not None:
                n_missing = int(normalize_text_series(df[col]).isna().sum())
                pct = (n_missing / len(df) * 100) if len(df) else 0
                if n_missing > 0:
                    missing_notes.append(f"{kind}:{col} missing {n_missing:,}/{len(df):,} ({pct:.1f}%)")

        if missing_notes:
            lines.append(f"[CHECK] {name}: " + "; ".join(missing_notes))

    return "\n".join(lines) + "\n"


# =============================================================================
# Report writer
# =============================================================================

def build_report(
    tables: Dict[str, pd.DataFrame],
    table_dir: Path,
    ligase_dir: Path,
    top_n: int,
    write_derived: bool = True,
) -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    master_map, master_columns, master_warnings = build_master_mapping(tables)
    ligand_conf_summary, enriched_recruiter_map = summarize_ligand_conformations(master_map, tables)
    overview = summarize_dataset_overview(tables, ligase_dir, master_map, ligand_conf_summary)

    if write_derived:
        out1 = table_dir / "Derived_Ligand_Conformation_Summary.csv"
        out2 = table_dir / "Derived_RecruiterCode_To_Ligand_Map.csv"
        ligand_conf_summary.to_csv(out1, index=False)
        enriched_recruiter_map.to_csv(out2, index=False)

    report = []

    report.append(line("="))
    report.append("E3 RECRUITER LIGANDALYZER DATA SUMMARY")
    report.append(line("="))
    report.append(f"Generated: {now}")
    report.append(f"Working directory: {Path.cwd()}")
    report.append(f"Table directory: {table_dir.resolve()}")
    report.append(f"Ligase directory: {ligase_dir.resolve()}")
    report.append("")

    report.append(wrap_paragraph(
        "This report summarizes the current E3 Recruiter Ligandalyzer data tables for manuscript, "
        "platform, and quality-control use. In this fixed version, unique ligands are counted using "
        "the original ligand/component names from Recruiter_Master_Map.csv. LR recruiter codes are "
        "reported separately as recruiter-code conformation/instance tags. This avoids the misleading "
        "interpretation that every LR##### code is a chemically unique ligand."
    ))

    report.append(section("0. MASTER MAP INTERPRETATION"))
    report.append(format_df(pd.DataFrame([
        ("Master mapping source", PREFERRED_MASTER_MAP),
        ("Recruiter-code column detected", master_columns.get("recruiter_code_column")),
        ("Ligand/component column detected", master_columns.get("ligand_name_column")),
        ("Unique original ligand/component names", f"{overview['unique_ligand_count']:,}"),
        ("Unique LR recruiter-code tags", f"{overview['lr_recruiter_code_count']:,}"),
        ("Ligands with >1 LR code/conformation", f"{overview['ligands_with_multiple_conformations']:,}"),
        ("Maximum LR codes for one ligand", f"{overview['max_conformations_for_single_ligand']:,}"),
    ], columns=["Metric", "Value"]), max_rows=50))

    if master_warnings:
        report.append(subsection("Master map warnings"))
        for w in master_warnings:
            report.append(w + "\n")

    report.append(section("1. MANUSCRIPT-FACING DATASET OVERVIEW"))
    overview_rows = [
        ("CSV tables loaded", f"{len(tables):,}"),
        ("Curated E3 ligase folders detected", f"{overview['ligase_folder_count']:,}"),
        ("Clean E3 ligases detected from CSVs", f"{overview['csv_ligase_count_cleaned']:,}"),
        ("Unique ligands / chemical component names", f"{overview['unique_ligand_count']:,}"),
        ("Unique ligand basis", str(overview["unique_ligand_basis"])),
        ("LR recruiter-code conformation tags", f"{overview['lr_recruiter_code_count']:,}"),
        ("LR recruiter-code basis", str(overview["lr_recruiter_code_basis"])),
        ("Ligands with multiple conformations", f"{overview['ligands_with_multiple_conformations']:,}"),
        ("Unique PDB / structure IDs detected", f"{overview['unique_pdb_count']:,}"),
        ("Unique scaffolds detected", f"{overview['unique_scaffold_count']:,}"),
        ("Unique scaffold superclusters detected", f"{overview['unique_supercluster_count']:,}"),
    ]
    report.append(format_df(pd.DataFrame(overview_rows, columns=["Metric", "Value"]), max_rows=100))

    report.append(subsection("Curated ligase folders detected"))
    report.append(", ".join(overview["ligase_folders"]) + "\n" if overview["ligase_folders"] else "No ligase folders detected.\n")

    report.append(subsection("Clean ligases detected from CSV tables"))
    report.append(", ".join(overview["csv_ligases_cleaned"]) + "\n" if overview["csv_ligases_cleaned"] else "No clean ligase names detected from CSV columns.\n")

    report.append(section("2. CSV TABLE INVENTORY"))
    report.append(format_df(summarize_table_inventory(tables), max_rows=250, max_colwidth=64))

    report.append(section("3. LIGAND CONFORMATION SUMMARY"))
    top_ligands, top_ligand_note = summarize_top_ligands(ligand_conf_summary, top_n=top_n)
    report.append(top_ligand_note + "\n\n")
    report.append(format_df(top_ligands, max_rows=top_n, max_colwidth=72))

    report.append(subsection("Conformation count distribution"))
    if not ligand_conf_summary.empty:
        dist = (
            ligand_conf_summary["Conformations_LR_Codes"]
            .value_counts()
            .sort_index()
            .reset_index()
        )
        dist.columns = ["LR codes / conformations per ligand", "Number of ligands"]
        report.append(format_df(dist, max_rows=100))
    else:
        report.append("No conformation distribution available.\n")

    ligase_summary, ligase_note = summarize_ligase_level(tables, master_map)
    report.append(section("4. LIGASE-LEVEL LIGAND, CONFORMATION, AND SCAFFOLD SUMMARY"))
    report.append(ligase_note + "\n\n")
    report.append(format_df(ligase_summary, max_rows=100, max_colwidth=56))

    if not ligase_summary.empty:
        report.append(subsection("Manuscript-style interpretation"))
        top_names = ", ".join(ligase_summary.head(5)["Ligase"].astype(str).tolist())
        report.append(wrap_paragraph(
            f"The most heavily represented ligases by unique original ligand/component count are: {top_names}. "
            "LR recruiter-code counts are reported separately because they capture distinct curated instances "
            "or conformations of ligands rather than necessarily distinct chemical components."
        ))

    scaffold_top, scaffold_stats, scaffold_note = summarize_scaffolds(tables, master_map, top_n=top_n)
    report.append(section("5. SCAFFOLD SUMMARY"))
    report.append(scaffold_note + "\n\n")
    if scaffold_stats:
        report.append(format_df(pd.DataFrame([
            ("Scaffold source table", scaffold_stats.get("source_table")),
            ("Scaffold column", scaffold_stats.get("scaffold_column")),
            ("Unique scaffolds in source", scaffold_stats.get("unique_scaffolds_in_source")),
        ], columns=["Metric", "Value"]), max_rows=20))
    report.append(subsection(f"Top {top_n} scaffolds by ligand/conformation coverage"))
    report.append(format_df(scaffold_top, max_rows=top_n, max_colwidth=72))

    supercluster_top, supercluster_note = summarize_superclusters(tables, master_map, top_n=top_n)
    report.append(section("6. SCAFFOLD SUPERCLUSTER SUMMARY"))
    report.append(supercluster_note + "\n\n")
    report.append(format_df(supercluster_top, max_rows=top_n, max_colwidth=72))

    descriptor_summary, descriptor_by_ligase, descriptor_note = summarize_descriptors(tables, master_map)
    report.append(section("7. PHYSICOCHEMICAL AND TOPOLOGICAL DESCRIPTOR SUMMARY"))
    report.append(descriptor_note + "\n\n")
    report.append(format_df(descriptor_summary, max_rows=100, max_colwidth=56))
    report.append(subsection("Descriptor medians/means by ligase"))
    report.append(format_df(descriptor_by_ligase, max_rows=100, max_colwidth=56))

    sasa_text, ligand_sasa_summary, atom_sasa_top = summarize_sasa(tables, master_map, top_n=top_n)
    report.append(section("8. SOLVENT ACCESSIBLE SURFACE AREA / SASA SUMMARY"))
    report.append(sasa_text + "\n")
    report.append(subsection("Ligand-level SASA by ligase"))
    report.append(format_df(ligand_sasa_summary, max_rows=100, max_colwidth=56))
    report.append(subsection(f"Top {top_n} most solvent-exposed atoms by SASA"))
    report.append(format_df(atom_sasa_top, max_rows=top_n, max_colwidth=72))

    report.append(section("9. DATA COMPLETENESS AND SANITY CHECKS"))
    report.append(data_completeness_checks(tables, master_map))

    report.append(section("10. DRAFTABLE SUMMARY LANGUAGE"))
    report.append(wrap_paragraph(
        f"The current E3 Recruiter Ligandalyzer dataset contains {overview['unique_ligand_count']:,} "
        f"unique ligand/component names represented by {overview['lr_recruiter_code_count']:,} curated "
        f"LR recruiter-code conformation tags across {overview['unique_pdb_count']:,} unique PDB/structure "
        f"identifiers and {overview['ligase_folder_count']:,} curated E3 ligases. Across scaffold annotations, "
        f"the dataset contains {overview['unique_scaffold_count']:,} unique scaffolds and "
        f"{overview['unique_supercluster_count']:,} scaffold superclusters, supporting ligand-, conformation-, "
        f"ligase-, and scaffold-centric analysis of E3 recruiter chemical space."
    ))
    report.append(wrap_paragraph(
        f"In this dataset, {overview['ligands_with_multiple_conformations']:,} ligand/component names are associated "
        "with more than one LR recruiter-code tag, indicating multiple curated structures, poses, or conformational "
        "representations for the same underlying ligand. These repeated entries are retained because they preserve "
        "structural context and binding-mode variability rather than being collapsed into a single representative."
    ))

    report.append(section("11. REPRODUCIBILITY NOTES"))
    report.append(wrap_paragraph(
        "This summary was generated directly from the local Ligase_Table CSV files. The authoritative ligand identity "
        "mapping is Recruiter_Master_Map.csv, where the original ligand/component name defines unique ligand identity "
        "and the LR recruiter code defines a curated recruiter/conformation/instance tag. If upstream cleaning, scaffold "
        "unification, SASA calculation, or descriptor-generation scripts are rerun, rerun this summary script afterward "
        "so manuscript numbers remain synchronized with the database."
    ))

    if write_derived:
        report.append(subsection("Derived files written"))
        report.append(str(table_dir / "Derived_Ligand_Conformation_Summary.csv") + "\n")
        report.append(str(table_dir / "Derived_RecruiterCode_To_Ligand_Map.csv") + "\n")

    return "\n".join(report)


# =============================================================================
# Main entry point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate a manuscript-ready summary of the E3 Recruiter Ligandalyzer dataset."
    )

    parser.add_argument(
        "--table-dir",
        default="Ligase_Table",
        help="Directory containing Ligase_Table CSV files. Default: Ligase_Table",
    )

    parser.add_argument(
        "--ligase-dir",
        default="Ligases",
        help="Directory containing per-ligase folders. Default: Ligases",
    )

    parser.add_argument(
        "--out",
        default="E3Recruiter_Data_Summary.txt",
        help="Output text report path. Default: E3Recruiter_Data_Summary.txt",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=25,
        help="Number of top ligands/scaffolds/atoms to include. Default: 25",
    )

    parser.add_argument(
        "--no-derived",
        action="store_true",
        help="Do not write derived CSV files into Ligase_Table.",
    )

    args = parser.parse_args()

    table_dir = Path(args.table_dir)
    ligase_dir = Path(args.ligase_dir)
    out_path = Path(args.out)

    print(f"[INFO] Loading CSV tables from: {table_dir.resolve()}")
    tables = load_tables(table_dir)

    if not tables:
        raise RuntimeError(f"No CSV tables found in {table_dir}")

    print(f"[INFO] Loaded {len(tables):,} CSV tables.")
    print(f"[INFO] Using {PREFERRED_MASTER_MAP} to count unique original ligand/component names.")
    print("[INFO] Building report...")

    report = build_report(
        tables=tables,
        table_dir=table_dir,
        ligase_dir=ligase_dir,
        top_n=args.top_n,
        write_derived=not args.no_derived,
    )

    out_path.write_text(report, encoding="utf-8")
    print(f"[DONE] Wrote summary report to: {out_path.resolve()}")

    master_map, _, _ = build_master_mapping(tables)
    ligand_conf_summary, _ = summarize_ligand_conformations(master_map, tables)
    overview = summarize_dataset_overview(tables, ligase_dir, master_map, ligand_conf_summary)

    print("")
    print("Key counts:")
    print(f"  Curated E3 ligase folders:          {overview['ligase_folder_count']:,}")
    print(f"  Clean CSV ligases detected:         {overview['csv_ligase_count_cleaned']:,}")
    print(f"  Unique ligand/component names:      {overview['unique_ligand_count']:,}")
    print(f"  LR recruiter-code conformation tags:{overview['lr_recruiter_code_count']:,}")
    print(f"  Ligands with >1 LR code:            {overview['ligands_with_multiple_conformations']:,}")
    print(f"  Unique PDB / structure IDs:         {overview['unique_pdb_count']:,}")
    print(f"  Unique scaffolds:                   {overview['unique_scaffold_count']:,}")
    print(f"  Unique scaffold superclusters:      {overview['unique_supercluster_count']:,}")
    print("")
    print("Open the report with:")
    print(f"  less {out_path}")
    print("")
    print("Or on macOS:")
    print(f"  open {out_path}")


if __name__ == "__main__":
    main()
