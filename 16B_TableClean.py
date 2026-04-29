#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
16B_TableClean.py
===============================================================================
Purpose:
    Replace ligase-local Scaffold_ID values with unified LR-SCAF IDs using:

        Ligase_Table/Scaffold_Unified_Map.csv

    This keeps scaffold IDs consistent across downstream database tables.

Primary target tables:
    Ligase_Table/Ligase_Scaffold_Data.csv
    Ligase_Table/Ligase_Scaffold_Frequency.csv

Also updates if present:
    Ligase_Table/Ligase_Recruiters_Scaffold.csv
    Ligase_Table/Ligase_Recruiters_Superclustered.csv

Mapping input:
    Scaffold_Unified_Map.csv

Expected columns:
    Unified_Scaffold_ID
    Scaffold_Hash
    Original_Scaffold_ID

Behavior:
    - Scaffold_ID is replaced with Unified_Scaffold_ID.
    - Original scaffold ID is preserved in Original_Scaffold_ID.
    - Frequency/data tables are regrouped by Ligase + unified Scaffold_ID.
    - Recruiter_Count is summed after unification.
    - Diversity metrics are recomputed where applicable.
    - Original scaffold IDs are kept in Original_Scaffold_ID_Group.

Usage:
    python 16B_TableClean.py --dry-run
    python 16B_TableClean.py --overwrite
===============================================================================
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict

import pandas as pd


# =============================================================================
# Defaults
# =============================================================================

TABLE_DIR = Path("Ligase_Table")

DEFAULT_UNIFIED_MAP = TABLE_DIR / "Scaffold_Unified_Map.csv"

SCAFFOLD_DATA = TABLE_DIR / "Ligase_Scaffold_Data.csv"
SCAFFOLD_FREQ = TABLE_DIR / "Ligase_Scaffold_Frequency.csv"

RECRUITER_SCAFFOLD = TABLE_DIR / "Ligase_Recruiters_Scaffold.csv"
RECRUITER_SUPERCLUSTERED = TABLE_DIR / "Ligase_Recruiters_Superclustered.csv"

AUDIT_OUT = TABLE_DIR / "16B_TableClean_Audit.csv"
LOG_OUT = TABLE_DIR / "16B_TableClean.log"


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


def join_unique(values, sep=";") -> str:
    vals = sorted({clean_str(v) for v in values if clean_str(v)})
    return sep.join(vals)


def first_nonempty(values) -> str:
    for v in values:
        s = clean_str(v)
        if s:
            return s
    return ""


def weighted_mean(group: pd.DataFrame, column: str, weight_col: str = "Recruiter_Count"):
    if column not in group.columns:
        return ""
    vals = pd.to_numeric(group[column], errors="coerce")
    if vals.dropna().empty:
        return ""
    if weight_col in group.columns:
        weights = pd.to_numeric(group[weight_col], errors="coerce").fillna(1)
        try:
            return float((vals.fillna(0) * weights).sum() / weights[vals.notna()].sum())
        except Exception:
            return float(vals.mean())
    return float(vals.mean())


def shannon_index(counts) -> float:
    total = sum(counts)
    if total <= 0:
        return 0.0
    return -sum((c / total) * math.log(c / total) for c in counts if c > 0)


# =============================================================================
# Load scaffold map
# =============================================================================

def load_unified_map(path: Path) -> pd.DataFrame:
    require_file(path)
    df = pd.read_csv(path)

    required = {"Unified_Scaffold_ID", "Original_Scaffold_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")

    if "Scaffold_Hash" not in df.columns:
        df["Scaffold_Hash"] = ""

    df = df.copy()
    df["Unified_Scaffold_ID"] = df["Unified_Scaffold_ID"].astype(str).str.strip()
    df["Original_Scaffold_ID"] = df["Original_Scaffold_ID"].astype(str).str.strip()
    df["Scaffold_Hash"] = df["Scaffold_Hash"].astype(str).str.strip()

    df = df[
        (df["Unified_Scaffold_ID"] != "") &
        (df["Original_Scaffold_ID"] != "")
    ].copy()

    # Verify each original scaffold maps to one unified scaffold.
    bad = (
        df.groupby("Original_Scaffold_ID")["Unified_Scaffold_ID"]
        .nunique()
        .reset_index(name="n")
    )
    bad = bad[bad["n"] > 1]
    if len(bad):
        raise ValueError(
            "Some Original_Scaffold_ID values map to multiple Unified_Scaffold_ID values. "
            "Fix Scaffold_Unified_Map.csv before continuing."
        )

    return df.drop_duplicates(subset=["Original_Scaffold_ID", "Unified_Scaffold_ID"])


def build_maps(unified: pd.DataFrame):
    original_to_unified = dict(zip(
        unified["Original_Scaffold_ID"],
        unified["Unified_Scaffold_ID"],
    ))

    original_to_hash = dict(zip(
        unified["Original_Scaffold_ID"],
        unified["Scaffold_Hash"],
    ))

    hash_to_unified = {}
    if "Scaffold_Hash" in unified.columns:
        tmp = unified[unified["Scaffold_Hash"].astype(str).str.strip() != ""]
        for scaffold_hash, group in tmp.groupby("Scaffold_Hash"):
            ids = sorted(set(group["Unified_Scaffold_ID"]))
            if len(ids) == 1:
                hash_to_unified[scaffold_hash] = ids[0]

    return original_to_unified, original_to_hash, hash_to_unified


# =============================================================================
# Row-level Scaffold_ID replacement
# =============================================================================

def apply_unified_scaffold_ids(
    df: pd.DataFrame,
    original_to_unified: Dict[str, str],
    original_to_hash: Dict[str, str],
    hash_to_unified: Dict[str, str],
    table_name: str,
):
    df = df.copy()

    if "Scaffold_ID" not in df.columns:
        raise ValueError(f"{table_name} does not contain Scaffold_ID column.")

    df["Scaffold_ID"] = df["Scaffold_ID"].astype(str).str.strip()

    if "Original_Scaffold_ID" not in df.columns:
        df["Original_Scaffold_ID"] = df["Scaffold_ID"]

    df["Original_Scaffold_ID"] = df["Original_Scaffold_ID"].astype(str).str.strip()

    audit_rows = []

    new_ids = []
    new_hashes = []

    for _, row in df.iterrows():
        old_id = clean_str(row["Original_Scaffold_ID"])
        current_id = clean_str(row["Scaffold_ID"])
        scaffold_hash = clean_str(row.get("Scaffold_Hash", ""))

        unified_id = ""

        if old_id in original_to_unified:
            unified_id = original_to_unified[old_id]
            scaffold_hash = original_to_hash.get(old_id, scaffold_hash)
        elif current_id in original_to_unified:
            unified_id = original_to_unified[current_id]
            scaffold_hash = original_to_hash.get(current_id, scaffold_hash)
            old_id = current_id
        elif scaffold_hash and scaffold_hash in hash_to_unified:
            unified_id = hash_to_unified[scaffold_hash]
        elif current_id.startswith("LR-SCAF"):
            unified_id = current_id
        else:
            unified_id = current_id
            audit_rows.append({
                "Table": table_name,
                "Issue": "NO_UNIFIED_MAP_MATCH",
                "Original_Scaffold_ID": old_id,
                "Scaffold_ID": current_id,
                "Scaffold_Hash": scaffold_hash,
                "Assigned_Scaffold_ID": unified_id,
                "Notes": "Kept existing Scaffold_ID because no mapping was found.",
            })

        new_ids.append(unified_id)
        new_hashes.append(scaffold_hash)

    df["Scaffold_ID"] = new_ids

    if "Scaffold_Hash" in df.columns:
        df["Scaffold_Hash"] = new_hashes
    else:
        df.insert(
            df.columns.get_loc("Scaffold_ID") + 1,
            "Scaffold_Hash",
            new_hashes,
        )

    audit = pd.DataFrame(
        audit_rows,
        columns=[
            "Table",
            "Issue",
            "Original_Scaffold_ID",
            "Scaffold_ID",
            "Scaffold_Hash",
            "Assigned_Scaffold_ID",
            "Notes",
        ],
    )

    return df, audit


# =============================================================================
# Aggregation for scaffold data/frequency tables
# =============================================================================

def aggregate_scaffold_frequency(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    required = {"Ligase", "Scaffold_ID", "Recruiter_Count"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Ligase_Scaffold_Frequency missing required columns: {missing}")

    df["Recruiter_Count"] = pd.to_numeric(df["Recruiter_Count"], errors="coerce").fillna(0).astype(int)

    grouped = (
        df.groupby(["Ligase", "Scaffold_ID"], dropna=False)
        .agg(
            Recruiter_Count=("Recruiter_Count", "sum"),
            Scaffold_Hash=("Scaffold_Hash", first_nonempty),
            Original_Scaffold_ID_Group=("Original_Scaffold_ID", join_unique),
        )
        .reset_index()
    )

    grouped = grouped.sort_values(
        ["Ligase", "Recruiter_Count", "Scaffold_ID"],
        ascending=[True, False, True],
    ).reset_index(drop=True)

    return grouped


def aggregate_scaffold_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    required = {"Ligase", "Scaffold_ID", "Recruiter_Count"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Ligase_Scaffold_Data missing required columns: {missing}")

    df["Recruiter_Count"] = pd.to_numeric(df["Recruiter_Count"], errors="coerce").fillna(0).astype(int)

    rows = []

    for (ligase, scaffold_id), group in df.groupby(["Ligase", "Scaffold_ID"], dropna=False):
        row = {
            "Ligase": ligase,
            "Scaffold_ID": scaffold_id,
            "Recruiter_Count": int(group["Recruiter_Count"].sum()),
            "Scaffold_Hash": first_nonempty(group.get("Scaffold_Hash", [])),
            "Original_Scaffold_ID_Group": join_unique(group["Original_Scaffold_ID"]),
        }

        # Preserve common descriptor-ish columns.
        for col in [
            "Scaffold_SMILES",
            "Murcko_SMILES",
            "Scaffold_Class",
        ]:
            if col in group.columns:
                row[col] = first_nonempty(group[col])

        # Weighted numeric columns.
        for col in [
            "Scaffold_Center_of_Mass_X",
            "Scaffold_Center_of_Mass_Y",
            "Scaffold_Center_of_Mass_Z",
        ]:
            if col in group.columns:
                row[col] = weighted_mean(group, col, weight_col="Recruiter_Count")

        rows.append(row)

    out = pd.DataFrame(rows)

    # Recompute total recruiters per ligase and density.
    totals = out.groupby("Ligase")["Recruiter_Count"].sum().to_dict()
    out["Total_Recruiters"] = out["Ligase"].map(totals)
    out["Recruiter_Density_Score"] = out["Recruiter_Count"] / out["Total_Recruiters"]

    # Recompute scaffold connectivity: number of ligases per unified scaffold.
    connectivity = out.groupby("Scaffold_ID")["Ligase"].nunique().to_dict()
    out["Ligase_Scaffold_Connectivity"] = out["Scaffold_ID"].map(connectivity)

    # Recompute Shannon and normalized diversity by ligase.
    shannon = {}
    normalized = {}

    for ligase, group in out.groupby("Ligase"):
        counts = group["Recruiter_Count"].astype(float).tolist()
        H = shannon_index(counts)
        Hmax = math.log(len(counts)) if len(counts) > 1 else 0.0
        shannon[ligase] = H
        normalized[ligase] = H / Hmax if Hmax > 0 else 0.0

    out["Shannon_Diversity_Index"] = out["Ligase"].map(shannon)
    out["Normalized_Diversity"] = out["Ligase"].map(normalized)

    preferred_order = [
        "Ligase",
        "Scaffold_ID",
        "Recruiter_Count",
        "Scaffold_Hash",
        "Original_Scaffold_ID_Group",
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

    cols = [c for c in preferred_order if c in out.columns] + [
        c for c in out.columns if c not in preferred_order
    ]

    out = out[cols].sort_values(
        ["Ligase", "Recruiter_Count", "Scaffold_ID"],
        ascending=[True, False, True],
    ).reset_index(drop=True)

    return out


def update_row_level_table(df: pd.DataFrame) -> pd.DataFrame:
    # For recruiter-level tables, no aggregation.
    preferred = [
        "Ligase",
        "Scaffold_ID",
        "Original_Scaffold_ID",
        "RECRUITER_CODE",
        "Scaffold_SMILES",
        "Scaffold_Hash",
    ]
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    return df[cols]


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Replace local Scaffold_IDs with unified LR-SCAF IDs in scaffold tables."
    )

    parser.add_argument("--table-dir", default=str(TABLE_DIR))
    parser.add_argument("--unified-map", default=str(DEFAULT_UNIFIED_MAP))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")

    args = parser.parse_args()

    table_dir = Path(args.table_dir)
    unified_map_path = Path(args.unified_map)

    scaffold_data_path = table_dir / SCAFFOLD_DATA.name
    scaffold_freq_path = table_dir / SCAFFOLD_FREQ.name
    recruiter_scaffold_path = table_dir / RECRUITER_SCAFFOLD.name
    recruiter_superclustered_path = table_dir / RECRUITER_SUPERCLUSTERED.name

    audit_out = table_dir / AUDIT_OUT.name
    log_out = table_dir / LOG_OUT.name

    targets = [
        scaffold_data_path,
        scaffold_freq_path,
    ]

    optional_targets = [
        recruiter_scaffold_path,
        recruiter_superclustered_path,
    ]

    existing_targets = [p for p in targets + optional_targets if p.exists()]

    print("=== 16B Scaffold Table Clean ===")
    print(f"Started: {now()}")
    print(f"Unified map: {unified_map_path}")
    print(f"Table dir:   {table_dir}")
    print(f"Dry run:     {args.dry_run}")

    for p in targets:
        require_file(p)

    unified = load_unified_map(unified_map_path)
    original_to_unified, original_to_hash, hash_to_unified = build_maps(unified)

    print(f"\n📥 Unified map rows: {len(unified)}")
    print(f"🔹 Original scaffold IDs mapped: {len(original_to_unified)}")
    print(f"🔹 Scaffold hashes mapped:       {len(hash_to_unified)}")

    if not args.dry_run:
        existing_outputs = existing_targets + [audit_out, log_out]
        existing_outputs = [p for p in existing_outputs if p.exists()]
        if existing_outputs and not args.overwrite:
            print("\n❌ Output/input files already exist and will be modified:")
            for p in existing_outputs:
                print(f"   - {p}")
            print("Run with --overwrite to back them up and update.")
            return 2

        if args.overwrite:
            for p in existing_outputs:
                b = backup_if_exists(p)
                if b:
                    print(f"🧷 Backed up {p} -> {b}")

    audit_frames = []
    output_counts = {}

    # -------------------------------------------------------------------------
    # Ligase_Scaffold_Frequency.csv
    # -------------------------------------------------------------------------
    freq_raw = pd.read_csv(scaffold_freq_path)
    freq_mapped, audit = apply_unified_scaffold_ids(
        freq_raw,
        original_to_unified,
        original_to_hash,
        hash_to_unified,
        table_name="Ligase_Scaffold_Frequency",
    )
    audit_frames.append(audit)
    freq_clean = aggregate_scaffold_frequency(freq_mapped)
    output_counts["Ligase_Scaffold_Frequency"] = {
        "input_rows": len(freq_raw),
        "output_rows": len(freq_clean),
        "unique_scaffold_ids": int(freq_clean["Scaffold_ID"].nunique()),
    }

    # -------------------------------------------------------------------------
    # Ligase_Scaffold_Data.csv
    # -------------------------------------------------------------------------
    data_raw = pd.read_csv(scaffold_data_path)
    data_mapped, audit = apply_unified_scaffold_ids(
        data_raw,
        original_to_unified,
        original_to_hash,
        hash_to_unified,
        table_name="Ligase_Scaffold_Data",
    )
    audit_frames.append(audit)
    data_clean = aggregate_scaffold_data(data_mapped)
    output_counts["Ligase_Scaffold_Data"] = {
        "input_rows": len(data_raw),
        "output_rows": len(data_clean),
        "unique_scaffold_ids": int(data_clean["Scaffold_ID"].nunique()),
    }

    # -------------------------------------------------------------------------
    # Optional recruiter-level tables
    # -------------------------------------------------------------------------
    optional_outputs = {}

    if recruiter_scaffold_path.exists():
        rec_raw = pd.read_csv(recruiter_scaffold_path)
        rec_clean, audit = apply_unified_scaffold_ids(
            rec_raw,
            original_to_unified,
            original_to_hash,
            hash_to_unified,
            table_name="Ligase_Recruiters_Scaffold",
        )
        audit_frames.append(audit)
        rec_clean = update_row_level_table(rec_clean)
        optional_outputs[recruiter_scaffold_path] = rec_clean
        output_counts["Ligase_Recruiters_Scaffold"] = {
            "input_rows": len(rec_raw),
            "output_rows": len(rec_clean),
            "unique_scaffold_ids": int(rec_clean["Scaffold_ID"].nunique()),
        }

    if recruiter_superclustered_path.exists():
        sup_raw = pd.read_csv(recruiter_superclustered_path)
        sup_clean, audit = apply_unified_scaffold_ids(
            sup_raw,
            original_to_unified,
            original_to_hash,
            hash_to_unified,
            table_name="Ligase_Recruiters_Superclustered",
        )
        audit_frames.append(audit)
        sup_clean = update_row_level_table(sup_clean)
        optional_outputs[recruiter_superclustered_path] = sup_clean
        output_counts["Ligase_Recruiters_Superclustered"] = {
            "input_rows": len(sup_raw),
            "output_rows": len(sup_clean),
            "unique_scaffold_ids": int(sup_clean["Scaffold_ID"].nunique()),
        }

    audit_df = pd.concat(
        [a for a in audit_frames if a is not None and len(a)],
        ignore_index=True,
    ) if any(len(a) for a in audit_frames) else pd.DataFrame(
        columns=[
            "Table",
            "Issue",
            "Original_Scaffold_ID",
            "Scaffold_ID",
            "Scaffold_Hash",
            "Assigned_Scaffold_ID",
            "Notes",
        ]
    )

    print("\n=== Planned Output Counts ===")
    for name, counts in output_counts.items():
        print(
            f"{name}: {counts['input_rows']} rows → {counts['output_rows']} rows | "
            f"{counts['unique_scaffold_ids']} unified scaffold IDs"
        )
    print(f"Audit rows: {len(audit_df)}")

    if args.dry_run:
        print("\n🧪 Dry run complete. Nothing written.")
        print("Run for real with:")
        print("  python 16B_TableClean.py --overwrite")
        return 0

    # Write back in-place.
    freq_clean.to_csv(scaffold_freq_path, index=False)
    data_clean.to_csv(scaffold_data_path, index=False)

    for path, df in optional_outputs.items():
        df.to_csv(path, index=False)

    audit_df.to_csv(audit_out, index=False)

    log_payload = {
        "started_finished": now(),
        "unified_map": str(unified_map_path),
        "counts": output_counts,
        "audit_rows": len(audit_df),
        "outputs_updated": [str(p) for p in existing_targets],
    }

    log_out.write_text(json.dumps(log_payload, indent=2) + "\n")

    print("\n✅ Updated scaffold IDs in-place.")
    print(f"✅ Updated {scaffold_freq_path}")
    print(f"✅ Updated {scaffold_data_path}")
    for path in optional_outputs:
        print(f"✅ Updated {path}")
    print(f"✅ Written audit → {audit_out}")
    print(f"🧾 Written log → {log_out}")

    if len(audit_df):
        print("\n⚠️ Some scaffold IDs had no unified-map match. Review:")
        print(f"  {audit_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
