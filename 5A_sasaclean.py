#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
5A_PDBCLEAN.py
===============================================================================
Purpose:
    Retire unwanted ligand PDB files before SASA / SDF rebuild / metadata steps.

Why:
    Some ligand codes correspond to modified peptide-linking residues or BRD-like
    fragment/modified-residue artifacts that are not useful for the current E3
    recruiter dataset.

    Examples:
        DAR, DPN, DTR, DAL, DAS, DGL, DLE, DPR

    These can appear as ligand-specific PDB files like:
        Ligases/MDM2/PDB/3IWY_DTR.pdb
        Ligases/MDM2/PDB/7KJM_DPN.pdb

What this script does:
    1. Scans Ligases/<Ligase>/PDB/*.pdb.
    2. Parses the ligand code from the filename.
    3. If the ligand code is in the target retirement list, moves the PDB to:
          Retired_PDBs/PDBCLEAN_<timestamp>/Ligases/<Ligase>/PDB/<file>.pdb
    4. Optionally also retires matching old SDF/SDF_4Download files.
    5. Writes:
          5A_PDBCLEAN_manifest.csv
          5A_PDBCLEAN_summary.txt

Usage:
    Preview:
        python 5A_PDBCLEAN.py --dry-run

    Run for real:
        python 5A_PDBCLEAN.py --overwrite

    Custom targets:
        python 5A_PDBCLEAN.py --targets DAR DPN DTR --dry-run

Default retired ligands:
    DAL DAR DAS DGL DLE DPN DPR DTR

Recommended next steps:
    python 7_LigandPDBchecker_LigaseAware.py
    python 5A_PDBCLEAN.py --overwrite
    python 7_LigandPDBchecker_LigaseAware.py
    python 5_SASA_LIGASES_AUDITED.py --overwrite
    python 6_CleanSASA.py
    python 7B_RenameDuplicateLigandCodes.py --all-duplicates --dry-run
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import shutil
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple


DEFAULT_TARGETS = {
    "DAL", "DAR", "DAS", "DGL", "DLE", "DPN", "DPR", "DTR", "OEH", "DTY", "TDF", "O4B", "BIF", "DHI", "WHL", "D0C", "OEH"
}

MANIFEST_COLUMNS = [
    "Action",
    "Reason",
    "Ligase",
    "pdb_id",
    "Ligand",
    "Variant",
    "Original_Path",
    "Retired_Path",
    "Status",
    "Notes",
]


def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def parse_pdb_filename(pdb_path: Path) -> Tuple[str, str, str]:
    """
    Parse:
        3IWY_DTR.pdb    -> pdb_id=3IWY, ligand=DTR, variant=1
        8F0Z_DPN_2.pdb  -> pdb_id=8F0Z, ligand=DPN, variant=2
    """
    parts = pdb_path.stem.split("_")
    pdb_id = parts[0].upper() if len(parts) >= 1 else pdb_path.stem.upper()
    ligand = parts[1].upper() if len(parts) >= 2 else "UNK"
    variant = parts[2] if len(parts) >= 3 and parts[2].isdigit() else "1"
    return pdb_id, ligand, variant


def infer_ligase(path: Path, ligase_root: Path) -> str:
    """
    Expected:
        Ligases/MDM2/PDB/3IWY_DTR.pdb
    """
    try:
        rel = path.relative_to(ligase_root)
        return rel.parts[0]
    except Exception:
        if path.parent.name.upper() in {"PDB", "SDF", "SDF_4DOWNLOAD"}:
            return path.parent.parent.name
        return path.parent.name


def unique_path(dest: Path) -> Path:
    if not dest.exists():
        return dest

    parent = dest.parent
    stem = dest.stem
    suffix = dest.suffix

    i = 1
    while True:
        candidate = parent / f"{stem}.dup{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1


def retire_file(src: Path, ligase_root: Path, retired_root: Path, dry_run: bool) -> Path:
    """
    Preserve path relative to project root if possible.
    Example:
        Ligases/MDM2/PDB/3IWY_DTR.pdb
    becomes:
        Retired_PDBs/PDBCLEAN_stamp/Ligases/MDM2/PDB/3IWY_DTR.pdb
    """
    try:
        rel = src.relative_to(Path("."))
    except Exception:
        rel = src

    dest = unique_path(retired_root / rel)

    if not dry_run:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dest))

    return dest


def discover_pdbs(ligase_root: Path) -> List[Path]:
    return sorted(ligase_root.rglob("PDB/*.pdb"))


def discover_matching_sdfs(ligase_root: Path, ligase: str, pdb_id: str, ligand: str, variant: str) -> List[Path]:
    """
    Finds old SDF files that are clearly associated with a retired PDB/ligand.

    This catches:
        Ligases/<Ligase>/SDF/<Ligase>_<Ligand>.sdf
        Ligases/<Ligase>/SDF_4Download/<PDB>_<Ligand>.sdf
        Ligases/<Ligase>/SDF_4Download/<PDB>_<Ligand>_<variant>.sdf
    """
    ligase_dir = ligase_root / ligase

    candidates = []

    # Representative SDF
    candidates.append(ligase_dir / "SDF" / f"{ligase}_{ligand}.sdf")

    # Per-PDB SDF
    candidates.append(ligase_dir / "SDF_4Download" / f"{pdb_id}_{ligand}.sdf")
    if variant != "1":
        candidates.append(ligase_dir / "SDF_4Download" / f"{pdb_id}_{ligand}_{variant}.sdf")

    return [p for p in candidates if p.exists()]


def write_manifest(path: Path, rows: List[Dict[str, object]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Retire unwanted ligand PDB files before SASA and metadata generation."
    )

    parser.add_argument(
        "--ligase-root",
        default="Ligases",
        help="Root folder containing Ligases/<Ligase>/PDB/*.pdb",
    )

    parser.add_argument(
        "--targets",
        nargs="+",
        default=sorted(DEFAULT_TARGETS),
        help="Ligand codes to retire. Default: DAL DAR DAS DGL DLE DPN DPR DTR",
    )

    parser.add_argument(
        "--retired-root",
        default=None,
        help="Retirement root. Default: Retired_PDBs/PDBCLEAN_<timestamp>",
    )

    parser.add_argument(
        "--include-sdfs",
        action="store_true",
        help="Also retire matching old SDF/SDF_4Download files.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview files that would be moved.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Required for real file moves.",
    )

    args = parser.parse_args()

    ligase_root = Path(args.ligase_root)
    targets = {x.upper().strip() for x in args.targets}

    run_stamp = timestamp()
    retired_root = Path(args.retired_root) if args.retired_root else Path("Retired_PDBs") / f"PDBCLEAN_{run_stamp}"

    if not ligase_root.exists():
        print(f"❌ Ligase root not found: {ligase_root}")
        return 2

    if not args.dry_run and not args.overwrite:
        print("❌ Refusing to move files without --overwrite.")
        print("Preview first:")
        print("  python 5A_PDBCLEAN.py --dry-run")
        print("Then run:")
        print("  python 5A_PDBCLEAN.py --overwrite")
        return 2

    all_pdbs = discover_pdbs(ligase_root)

    matched_pdbs = []
    for pdb in all_pdbs:
        pdb_id, ligand, variant = parse_pdb_filename(pdb)
        if ligand in targets:
            ligase = infer_ligase(pdb, ligase_root)
            matched_pdbs.append((pdb, ligase, pdb_id, ligand, variant))

    print("=== 5A PDB Clean ===")
    print(f"Started: {now()}")
    print(f"Ligase root: {ligase_root}")
    print(f"Total PDB files discovered: {len(all_pdbs)}")
    print(f"Target ligand codes: {', '.join(sorted(targets))}")
    print(f"Matching PDB files to retire: {len(matched_pdbs)}")
    print(f"Retired root: {retired_root}")
    print(f"Also retire matching SDF files: {args.include_sdfs}")
    print(f"Dry run: {args.dry_run}")

    manifest_rows: List[Dict[str, object]] = []
    ligand_counts = Counter()
    ligase_counts = Counter()
    pdb_by_ligand = defaultdict(list)

    for pdb, ligase, pdb_id, ligand, variant in matched_pdbs:
        reason = f"TARGET_LIGAND_{ligand}_MODIFIED_OR_UNWANTED_RESIDUE"
        retired_path = retire_file(
            src=pdb,
            ligase_root=ligase_root,
            retired_root=retired_root,
            dry_run=args.dry_run,
        )

        ligand_counts[ligand] += 1
        ligase_counts[ligase] += 1
        pdb_by_ligand[ligand].append(pdb_id)

        manifest_rows.append({
            "Action": "RETIRE_PDB",
            "Reason": reason,
            "Ligase": ligase,
            "pdb_id": pdb_id,
            "Ligand": ligand,
            "Variant": variant,
            "Original_Path": str(pdb),
            "Retired_Path": str(retired_path),
            "Status": "DRY_RUN" if args.dry_run else "MOVED",
            "Notes": "",
        })

        if args.include_sdfs:
            sdfs = discover_matching_sdfs(ligase_root, ligase, pdb_id, ligand, variant)
            for sdf in sdfs:
                retired_sdf = retire_file(
                    src=sdf,
                    ligase_root=ligase_root,
                    retired_root=retired_root,
                    dry_run=args.dry_run,
                )
                manifest_rows.append({
                    "Action": "RETIRE_SDF",
                    "Reason": f"ASSOCIATED_WITH_RETIRED_PDB_{pdb_id}_{ligand}",
                    "Ligase": ligase,
                    "pdb_id": pdb_id,
                    "Ligand": ligand,
                    "Variant": variant,
                    "Original_Path": str(sdf),
                    "Retired_Path": str(retired_sdf),
                    "Status": "DRY_RUN" if args.dry_run else "MOVED",
                    "Notes": "",
                })

    manifest_path = Path("5A_PDBCLEAN_manifest.csv")
    summary_path = Path("5A_PDBCLEAN_summary.txt")

    write_manifest(manifest_path, manifest_rows)

    summary_lines = []
    summary_lines.append("=== 5A PDBCLEAN Summary ===")
    summary_lines.append(f"Started/finished: {now()}")
    summary_lines.append(f"Ligase root: {ligase_root}")
    summary_lines.append(f"Total PDB files discovered: {len(all_pdbs)}")
    summary_lines.append(f"Target ligand codes: {', '.join(sorted(targets))}")
    summary_lines.append(f"Matching PDB files retired: {len(matched_pdbs)}")
    summary_lines.append(f"Manifest rows: {len(manifest_rows)}")
    summary_lines.append(f"Retired root: {retired_root}")
    summary_lines.append(f"Dry run: {args.dry_run}")
    summary_lines.append("")
    summary_lines.append("Counts by ligand:")
    for lig, count in sorted(ligand_counts.items()):
        pdbs = ", ".join(sorted(set(pdb_by_ligand[lig])))
        summary_lines.append(f"  {lig}: {count} PDB files -> {pdbs}")
    summary_lines.append("")
    summary_lines.append("Counts by ligase:")
    for ligase, count in sorted(ligase_counts.items()):
        summary_lines.append(f"  {ligase}: {count}")

    summary_path.write_text("\n".join(summary_lines) + "\n")

    print("\n=== Retirement Preview/Summary ===")
    for lig, count in sorted(ligand_counts.items()):
        pdbs = ", ".join(sorted(set(pdb_by_ligand[lig])))
        print(f"  {lig}: {count} PDB files -> {pdbs}")

    print(f"\n✅ Written manifest: {manifest_path}")
    print(f"✅ Written summary:  {summary_path}")

    if args.dry_run:
        print("\n🧪 Dry run complete. Nothing was moved.")
        print("Run for real with:")
        print("  python 5A_PDBCLEAN.py --overwrite")
        print("")
        print("If you also want to move old matching SDFs:")
        print("  python 5A_PDBCLEAN.py --overwrite --include-sdfs")
    else:
        print(f"\n🧾 Retired files moved under: {retired_root}")
        print("\nNext recommended commands:")
        print("  python 7_LigandPDBchecker_LigaseAware.py")
        print("  python 5_SASA_LIGASES_AUDITED.py --overwrite")
        print("  python 6_CleanSASA.py")
        print("  python 7B_RenameDuplicateLigandCodes.py --all-duplicates --dry-run")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


