#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
7B_RenameDuplicateLigandCodes.py
===============================================================================
Purpose:
    Rename duplicate ligand codes in active PDB files so each ligand instance has
    a unique 3-character PDB-safe residue code.

Why:
    Old ligand codes like A1A, A36, A67, etc. may have been reused during the
    old 5-letter -> 3-letter renaming process. That means the same code can
    represent different molecules across different PDB files.

    This script fixes that upstream by:
      1. Renaming PDB filenames:
            9CUO_A1A.pdb -> 9CUO_E00.pdb
      2. Renaming matching HETATM residue names inside the PDB:
            residue name A1A -> E00
      3. Writing a mapping CSV:
            Ligand_Code_Rename_Map.csv
      4. Backing up original PDB files before changing anything.

Typical use:
    # Preview only A1A renaming:
    python 7B_RenameDuplicateLigandCodes.py --targets A1A --dry-run

    # Actually rename only A1A instances:
    python 7B_RenameDuplicateLigandCodes.py --targets A1A --overwrite

    # Preview all duplicated ligand codes:
    python 7B_RenameDuplicateLigandCodes.py --all-duplicates --dry-run

    # Rename all duplicated ligand codes:
    python 7B_RenameDuplicateLigandCodes.py --all-duplicates --overwrite

After running:
    python 5_SASA_LIGASES_AUDITED.py --overwrite --include-zero-atoms
    python 6_CleanSASA.py
    python 7_LigandPDBchecker_LigaseAware.py
    python 8B_MetaDataClean.py --overwrite
    python 8_Ligand_Metadata_From_SDF_HEAVYSAFE.py --overwrite
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


MAP_COLUMNS = [
    "Ligase",
    "pdb_id",
    "Variant",
    "Old_Ligand",
    "New_Ligand",
    "Original_PDB",
    "New_PDB",
    "Backup_PDB",
    "Atoms_Renamed",
    "Status",
    "Notes",
]


def now_stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def parse_pdb_filename(pdb_path: Path) -> Tuple[str, str, str]:
    """
    Parse filenames like:
        9CUO_A1A.pdb      -> pdb_id=9CUO, ligand=A1A, variant=1
        8GCG_A1A_3.pdb   -> pdb_id=8GCG, ligand=A1A, variant=3
    """
    parts = pdb_path.stem.split("_")
    pdb_id = parts[0].upper() if len(parts) >= 1 else pdb_path.stem.upper()
    ligand = parts[1].upper() if len(parts) >= 2 else "UNK"
    variant = parts[2] if len(parts) >= 3 and parts[2].isdigit() else "1"
    return pdb_id, ligand, variant


def infer_ligase(pdb_path: Path, ligase_root: Path) -> str:
    try:
        rel = pdb_path.relative_to(ligase_root)
        return rel.parts[0]
    except Exception:
        if pdb_path.parent.name.upper() == "PDB":
            return pdb_path.parent.parent.name
        return pdb_path.parent.name


def residue_name_from_pdb_line(line: str) -> str:
    return line[17:20].strip().upper()


def replace_residue_name(line: str, new_resname: str) -> str:
    """
    PDB residue name occupies columns 18-20, Python slice [17:20].
    New residue name must be exactly 3 characters.
    """
    return line[:17] + f"{new_resname:>3}" + line[20:]


def discover_pdbs(ligase_root: Path) -> List[Path]:
    return sorted(ligase_root.rglob("PDB/*.pdb"))


def build_records(ligase_root: Path) -> List[Dict[str, object]]:
    records = []
    for pdb_path in discover_pdbs(ligase_root):
        pdb_id, ligand, variant = parse_pdb_filename(pdb_path)
        ligase = infer_ligase(pdb_path, ligase_root)

        records.append({
            "Ligase": ligase,
            "pdb_id": pdb_id,
            "Ligand": ligand,
            "Variant": variant,
            "Path": pdb_path,
        })

    return records


def generate_code_pool(prefixes: List[str], start: int = 0, stop: int = 99) -> List[str]:
    codes = []
    for prefix in prefixes:
        prefix = prefix.upper().strip()
        if len(prefix) != 1:
            raise ValueError(f"Prefix must be one character for PDB 3-letter codes: {prefix}")
        for i in range(start, stop + 1):
            codes.append(f"{prefix}{i:02d}")
    return codes


def count_atoms_to_rename(pdb_path: Path, old_ligand: str) -> int:
    count = 0
    with pdb_path.open(errors="replace") as f:
        for line in f:
            if line.startswith("HETATM") and residue_name_from_pdb_line(line) == old_ligand.upper():
                count += 1
    return count


def make_new_filename(pdb_path: Path, old_ligand: str, new_ligand: str) -> Path:
    """
    Replace the ligand token in filename, preserving pdb_id and variant.
    """
    pdb_id, ligand, variant = parse_pdb_filename(pdb_path)

    if ligand != old_ligand:
        raise ValueError(f"Filename ligand {ligand} does not match expected old ligand {old_ligand}")

    if variant and variant != "1":
        new_name = f"{pdb_id}_{new_ligand}_{variant}.pdb"
    else:
        new_name = f"{pdb_id}_{new_ligand}.pdb"

    return pdb_path.with_name(new_name)


def rewrite_pdb_with_new_ligand(
    old_path: Path,
    new_path: Path,
    old_ligand: str,
    new_ligand: str,
) -> int:
    """
    Write a new PDB file with HETATM residue names changed from old_ligand to new_ligand.
    Return number of atoms renamed.
    """
    atoms_renamed = 0
    new_lines = []

    with old_path.open(errors="replace") as f:
        for line in f:
            if line.startswith("HETATM") and residue_name_from_pdb_line(line) == old_ligand.upper():
                new_lines.append(replace_residue_name(line, new_ligand))
                atoms_renamed += 1
            else:
                new_lines.append(line)

    new_path.write_text("".join(new_lines))
    return atoms_renamed


def write_mapping_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MAP_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rename duplicate ligand codes in PDB filenames and HETATM residue names."
    )

    parser.add_argument("--ligase-root", default="Ligases", help="Root folder containing <Ligase>/PDB/*.pdb")
    parser.add_argument("--targets", nargs="*", default=None, help="Specific ligand codes to rename, e.g. A1A A36 A67")
    parser.add_argument("--all-duplicates", action="store_true", help="Rename every ligand code that appears in >1 PDB file")
    parser.add_argument("--prefixes", nargs="+", default=["E"], help="One-letter prefixes for new codes. Default: E")
    parser.add_argument("--start", type=int, default=0, help="Starting number for code generation. Default: 0")
    parser.add_argument("--stop", type=int, default=99, help="Ending number for code generation. Default: 99")
    parser.add_argument("--mapping-out", default="Ligand_Code_Rename_Map.csv")
    parser.add_argument("--backup-root", default="Renamed_PDB_Backups")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite", action="store_true", help="Required for real file changes")

    args = parser.parse_args()

    ligase_root = Path(args.ligase_root)
    mapping_out = Path(args.mapping_out)
    backup_root = Path(args.backup_root) / f"Rename_{now_stamp()}"

    if not ligase_root.exists():
        print(f"❌ Ligase root not found: {ligase_root}")
        return 2

    if not args.dry_run and not args.overwrite:
        print("❌ Refusing to modify files without --overwrite.")
        print("Preview first:")
        print("  python 7B_RenameDuplicateLigandCodes.py --targets A1A --dry-run")
        print("Then run:")
        print("  python 7B_RenameDuplicateLigandCodes.py --targets A1A --overwrite")
        return 2

    if not args.targets and not args.all_duplicates:
        print("❌ Choose what to rename:")
        print("  --targets A1A")
        print("or")
        print("  --all-duplicates")
        return 2

    records = build_records(ligase_root)

    by_ligand: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for r in records:
        by_ligand[str(r["Ligand"]).upper()].append(r)

    if args.targets:
        target_codes = {x.upper().strip() for x in args.targets}
    else:
        target_codes = {lig for lig, rows in by_ligand.items() if len(rows) > 1}

    target_records = []
    for lig in sorted(target_codes):
        target_records.extend(by_ligand.get(lig, []))

    target_records = sorted(
        target_records,
        key=lambda r: (str(r["Ligand"]), str(r["Ligase"]), str(r["pdb_id"]), str(r["Variant"]), str(r["Path"]))
    )

    if not target_records:
        print("❌ No matching PDB files found for requested target ligand codes.")
        return 2

    existing_ligand_codes = {str(r["Ligand"]).upper() for r in records}
    old_target_codes = {str(r["Ligand"]).upper() for r in target_records}

    # Avoid assigning a new code that already exists outside the target set.
    reserved_codes = existing_ligand_codes - old_target_codes

    code_pool = generate_code_pool(args.prefixes, args.start, args.stop)
    available_codes = [c for c in code_pool if c not in reserved_codes]

    if len(available_codes) < len(target_records):
        print("❌ Not enough available new ligand codes.")
        print(f"Need: {len(target_records)}")
        print(f"Available: {len(available_codes)}")
        print("Use more prefixes, e.g.: --prefixes E F G H")
        return 2

    print("=== Duplicate Ligand Code Renamer ===")
    print(f"Ligase root: {ligase_root}")
    print(f"Total PDB files discovered: {len(records)}")
    print(f"Target ligand codes: {', '.join(sorted(target_codes))}")
    print(f"Target PDB files to rename: {len(target_records)}")
    print(f"Code range: {available_codes[0]} ... {available_codes[len(target_records)-1]}")
    print(f"Dry run: {args.dry_run}")

    if args.dry_run:
        print("\n🧪 DRY RUN: no files will be changed.")

    rows = []

    for idx, rec in enumerate(target_records):
        old_path: Path = rec["Path"]  # type: ignore[assignment]
        old_ligand = str(rec["Ligand"]).upper()
        new_ligand = available_codes[idx]

        new_path = make_new_filename(old_path, old_ligand, new_ligand)

        atoms_to_rename = count_atoms_to_rename(old_path, old_ligand)

        rel = old_path.relative_to(ligase_root)
        backup_path = backup_root / rel

        status = "DRY_RUN" if args.dry_run else "RENAMED"
        notes = ""

        if new_path.exists() and new_path != old_path:
            status = "SKIPPED"
            notes = f"Destination already exists: {new_path}"
        elif atoms_to_rename == 0:
            status = "SKIPPED"
            notes = f"No HETATM atoms found with residue name {old_ligand}"

        if not args.dry_run and status == "RENAMED":
            backup_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(old_path, backup_path)

            atoms_renamed = rewrite_pdb_with_new_ligand(
                old_path=old_path,
                new_path=new_path,
                old_ligand=old_ligand,
                new_ligand=new_ligand,
            )

            if new_path != old_path:
                old_path.unlink()

            atoms_to_rename = atoms_renamed

        rows.append({
            "Ligase": rec["Ligase"],
            "pdb_id": rec["pdb_id"],
            "Variant": rec["Variant"],
            "Old_Ligand": old_ligand,
            "New_Ligand": new_ligand,
            "Original_PDB": str(old_path),
            "New_PDB": str(new_path),
            "Backup_PDB": str(backup_path) if not args.dry_run else "",
            "Atoms_Renamed": atoms_to_rename,
            "Status": status,
            "Notes": notes,
        })

    write_mapping_csv(mapping_out, rows)

    print(f"\n✅ Written mapping: {mapping_out}")

    print("\nPreview:")
    for row in rows[:25]:
        print(
            f"  {row['Ligase']} {row['pdb_id']} "
            f"{row['Old_Ligand']} -> {row['New_Ligand']} "
            f"atoms={row['Atoms_Renamed']} status={row['Status']}"
        )

    if len(rows) > 25:
        print(f"  ... {len(rows) - 25} more rows")

    if args.dry_run:
        print("\nNext, run for real, for example:")
        if args.targets:
            print(f"  python 7B_RenameDuplicateLigandCodes.py --targets {' '.join(sorted(target_codes))} --overwrite")
        else:
            print("  python 7B_RenameDuplicateLigandCodes.py --all-duplicates --overwrite")
    else:
        print(f"\n🧾 Backups written under: {backup_root}")
        print("\nNow rerun the downstream chain:")
        print("  python 5_SASA_LIGASES_AUDITED.py --overwrite --include-zero-atoms")
        print("  python 6_CleanSASA.py")
        print("  python 7_LigandPDBchecker_LigaseAware.py")
        print("  python 8B_MetaDataClean.py --overwrite")
        print("  python 8_Ligand_Metadata_From_SDF_HEAVYSAFE.py --overwrite")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
