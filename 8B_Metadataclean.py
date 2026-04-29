#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
8B_MetaDataClean.py
===============================================================================
Purpose:
    Rebuild active ligand SDF files directly from the current PDB files.

Why:
    Old Ligases/<Ligase>/SDF/*.sdf files may be stale because old 5-letter ligand
    names were renamed into 3-letter A## style names. This can cause cases where:

        Ligases/TRIM21/SDF/TRIM21_A03.sdf

    is chemically different from the active PDB ligand residue A03 in:

        Ligases/TRIM21/PDB/<PDB>_A03.pdb

    This script retires the old SDF files and regenerates new SDF files from the
    actual PDB ligand residue blocks that were used for SASA.

What it does:
    1. Scans Ligases/<Ligase>/PDB/*.pdb.
    2. Parses pdb_id, ligand, variant from filenames like:
          7HLV_A28_2.pdb
          9QEU_A90.pdb
    3. Extracts exact HETATM residue matching the ligand code.
    4. Uses RDKit to build a molecule from that PDB ligand block.
    5. Writes fresh PDB-derived SDFs to:
          Ligases/<Ligase>/SDF_4Download/<pdb_id>_<ligand>[_variant].sdf
    6. Writes representative current SDFs to:
          Ligases/<Ligase>/SDF/<Ligase>_<Ligand>.sdf
    7. Moves old SDF/SDF_4Download files to:
          Retired_SDFs/Rebuild_<timestamp>/...
    8. Writes audit files:
          8B_MetadataClean_Rebuild_Manifest.csv
          8B_MetadataClean_Group_Collisions.csv
          8B_MetadataClean_Failures.csv
          8B_MetadataClean.log

Important:
    A single Ligase+Ligand code can still represent multiple PDB-instance
    molecules. Example: CRBN_A1A may occur across several PDBs with different
    atom counts. This script will still create per-PDB SDFs in SDF_4Download,
    and one representative in SDF. The collision report tells you where this
    happened.

Usage:
    python 8B_MetaDataClean.py --dry-run
    python 8B_MetaDataClean.py --overwrite

Then rerun:
    python 8_Ligand_Metadata_From_SDF_HEAVYSAFE.py --overwrite
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import traceback
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


# =============================================================================
# Output files
# =============================================================================

MANIFEST_CSV = Path("8B_MetadataClean_Rebuild_Manifest.csv")
COLLISION_CSV = Path("8B_MetadataClean_Group_Collisions.csv")
FAILURE_CSV = Path("8B_MetadataClean_Failures.csv")
LOG_FILE = Path("8B_MetadataClean.log")

MANIFEST_COLUMNS = [
    "Ligase",
    "pdb_id",
    "Ligand",
    "Variant",
    "RECRUITER_CODE",
    "Input_PDB",
    "Output_Instance_SDF",
    "Output_Representative_SDF",
    "Total_Atoms",
    "Heavy_Atoms",
    "Hydrogen_Atoms",
    "Formula",
    "MW",
    "InChIKey",
    "SMILES",
    "Status",
    "Warnings",
]

COLLISION_COLUMNS = [
    "Ligase",
    "Ligand",
    "RECRUITER_CODE",
    "Instance_Count",
    "Unique_Heavy_Atom_Counts",
    "Unique_InChIKeys",
    "Representative_PDB",
    "Representative_SDF",
    "All_Instance_SDFs",
    "Collision_Type",
]

FAILURE_COLUMNS = [
    "Ligase",
    "pdb_id",
    "Ligand",
    "Variant",
    "Input_PDB",
    "Error",
    "Traceback",
]

WATER_RESNAMES = {"HOH", "WAT", "DOD"}


# =============================================================================
# Helpers
# =============================================================================

def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def clean_str(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and math.isnan(x):
        return ""
    return str(x).strip()


def parse_pdb_filename(pdb_path: Path) -> Tuple[str, str, str]:
    """
    Parse:
        7HLV_A28_2.pdb -> pdb_id=7HLV, ligand=A28, variant=2
        9QEU_A90.pdb   -> pdb_id=9QEU, ligand=A90, variant=1
    """
    parts = pdb_path.stem.split("_")
    pdb_id = parts[0].upper() if len(parts) >= 1 else pdb_path.stem.upper()
    ligand = parts[1].upper() if len(parts) >= 2 else "UNK"
    variant = parts[2] if len(parts) >= 3 and parts[2].isdigit() else "1"
    return pdb_id, ligand, variant


def infer_ligase(pdb_path: Path, ligase_root: Path) -> str:
    """
    Expected path:
        Ligases/TRIM21/PDB/9QEU_A90.pdb
    """
    try:
        rel = pdb_path.relative_to(ligase_root)
        return rel.parts[0]
    except Exception:
        if pdb_path.parent.name.upper() == "PDB":
            return pdb_path.parent.parent.name
        return pdb_path.parent.name


def residue_name_from_pdb_line(line: str) -> str:
    return line[17:20].strip().upper()


def atom_serial_from_pdb_line(line: str) -> Optional[int]:
    try:
        return int(line[6:11])
    except Exception:
        return None


def ensure_empty_or_create_dir(path: Path, dry_run: bool = False) -> None:
    if dry_run:
        return
    path.mkdir(parents=True, exist_ok=True)


def unique_path(path: Path) -> Path:
    """
    Avoid overwriting if a destination already exists.
    """
    if not path.exists():
        return path

    parent = path.parent
    stem = path.stem
    suffix = path.suffix

    i = 1
    while True:
        candidate = parent / f"{stem}.dup{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1


def write_csv(path: Path, columns: List[str], rows: List[Dict[str, object]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# =============================================================================
# Retirement of old SDFs
# =============================================================================

def collect_existing_sdfs(ligase_root: Path) -> List[Path]:
    """
    Collect current old SDF files from active SDF and SDF_4Download folders.
    """
    paths = []
    for folder_name in ["SDF", "SDF_4Download"]:
        paths.extend(ligase_root.rglob(f"{folder_name}/*.sdf"))
    return sorted(set(paths))


def retire_existing_sdfs(
    ligase_root: Path,
    retired_root: Path,
    run_stamp: str,
    dry_run: bool,
) -> List[Dict[str, object]]:
    """
    Move old SDF and SDF_4Download files to Retired_SDFs/Rebuild_<timestamp>/...
    preserving relative paths.
    """
    old_sdfs = collect_existing_sdfs(ligase_root)
    retired_rows = []

    for src in old_sdfs:
        try:
            rel = src.relative_to(ligase_root)
        except Exception:
            rel = Path(src.name)

        dest = retired_root / f"Rebuild_{run_stamp}" / rel
        dest = unique_path(dest)

        retired_rows.append({
            "Original": str(src),
            "Retired": str(dest),
        })

        if not dry_run:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dest))

    return retired_rows


def clean_active_sdf_dirs(ligase_root: Path, dry_run: bool) -> None:
    """
    Recreate empty SDF and SDF_4Download folders for each ligase that has PDBs.
    """
    ligase_dirs = sorted({p.parent.parent for p in ligase_root.rglob("PDB/*.pdb")})

    for ligase_dir in ligase_dirs:
        for folder in ["SDF", "SDF_4Download"]:
            outdir = ligase_dir / folder
            if not dry_run:
                outdir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# PDB ligand extraction
# =============================================================================

def read_ligand_pdb_block(pdb_path: Path, ligand: str) -> Tuple[str, List[str], List[int]]:
    """
    Extract exact HETATM lines for the ligand residue name, plus relevant CONECT lines.
    """
    hetatm_lines: List[str] = []
    conect_lines: List[str] = []
    ligand_serials: List[int] = []

    with pdb_path.open(errors="replace") as f:
        all_lines = f.readlines()

    for line in all_lines:
        if not line.startswith("HETATM"):
            continue
        resname = residue_name_from_pdb_line(line)
        if resname == ligand.upper():
            hetatm_lines.append(line)
            serial = atom_serial_from_pdb_line(line)
            if serial is not None:
                ligand_serials.append(serial)

    ligand_serial_set = set(ligand_serials)

    for line in all_lines:
        if not line.startswith("CONECT"):
            continue

        nums = []
        for i in range(6, len(line), 5):
            chunk = line[i:i + 5].strip()
            if chunk.isdigit():
                nums.append(int(chunk))

        if nums and nums[0] in ligand_serial_set:
            # Keep only CONECT records originating from ligand atoms.
            # This avoids pulling in protein connectivity.
            conect_lines.append(line)

    pdb_block = "".join(hetatm_lines + conect_lines)
    return pdb_block, hetatm_lines, ligand_serials


def mol_from_pdb_block_tolerant(pdb_block: str) -> Tuple[Optional[Chem.Mol], str]:
    """
    Build RDKit mol from a PDB ligand block.
    """
    if not pdb_block.strip():
        return None, "EMPTY_PDB_BLOCK"

    errors = []

    # Try with sanitization first.
    try:
        mol = Chem.MolFromPDBBlock(
            pdb_block,
            sanitize=True,
            removeHs=False,
            proximityBonding=True,
        )
        if mol is not None:
            return mol, ""
    except Exception as e:
        errors.append(f"sanitize=True failed: {type(e).__name__}: {e}")

    # Try without sanitization.
    try:
        mol = Chem.MolFromPDBBlock(
            pdb_block,
            sanitize=False,
            removeHs=False,
            proximityBonding=True,
        )
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except Exception as e:
                errors.append(f"post-sanitize failed: {type(e).__name__}: {e}")
            return mol, "; ".join(errors)
    except Exception as e:
        errors.append(f"sanitize=False failed: {type(e).__name__}: {e}")

    return None, "; ".join(errors) if errors else "RDKit returned None"


def compute_mol_info(mol: Chem.Mol) -> Dict[str, object]:
    total_atoms = mol.GetNumAtoms()
    heavy_atoms = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    hydrogen_atoms = total_atoms - heavy_atoms

    try:
        formula = rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        formula = ""

    try:
        mw = round(float(Descriptors.MolWt(mol)), 3)
    except Exception:
        mw = ""

    try:
        mol_no_h = Chem.RemoveHs(mol, sanitize=False)
    except Exception:
        mol_no_h = mol

    try:
        smiles = Chem.MolToSmiles(mol_no_h, canonical=True, isomericSmiles=True)
    except Exception:
        smiles = ""

    try:
        inchi = Chem.MolToInchi(mol_no_h)
        inchikey = Chem.InchiToInchiKey(inchi) if inchi else ""
    except Exception:
        inchikey = ""

    return {
        "Total_Atoms": total_atoms,
        "Heavy_Atoms": heavy_atoms,
        "Hydrogen_Atoms": hydrogen_atoms,
        "Formula": formula,
        "MW": mw,
        "SMILES": smiles,
        "InChIKey": inchikey,
    }


def write_sdf(mol: Chem.Mol, out_path: Path, name: str, dry_run: bool) -> None:
    if dry_run:
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)

    mol_to_write = Chem.Mol(mol)
    mol_to_write.SetProp("_Name", name)

    writer = Chem.SDWriter(str(out_path))
    writer.write(mol_to_write)
    writer.close()


# =============================================================================
# Main rebuild logic
# =============================================================================

def discover_active_pdbs(ligase_root: Path) -> List[Path]:
    return sorted(ligase_root.rglob("PDB/*.pdb"))


def instance_sdf_name(pdb_id: str, ligand: str, variant: str) -> str:
    if str(variant) and str(variant) != "1":
        return f"{pdb_id}_{ligand}_{variant}.sdf"
    return f"{pdb_id}_{ligand}.sdf"


def representative_sdf_name(ligase: str, ligand: str) -> str:
    return f"{ligase}_{ligand}.sdf"


def choose_representative(instances: List[Dict[str, object]]) -> Dict[str, object]:
    """
    Choose representative for Ligases/<Ligase>/SDF/<Ligase>_<Ligand>.sdf.

    Priority:
        1. Most common Heavy_Atoms.
        2. First sorted PDB ID.
    """
    heavy_counts = [row.get("Heavy_Atoms") for row in instances if row.get("Heavy_Atoms") != ""]
    if heavy_counts:
        most_common_heavy = Counter(heavy_counts).most_common(1)[0][0]
        candidates = [r for r in instances if r.get("Heavy_Atoms") == most_common_heavy]
    else:
        candidates = instances

    candidates = sorted(candidates, key=lambda r: (str(r.get("pdb_id", "")), str(r.get("Variant", ""))))
    return candidates[0]


def rebuild_sdfs(
    ligase_root: Path,
    dry_run: bool,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, object]]]:
    pdbs = discover_active_pdbs(ligase_root)

    manifest_rows: List[Dict[str, object]] = []
    failure_rows: List[Dict[str, object]] = []

    # Store actual molecules temporarily so we can write representatives after grouping.
    mol_cache: Dict[str, Chem.Mol] = {}

    for pdb_path in pdbs:
        ligase = infer_ligase(pdb_path, ligase_root)
        pdb_id, ligand, variant = parse_pdb_filename(pdb_path)
        recruiter_code = f"{ligase}_{ligand}"

        try:
            pdb_block, hetatm_lines, ligand_serials = read_ligand_pdb_block(pdb_path, ligand)

            if not hetatm_lines:
                raise ValueError(f"No exact HETATM residue found for ligand {ligand}")

            mol, warning = mol_from_pdb_block_tolerant(pdb_block)

            if mol is None:
                raise ValueError(f"RDKit could not build ligand mol: {warning}")

            info = compute_mol_info(mol)

            instance_name = instance_sdf_name(pdb_id, ligand, variant)
            instance_sdf = ligase_root / ligase / "SDF_4Download" / instance_name

            sdf_name = f"{pdb_id}_{ligand}" + (f"_{variant}" if variant != "1" else "")
            write_sdf(mol, instance_sdf, sdf_name, dry_run=dry_run)

            cache_key = str(instance_sdf)
            mol_cache[cache_key] = mol

            row = {
                "Ligase": ligase,
                "pdb_id": pdb_id,
                "Ligand": ligand,
                "Variant": variant,
                "RECRUITER_CODE": recruiter_code,
                "Input_PDB": str(pdb_path),
                "Output_Instance_SDF": str(instance_sdf),
                "Output_Representative_SDF": "",
                "Total_Atoms": info["Total_Atoms"],
                "Heavy_Atoms": info["Heavy_Atoms"],
                "Hydrogen_Atoms": info["Hydrogen_Atoms"],
                "Formula": info["Formula"],
                "MW": info["MW"],
                "InChIKey": info["InChIKey"],
                "SMILES": info["SMILES"],
                "Status": "OK" if not warning else "CHECK",
                "Warnings": warning,
            }
            manifest_rows.append(row)

        except Exception as e:
            failure_rows.append({
                "Ligase": ligase,
                "pdb_id": pdb_id,
                "Ligand": ligand,
                "Variant": variant,
                "Input_PDB": str(pdb_path),
                "Error": f"{type(e).__name__}: {e}",
                "Traceback": traceback.format_exc(),
            })

    # Group by Ligase+Ligand and write representative SDFs.
    grouped: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for row in manifest_rows:
        grouped[(str(row["Ligase"]), str(row["Ligand"]))].append(row)

    collision_rows: List[Dict[str, object]] = []

    for (ligase, ligand), rows in grouped.items():
        recruiter_code = f"{ligase}_{ligand}"
        representative = choose_representative(rows)

        rep_path = ligase_root / ligase / "SDF" / representative_sdf_name(ligase, ligand)
        instance_path = str(representative["Output_Instance_SDF"])

        mol = mol_cache.get(instance_path)

        if mol is not None:
            write_sdf(mol, rep_path, recruiter_code, dry_run=dry_run)

        # Update manifest rows with representative path.
        for row in rows:
            row["Output_Representative_SDF"] = str(rep_path)

        heavy_set = sorted({str(r.get("Heavy_Atoms", "")) for r in rows if clean_str(r.get("Heavy_Atoms"))})
        inchikey_set = sorted({str(r.get("InChIKey", "")) for r in rows if clean_str(r.get("InChIKey"))})
        instance_sdfs = sorted({str(r.get("Output_Instance_SDF", "")) for r in rows})

        collision_type = "NONE"
        if len(rows) > 1 and (len(heavy_set) > 1 or len(inchikey_set) > 1):
            collision_type = "MULTIPLE_DISTINCT_PDB_INSTANCE_MOLECULES"
        elif len(rows) > 1:
            collision_type = "MULTIPLE_PDB_INSTANCES_SAME_APPARENT_MOLECULE"

        if collision_type != "NONE":
            collision_rows.append({
                "Ligase": ligase,
                "Ligand": ligand,
                "RECRUITER_CODE": recruiter_code,
                "Instance_Count": len(rows),
                "Unique_Heavy_Atom_Counts": ";".join(heavy_set),
                "Unique_InChIKeys": ";".join(inchikey_set),
                "Representative_PDB": str(representative.get("pdb_id", "")),
                "Representative_SDF": str(rep_path),
                "All_Instance_SDFs": ";".join(instance_sdfs),
                "Collision_Type": collision_type,
            })

    return manifest_rows, collision_rows, failure_rows


# =============================================================================
# CLI
# =============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Retire stale SDFs and rebuild fresh active ligand SDFs from current PDB files."
    )

    parser.add_argument("--ligase-root", default="Ligases", help="Root folder containing Ligases/<Ligase>/PDB/*.pdb")
    parser.add_argument("--retired-root", default="Retired_SDFs", help="Where old SDFs are moved")
    parser.add_argument("--overwrite", action="store_true", help="Required for real run unless --dry-run is used")
    parser.add_argument("--dry-run", action="store_true", help="Preview actions without moving/writing SDF files")

    args = parser.parse_args()

    ligase_root = Path(args.ligase_root)
    retired_root = Path(args.retired_root)
    run_stamp = timestamp()

    if not ligase_root.exists():
        print(f"❌ Ligase root not found: {ligase_root}")
        return 2

    if not args.dry_run and not args.overwrite:
        print("❌ Refusing to modify files without --overwrite.")
        print("Run first:")
        print("  python 8B_MetaDataClean.py --dry-run")
        print("Then:")
        print("  python 8B_MetaDataClean.py --overwrite")
        return 2

    pdbs = discover_active_pdbs(ligase_root)
    old_sdfs = collect_existing_sdfs(ligase_root)

    print("=== 8B Metadata/SDF Clean Rebuild ===")
    print(f"Started: {now()}")
    print(f"Ligase root: {ligase_root}")
    print(f"Active PDB files discovered: {len(pdbs)}")
    print(f"Existing SDF/SDF_4Download files to retire: {len(old_sdfs)}")
    print(f"Retired root: {retired_root / ('Rebuild_' + run_stamp)}")
    print(f"Dry run: {args.dry_run}")

    if args.dry_run:
        print("\n🧪 DRY RUN ONLY. No files will be moved or written.")

    # Retire old SDFs.
    retired_rows = retire_existing_sdfs(
        ligase_root=ligase_root,
        retired_root=retired_root,
        run_stamp=run_stamp,
        dry_run=args.dry_run,
    )

    clean_active_sdf_dirs(ligase_root, dry_run=args.dry_run)

    # Rebuild new SDFs from current PDBs.
    manifest_rows, collision_rows, failure_rows = rebuild_sdfs(
        ligase_root=ligase_root,
        dry_run=args.dry_run,
    )

    if not args.dry_run:
        write_csv(MANIFEST_CSV, MANIFEST_COLUMNS, manifest_rows)
        write_csv(COLLISION_CSV, COLLISION_COLUMNS, collision_rows)
        write_csv(FAILURE_CSV, FAILURE_COLUMNS, failure_rows)

        retired_manifest = retired_root / f"Rebuild_{run_stamp}" / "Retired_SDF_Rebuild_Manifest.csv"
        retired_manifest.parent.mkdir(parents=True, exist_ok=True)
        write_csv(
            retired_manifest,
            ["Original", "Retired"],
            retired_rows,
        )

        log_payload = {
            "started_finished": now(),
            "ligase_root": str(ligase_root),
            "active_pdb_files": len(pdbs),
            "old_sdfs_retired": len(retired_rows),
            "rebuilt_instance_sdfs": len(manifest_rows),
            "group_collisions": len(collision_rows),
            "failures": len(failure_rows),
            "manifest_csv": str(MANIFEST_CSV),
            "collision_csv": str(COLLISION_CSV),
            "failure_csv": str(FAILURE_CSV),
            "retired_manifest": str(retired_manifest),
        }
        LOG_FILE.write_text(json.dumps(log_payload, indent=2) + "\n")

    print("\n=== Rebuild Summary ===")
    print(f"Old SDFs retired: {len(retired_rows)}")
    print(f"PDB-derived instance SDFs rebuilt: {len(manifest_rows)}")
    print(f"Ligase+Ligand groups with multiple PDB instances/collisions: {len(collision_rows)}")
    print(f"Failures: {len(failure_rows)}")

    if args.dry_run:
        print("\n🧪 Dry run complete. Nothing was moved/written.")
        print("Run for real with:")
        print("  python 8B_MetaDataClean.py --overwrite")
    else:
        print(f"\n✅ Written: {MANIFEST_CSV}")
        print(f"✅ Written: {COLLISION_CSV}")
        print(f"✅ Written: {FAILURE_CSV}")
        print(f"🧾 Written: {LOG_FILE}")
        print("\nNext:")
        print("  python 8_Ligand_Metadata_From_SDF_HEAVYSAFE.py --overwrite")

    if failure_rows:
        print(f"\n⚠️ Review failures in {FAILURE_CSV}")

    if collision_rows:
        print(f"\n⚠️ Review ligand-code collisions in {COLLISION_CSV}")
        print("   These are cases where one Ligase+Ligand code maps to multiple PDB-instance molecules.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())