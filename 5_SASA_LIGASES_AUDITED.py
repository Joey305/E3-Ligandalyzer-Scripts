#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
5_SASA_LIGASES_AUDITED.py
===============================================================================
Purpose:
    Safer replacement for 5_SASA_LIGASES.py.

Major fixes vs old version:
    1. No silent `except: pass` for ligand MW / ligand discovery failures.
    2. Uses exact residue-name matching, not substring matching.
    3. Writes CSVs from the parent process only, avoiding parallel write races.
    4. Adds Heavy_atoms / Hydrogen_atoms / Exposed_heavy_atoms to the summary.
    5. Logs every warning/failure to explicit audit CSVs.
    6. Refuses to append to old outputs unless --overwrite is supplied.

Outputs:
    Ligand_SASA_atoms.csv
    Ligand_SASA_summary.csv
    Ligand_SASA_Audit.csv
    Ligand_SASA_failures.csv
    SASA_Report.log
    BRD_Candidates.txt

Usage:
    python 5_SASA_LIGASES_AUDITED.py --overwrite
    python 5_SASA_LIGASES_AUDITED.py --ligase-root Ligases --overwrite
    python 5_SASA_LIGASES_AUDITED.py --ligase VHL --overwrite
===============================================================================
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import multiprocessing
import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

# ---------- CONFIG ----------
DEFAULT_PROBE = 1.4
DEFAULT_THRESHOLD = 0.1
WATER_RESNAMES = {"HOH", "WAT", "DOD"}

ATOM_CSV = Path("Ligand_SASA_atoms.csv")
SUMMARY_CSV = Path("Ligand_SASA_summary.csv")
AUDIT_CSV = Path("Ligand_SASA_Audit.csv")
FAIL_CSV = Path("Ligand_SASA_failures.csv")
LOG_FILE = Path("SASA_Report.log")
BRD_REPORT = Path("BRD_Candidates.txt")

ATOM_COLUMNS = [
    "Ligase", "pdb_id", "Ligand", "Residue_ID", "Variant", "MW", "Recruiter_Class",
    "Chain", "atom_id", "exact_atom", "atom_type", "is_heavy_atom", "x", "y", "z", "Exposure_A2",
]

SUMMARY_COLUMNS = [
    # Original / downstream-compatible columns
    "Ligase", "pdb_id", "Ligand", "Residue_ID", "Variant", "MW", "Recruiter_Class",
    "Total_atoms", "Exposed_atoms", "SASA_in_complex_A2", "%Exposed", "%Buried",
    # Added validation-safe columns
    "Heavy_atoms", "Hydrogen_atoms", "Exposed_heavy_atoms", "SASA_heavy_A2", "%Exposed_Heavy", "%Buried_Heavy",
    "Ligand_Residue_Count", "Ligand_Chains", "MW_Status", "Warnings",
]

AUDIT_COLUMNS = [
    "Status", "pdb_file", "Ligase", "pdb_id", "Ligand", "Variant",
    "Total_atoms", "Heavy_atoms", "Hydrogen_atoms", "Exposed_atoms", "Exposed_heavy_atoms",
    "Ligand_Residue_Count", "Ligand_Chains", "MW", "Recruiter_Class", "MW_Status", "Warnings",
]

FAIL_COLUMNS = [
    "pdb_file", "Ligase", "pdb_id", "Ligand", "Variant", "Error", "Traceback",
]


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def backup_if_exists(path: Path) -> Optional[Path]:
    if not path.exists():
        return None
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = path.with_name(f"{path.stem}.backup_{stamp}{path.suffix}")
    shutil.copy2(path, backup)
    return backup


def prepare_outputs(overwrite: bool) -> None:
    outputs = [ATOM_CSV, SUMMARY_CSV, AUDIT_CSV, FAIL_CSV, LOG_FILE, BRD_REPORT]
    existing = [p for p in outputs if p.exists()]
    if existing and not overwrite:
        joined = "\n  ".join(str(p) for p in existing)
        raise SystemExit(
            "❌ Output files already exist. This audited version will not append into old results.\n"
            f"Existing files:\n  {joined}\n\n"
            "Run with --overwrite to back them up and regenerate clean outputs."
        )
    if overwrite:
        for p in existing:
            backup = backup_if_exists(p)
            if backup:
                p.unlink()
                print(f"🧾 Backed up {p} → {backup.name}")

    for path, cols in [(ATOM_CSV, ATOM_COLUMNS), (SUMMARY_CSV, SUMMARY_COLUMNS), (AUDIT_CSV, AUDIT_COLUMNS), (FAIL_CSV, FAIL_COLUMNS)]:
        with path.open("w", newline="") as f:
            csv.DictWriter(f, fieldnames=cols).writeheader()

    LOG_FILE.write_text(f"=== Audited SASA run started {now()} ===\n", encoding="utf-8")
    if BRD_REPORT.exists():
        BRD_REPORT.unlink()


def parse_filename(pdb_path: Path) -> Tuple[str, str, str]:
    """
    Parse filenames like:
        7HLV_A28_2.pdb -> pdb_id=7HLV, ligand=A28, variant=2
        8GCG_A1A_3.pdb -> pdb_id=8GCG, ligand=A1A, variant=3
        5XYZ_ABC.pdb   -> pdb_id=5XYZ, ligand=ABC, variant=1
    """
    parts = pdb_path.stem.split("_")
    pdb_id = parts[0] if parts else pdb_path.stem
    ligand = parts[1] if len(parts) >= 2 else "UNK"
    variant = parts[2] if len(parts) >= 3 and parts[2].isdigit() else "1"
    return pdb_id.upper(), ligand.upper(), variant


def infer_ligase(pdb_path: Path) -> str:
    parts = list(pdb_path.parts)
    if "Ligases" in parts:
        i = parts.index("Ligases")
        if i + 1 < len(parts):
            return parts[i + 1]
    # fallback for paths like Ligases/VHL/PDB/x.pdb or VHL/PDB/x.pdb
    if pdb_path.parent.name.upper() == "PDB":
        return pdb_path.parent.parent.name
    return pdb_path.parent.name


def residue_name_from_pdb_line(line: str) -> str:
    return line[17:20].strip().upper()


def atom_element(atom) -> str:
    elem = getattr(atom, "element", "") or ""
    elem = str(elem).strip().upper()
    if elem:
        return elem
    # Fallback from atom name, e.g. " C12" or "H1"
    name = atom.get_name().strip().upper()
    letters = "".join(c for c in name if c.isalpha())
    return letters[0] if letters else ""


def is_hydrogen_atom(atom) -> bool:
    return atom_element(atom) == "H"


def pdb_ligand_lines_exact(pdb_path: Path, ligand_tag: str) -> List[str]:
    lines = []
    with pdb_path.open(errors="replace") as f:
        for line in f:
            if line.startswith("HETATM") and residue_name_from_pdb_line(line) == ligand_tag.upper():
                lines.append(line)
    return lines


def pdb_all_nonwater_hetatm_by_resname(pdb_path: Path) -> Dict[str, List[str]]:
    grouped: Dict[str, List[str]] = defaultdict(list)
    with pdb_path.open(errors="replace") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            resname = residue_name_from_pdb_line(line)
            if resname in WATER_RESNAMES or not resname:
                continue
            grouped[resname].append(line)
    return grouped


def mol_from_pdb_lines(lines: List[str]) -> Optional[Chem.Mol]:
    if not lines:
        return None
    pdb_block = "".join(lines)
    try:
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=True)
        if mol:
            return mol
    except Exception:
        pass
    try:
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
        if mol:
            return mol
    except Exception:
        pass
    return None


def compute_ligand_mw(pdb_path: Path, ligand_tag: str) -> Tuple[Optional[float], str]:
    """Return (MW, status). No silent failure."""
    lines = pdb_ligand_lines_exact(pdb_path, ligand_tag)
    if not lines:
        return None, "NO_EXACT_HETATM_LINES"
    mol = mol_from_pdb_lines(lines)
    if mol is None:
        return None, "RDKIT_PDB_MOL_FAILED"
    try:
        return round(float(Descriptors.MolWt(mol)), 3), "OK"
    except Exception as e:
        return None, f"MW_DESCRIPTOR_FAILED: {type(e).__name__}: {e}"


def get_all_ligands_for_brd(pdb_path: Path) -> Tuple[Dict[str, float], List[str]]:
    """Return {resname: MW} for all non-water HETATM residues with exact resname grouping."""
    warnings: List[str] = []
    ligands: Dict[str, float] = {}
    grouped = pdb_all_nonwater_hetatm_by_resname(pdb_path)
    for resname, lines in grouped.items():
        mol = mol_from_pdb_lines(lines)
        if mol is None:
            warnings.append(f"BRD_MW_PARSE_FAILED:{resname}")
            continue
        try:
            ligands[resname] = round(float(Descriptors.MolWt(mol)), 3)
        except Exception as e:
            warnings.append(f"BRD_MW_DESCRIPTOR_FAILED:{resname}:{type(e).__name__}")
    return ligands, warnings


def classify_recruiter(mw: Optional[float], brd_flag: bool) -> str:
    if brd_flag:
        return "BRD Recruiter"
    if mw is None or (isinstance(mw, float) and math.isnan(mw)):
        return "Unknown"
    if mw < 200:
        return "Fragment Recruiter"
    if mw <= 500:
        return "Drug-like Recruiter"
    return "Peptide-like Recruiter"


def analyze_pdb(pdb_path_str: str, probe_radius: float, threshold: float, include_zero_atoms: bool) -> Dict[str, object]:
    pdb_path = Path(pdb_path_str)
    ligase = infer_ligase(pdb_path)
    pdb_id, ligand_tag, variant = parse_filename(pdb_path)
    warnings: List[str] = []

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", str(pdb_path))

        # Remove waters before SASA
        for model in structure:
            for chain in model:
                residues_to_remove = [r for r in chain if r.resname.strip().upper() in WATER_RESNAMES]
                for res in residues_to_remove:
                    chain.detach_child(res.id)

        sr = ShrakeRupley(probe_radius=probe_radius)
        sr.compute(structure, level="A")

        # MW + BRD
        all_ligs, brd_warnings = get_all_ligands_for_brd(pdb_path)
        warnings.extend(brd_warnings)
        small_fragments = [mw for mw in all_ligs.values() if mw < 200]
        brd_flag = len(small_fragments) >= 4
        mw, mw_status = compute_ligand_mw(pdb_path, ligand_tag)
        if mw_status != "OK":
            warnings.append(mw_status)
        recruiter_class = classify_recruiter(mw, brd_flag)

        # Collect matching ligand residues by exact residue name
        ligand_atoms: List[Dict[str, object]] = []
        residue_ids: List[str] = []
        chains: List[str] = []
        total_atoms = 0
        heavy_atoms = 0
        hydrogen_atoms = 0
        exposed_atoms = 0
        exposed_heavy_atoms = 0
        sasa_total = 0.0
        sasa_heavy = 0.0

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.resname.strip().upper() != ligand_tag:
                        continue
                    resnum = residue.id[1]
                    icode = residue.id[2].strip() if len(residue.id) > 2 else ""
                    residue_id = f"{chain.id}:{resnum}{icode}"
                    residue_ids.append(residue_id)
                    chains.append(str(chain.id))
                    for atom in residue.get_atoms():
                        elem = atom_element(atom)
                        is_heavy = elem != "H"
                        sasa_val = float(getattr(atom, "sasa", 0.0) or 0.0)

                        total_atoms += 1
                        if is_heavy:
                            heavy_atoms += 1
                        else:
                            hydrogen_atoms += 1

                        if sasa_val > threshold:
                            exposed_atoms += 1
                            sasa_total += sasa_val
                            if is_heavy:
                                exposed_heavy_atoms += 1
                                sasa_heavy += sasa_val

                        if include_zero_atoms or sasa_val > threshold:
                            ligand_atoms.append({
                                "Ligase": ligase,
                                "pdb_id": pdb_id,
                                "Ligand": ligand_tag,
                                "Residue_ID": residue_id,
                                "Variant": variant,
                                "MW": mw if mw is not None else "NA",
                                "Recruiter_Class": recruiter_class,
                                "Chain": chain.id,
                                "atom_id": atom.serial_number,
                                "exact_atom": atom.name,
                                "atom_type": elem,
                                "is_heavy_atom": bool(is_heavy),
                                "x": round(float(atom.coord[0]), 3),
                                "y": round(float(atom.coord[1]), 3),
                                "z": round(float(atom.coord[2]), 3),
                                "Exposure_A2": round(sasa_val, 3),
                            })

        if total_atoms == 0:
            raise ValueError(f"No atoms found for exact ligand residue '{ligand_tag}' in {pdb_path.name}")

        percent_exposed = exposed_atoms / total_atoms if total_atoms else 0.0
        percent_exposed_heavy = exposed_heavy_atoms / heavy_atoms if heavy_atoms else 0.0

        summary = {
            "Ligase": ligase,
            "pdb_id": pdb_id,
            "Ligand": ligand_tag,
            "Residue_ID": ";".join(sorted(set(residue_ids))) if residue_ids else "NA",
            "Variant": variant,
            "MW": mw if mw is not None else "NA",
            "Recruiter_Class": recruiter_class,
            "Total_atoms": total_atoms,
            "Exposed_atoms": exposed_atoms,
            "SASA_in_complex_A2": round(sasa_total, 3),
            "%Exposed": round(percent_exposed, 3),
            "%Buried": round(1.0 - percent_exposed, 3),
            "Heavy_atoms": heavy_atoms,
            "Hydrogen_atoms": hydrogen_atoms,
            "Exposed_heavy_atoms": exposed_heavy_atoms,
            "SASA_heavy_A2": round(sasa_heavy, 3),
            "%Exposed_Heavy": round(percent_exposed_heavy, 3),
            "%Buried_Heavy": round(1.0 - percent_exposed_heavy, 3),
            "Ligand_Residue_Count": len(set(residue_ids)),
            "Ligand_Chains": ";".join(sorted(set(chains))),
            "MW_Status": mw_status,
            "Warnings": "; ".join(sorted(set(warnings))),
        }

        audit = {
            "Status": "OK" if not warnings else "CHECK",
            "pdb_file": str(pdb_path),
            "Ligase": ligase,
            "pdb_id": pdb_id,
            "Ligand": ligand_tag,
            "Variant": variant,
            "Total_atoms": total_atoms,
            "Heavy_atoms": heavy_atoms,
            "Hydrogen_atoms": hydrogen_atoms,
            "Exposed_atoms": exposed_atoms,
            "Exposed_heavy_atoms": exposed_heavy_atoms,
            "Ligand_Residue_Count": len(set(residue_ids)),
            "Ligand_Chains": ";".join(sorted(set(chains))),
            "MW": mw if mw is not None else "NA",
            "Recruiter_Class": recruiter_class,
            "MW_Status": mw_status,
            "Warnings": summary["Warnings"],
        }

        brd_info = None
        if brd_flag:
            brd_info = {
                "pdb_id": pdb_id,
                "ligase": ligase,
                "fragment_count": len(small_fragments),
                "ligands": ", ".join([f"{k}({v} Da)" for k, v in sorted(all_ligs.items())]),
            }

        msg = (
            f"✅ {pdb_path.name}: {exposed_atoms}/{total_atoms} atoms exposed "
            f"({percent_exposed:.2f}); heavy {exposed_heavy_atoms}/{heavy_atoms} ({percent_exposed_heavy:.2f})"
        )
        if warnings:
            msg += f" ⚠️ {'; '.join(sorted(set(warnings)))}"

        return {
            "ok": True,
            "pdb_file": str(pdb_path),
            "atom_rows": ligand_atoms,
            "summary": summary,
            "audit": audit,
            "brd_info": brd_info,
            "message": msg,
        }

    except Exception as e:
        tb = traceback.format_exc()
        fail = {
            "pdb_file": str(pdb_path),
            "Ligase": ligase,
            "pdb_id": pdb_id,
            "Ligand": ligand_tag,
            "Variant": variant,
            "Error": f"{type(e).__name__}: {e}",
            "Traceback": tb,
        }
        return {
            "ok": False,
            "pdb_file": str(pdb_path),
            "failure": fail,
            "message": f"❌ {pdb_path.name} failed: {type(e).__name__}: {e}",
        }


def write_dict_rows(path: Path, fieldnames: List[str], rows: Iterable[Dict[str, object]]) -> int:
    count = 0
    with path.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def write_dict_row(path: Path, fieldnames: List[str], row: Dict[str, object]) -> None:
    with path.open("a", newline="") as f:
        csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore").writerow(row)


def clean_and_sort_outputs() -> None:
    if SUMMARY_CSV.exists():
        df = pd.read_csv(SUMMARY_CSV)
        if not df.empty and {"Ligase", "%Exposed_Heavy"}.issubset(df.columns):
            df = df.sort_values(by=["Ligase", "%Exposed_Heavy", "%Exposed"], ascending=[True, False, False])
        elif not df.empty and {"Ligase", "%Exposed"}.issubset(df.columns):
            df = df.sort_values(by=["Ligase", "%Exposed"], ascending=[True, False])
        df.to_csv(SUMMARY_CSV, index=False)

    if ATOM_CSV.exists():
        df = pd.read_csv(ATOM_CSV)
        if not df.empty and {"Ligase", "pdb_id", "Ligand", "Exposure_A2"}.issubset(df.columns):
            df = df.sort_values(by=["Ligase", "pdb_id", "Ligand", "Exposure_A2"], ascending=[True, True, True, False])
        df.to_csv(ATOM_CSV, index=False)

    if AUDIT_CSV.exists():
        df = pd.read_csv(AUDIT_CSV)
        if not df.empty and {"Status", "Ligase", "pdb_id", "Ligand"}.issubset(df.columns):
            df = df.sort_values(by=["Status", "Ligase", "pdb_id", "Ligand"])
        df.to_csv(AUDIT_CSV, index=False)


def discover_pdbs(ligase_root: Path, ligase: Optional[str]) -> List[Path]:
    if ligase:
        candidates = [ligase_root / ligase, Path(ligase)]
        roots = [p for p in candidates if p.exists()]
        if not roots:
            raise FileNotFoundError(f"Could not find ligase folder '{ligase}' under {ligase_root} or current path.")
        pdbs = []
        for root in roots:
            pdbs.extend(root.rglob("PDB/*.pdb"))
        return sorted(set(pdbs))
    return sorted(ligase_root.rglob("PDB/*.pdb"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute ligand SASA with strict auditing and heavy-atom validation.")
    parser.add_argument("--ligase-root", default="Ligases", help="Root folder containing Ligases/<Ligase>/PDB/*.pdb")
    parser.add_argument("--ligase", default=None, help="Process only one ligase folder, e.g. VHL")
    parser.add_argument("--probe", type=float, default=DEFAULT_PROBE, help="Probe radius for SASA calculation")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, help="Exposure threshold for atom exposure")
    parser.add_argument("--workers", type=int, default=max(1, multiprocessing.cpu_count() - 1), help="Number of worker processes")
    parser.add_argument("--overwrite", action="store_true", help="Backup and replace existing output CSV/log files")
    parser.add_argument(
        "--exclude-zero-atoms",
        action="store_true",
        help="Drop zero-exposure atoms from atom-level CSV. Default keeps all ligand atoms."
    )
    args = parser.parse_args()

    ligase_root = Path(args.ligase_root)
    all_pdbs = discover_pdbs(ligase_root, args.ligase)
    total_files = len(all_pdbs)
    if total_files == 0:
        raise SystemExit(f"❌ No PDB files found under {ligase_root}")

    print(f"🔍 Scanning {total_files} PDB files with probe radius {args.probe} Å")
    print(f"🧪 Exact ligand matching + heavy-atom audit enabled")
    print(f"🧠 Using {args.workers} worker(s)\n")

    prepare_outputs(args.overwrite)

    ok_count = 0
    check_count = 0
    fail_count = 0
    atom_row_count = 0
    brd_hits = 0

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(analyze_pdb, str(p), args.probe, args.threshold, not args.exclude_zero_atoms): p
            for p in all_pdbs
        }
        for i, fut in enumerate(as_completed(futures), 1):
            result = fut.result()
            print(result["message"])
            with LOG_FILE.open("a", encoding="utf-8") as log:
                log.write(str(result["message"]) + "\n")

            if result.get("ok"):
                ok_count += 1
                summary = result["summary"]
                audit = result["audit"]
                if audit.get("Status") == "CHECK":
                    check_count += 1
                atom_row_count += write_dict_rows(ATOM_CSV, ATOM_COLUMNS, result.get("atom_rows", []))
                write_dict_row(SUMMARY_CSV, SUMMARY_COLUMNS, summary)
                write_dict_row(AUDIT_CSV, AUDIT_COLUMNS, audit)
                brd = result.get("brd_info")
                if brd:
                    brd_hits += 1
                    with BRD_REPORT.open("a", encoding="utf-8") as log:
                        log.write(
                            f"{brd['ligase']}: {brd['pdb_id']} "
                            f"({brd['fragment_count']} fragments) → {brd['ligands']}\n"
                        )
            else:
                fail_count += 1
                write_dict_row(FAIL_CSV, FAIL_COLUMNS, result["failure"])

            if i % 10 == 0 or i == total_files:
                print(f"📊 Progress: {i}/{total_files} processed...")

    print("\n🧹 Sorting outputs...")
    clean_and_sort_outputs()

    summary_lines = [
        "\n=== Audited SASA summary ===",
        f"PDB files discovered: {total_files}",
        f"Successful PDBs: {ok_count}",
        f"Successful with warnings/checks: {check_count}",
        f"Failed PDBs: {fail_count}",
        f"Atom-level rows written: {atom_row_count}",
        f"BRD candidate structures: {brd_hits}",
        f"Finished: {now()}",
    ]
    print("\n".join(summary_lines))
    with LOG_FILE.open("a", encoding="utf-8") as log:
        log.write("\n".join(summary_lines) + "\n")

    if fail_count:
        print(f"\n⚠️ Review failures in {FAIL_CSV}")
    if check_count:
        print(f"⚠️ Review warnings/checks in {AUDIT_CSV}")
    print("✅ Done.")


if __name__ == "__main__":
    main()
