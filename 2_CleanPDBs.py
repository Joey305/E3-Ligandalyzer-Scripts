#!/usr/bin/env python3
"""
cleanup_pdbs.py — removes intermediate .pdb files after LigaseBatchFetch
and logs a per-ligase summary.

Removes:
  • All .pdb files inside CIF/ directories
  • All plain ####.pdb files inside PDB/ directories
Keeps:
  • All ligand-bound PDBs named ####_LIG.pdb
Logs:
  • Every action (deleted / kept) and per-ligase summary
"""

import csv
from pathlib import Path
from datetime import datetime

def cleanup_ligase_dirs(base_dir="."):
    base = Path(base_dir)
    log_dir = base / "logs"
    log_dir.mkdir(exist_ok=True)

    log_path = log_dir / "PDB_FILE_REMOVAL.log"

    with open(log_path, "w", newline="") as logfile:
        writer = csv.writer(logfile)
        writer.writerow(["timestamp", "ligase", "directory", "filename", "action", "status"])

        ligase_dirs = [d for d in base.iterdir() if d.is_dir() and d.name not in {"logs", "__pycache__"}]
        summary = []  # store per-ligase counts

        for ligase_dir in ligase_dirs:
            ligase_name = ligase_dir.name
            pdb_kept = 0
            ligands_found = 0

            # 1️⃣ Remove all .pdb files in CIF/
            cif_dir = ligase_dir / "CIF"
            if cif_dir.exists():
                for pdb_file in cif_dir.glob("*.pdb"):
                    try:
                        pdb_file.unlink()
                        writer.writerow([
                            datetime.now().isoformat(),
                            ligase_name, "CIF", pdb_file.name, "delete", "success"
                        ])
                    except Exception as e:
                        writer.writerow([
                            datetime.now().isoformat(),
                            ligase_name, "CIF", pdb_file.name, "delete", f"error: {e}"
                        ])

            # 2️⃣ Remove unliganded ####.pdb in PDB/
            pdb_dir = ligase_dir / "PDB"
            if pdb_dir.exists():
                for pdb_file in pdb_dir.glob("*.pdb"):
                    name = pdb_file.stem
                    if "_" not in name:
                        try:
                            pdb_file.unlink()
                            writer.writerow([
                                datetime.now().isoformat(),
                                ligase_name, "PDB", pdb_file.name, "delete", "success"
                            ])
                        except Exception as e:
                            writer.writerow([
                                datetime.now().isoformat(),
                                ligase_name, "PDB", pdb_file.name, "delete", f"error: {e}"
                            ])
                    else:
                        pdb_kept += 1
                        writer.writerow([
                            datetime.now().isoformat(),
                            ligase_name, "PDB", pdb_file.name, "keep", "ligand_bound"
                        ])

            # 3️⃣ Count unique ligands in SDF/
            sdf_dir = ligase_dir / "SDF"
            if sdf_dir.exists():
                ligands_found = len(list(sdf_dir.glob("*.sdf")))

            summary.append({
                "ligase": ligase_name,
                "pdb_kept": pdb_kept,
                "ligands_found": ligands_found
            })

        # 4️⃣ Append summary section
        writer.writerow([])
        writer.writerow(["──────────── Summary per Ligase ────────────"])
        writer.writerow(["ligase", "pdb_files_kept", "unique_ligands_found"])
        for s in summary:
            writer.writerow([s["ligase"], s["pdb_kept"], s["ligands_found"]])

    # 5️⃣ Print summary to console too
    print("\n📊 Cleanup Summary:")
    for s in summary:
        print(f"   • {s['ligase']}: {s['pdb_kept']} kept PDBs, {s['ligands_found']} ligands")

    print(f"\n✅ Cleanup complete! Detailed log saved to {log_path}")


if __name__ == "__main__":
    cleanup_ligase_dirs(".")
