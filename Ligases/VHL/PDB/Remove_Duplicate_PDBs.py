#!/usr/bin/env python3
"""
🧹 Remove duplicate PDB variants (_2, _3, etc.) in current directory.

Run this script *inside the ligand subdirectory* (e.g., VHL/PDB/).
It will delete any PDB files whose names end with _2.pdb, _3.pdb, etc.
Only the base ligand-bound file (e.g., 5NW1_9BH.pdb) will remain.
"""

from pathlib import Path
import re

def remove_duplicates():
    current_dir = Path.cwd()
    pdb_files = list(current_dir.glob("*.pdb"))
    pattern = re.compile(r"_(\d+)\.pdb$", re.IGNORECASE)
    removed = []

    for pdb_file in pdb_files:
        match = pattern.search(pdb_file.name)
        if match and int(match.group(1)) >= 2:
            pdb_file.unlink()
            removed.append(pdb_file.name)

    if removed:
        print("🗑️  Removed duplicate PDBs:")
        for name in removed:
            print(f"   • {name}")
    else:
        print("✅ No duplicate (_2, _3, ...) PDB files found.")

if __name__ == "__main__":
    remove_duplicates()
