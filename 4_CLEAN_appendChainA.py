#!/usr/bin/env python3
"""
🧬 Force all residues (ATOM, HETATM, ANISOU, etc.) to Chain A
across all PDB files within */PDB/ directories.
Overwrites files in place and logs all processed ligases and PDBs.
"""

from pathlib import Path

def force_chain_a(pdb_path, ligase_name, report):
    """Force all atom records in a PDB file to Chain A."""
    modified = False
    updated_lines = []

    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "ANISOU")):
                # Chain ID is column 22 (index 21)
                if line[21] != "A":
                    line = line[:21] + "A" + line[22:]
                    modified = True
            updated_lines.append(line)

    # Overwrite file in place
    pdb_path.write_text("".join(updated_lines))

    if modified:
        report.append(f"🔹 [{ligase_name}] {pdb_path.name}: Chain IDs set to A.")
    else:
        report.append(f"✔ [{ligase_name}] {pdb_path.name}: Already all Chain A.")


def main():
    root = Path.cwd()
    pdb_dirs = [p for p in root.glob("*/PDB") if p.is_dir()]
    report = []

    for pdb_dir in pdb_dirs:
        ligase_name = pdb_dir.parent.name
        report.append(f"\n=== 🧬 LIGASE: {ligase_name} ===")
        pdb_files = sorted(pdb_dir.glob("*.pdb"))
        if not pdb_files:
            report.append("⚠️  No PDB files found.")
            continue

        for pdb_file in pdb_files:
            force_chain_a(pdb_file, ligase_name, report)

    log_file = root / "Final_ForceChainA_Report.log"
    log_file.write_text("\n".join(report))
    print(f"📄 Chain A enforcement complete. Log saved to {log_file}")


if __name__ == "__main__":
    main()
