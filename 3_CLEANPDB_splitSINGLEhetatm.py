#!/usr/bin/env python3
"""

🧩 Detects and splits PDBs that contain multiple ligand copies
into separate _1.pdb, _2.pdb variants (keeping proteins intact).
Now logs with ligase folder context (e.g., PARKIN, CRBN, etc.).
"""

from pathlib import Path
from collections import defaultdict

def parse_pdb(pdb_path):
    """Parse a PDB into atoms grouped by residue ID (chain+resnum)."""
    residues = defaultdict(list)
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip() or "A"
                resid = line[22:26].strip()
                resid_full = f"{chain}{resid}"
                residues[resid_full].append(line)
    return residues


def split_ligands(pdb_path, ligase_name, report):
    ligand_code = pdb_path.stem.split("_")[-1].upper()[:3]
    residues = parse_pdb(pdb_path)

    # Identify which residue IDs contain the ligand
    ligand_resids = [rid for rid, lines in residues.items() if any(
        ln.startswith("HETATM") and ln[17:20].strip().upper() == ligand_code for ln in lines
    )]

    if len(ligand_resids) <= 1:
        report.append(f"✔ [{ligase_name}] {pdb_path.name}: Single ligand ({ligand_code}) → no split needed.")
        return

    # Collect all protein ATOM lines once
    protein_lines = [ln for rid, lines in residues.items() for ln in lines if ln.startswith("ATOM")]

    for i, resid in enumerate(ligand_resids, start=1):
        out_file = pdb_path.with_name(f"{pdb_path.stem}_{i}.pdb")
        with open(out_file, "w") as f:
            for ln in protein_lines:
                f.write(ln)
            for ln in residues[resid]:
                if ln.startswith("HETATM"):
                    f.write(ln)
        report.append(f"🔸 [{ligase_name}] {pdb_path.name}: Multiple {ligand_code} residues → wrote {out_file.name} ({resid})")


def main():
    root = Path.cwd()
    pdb_dirs = [p for p in root.glob("*/PDB") if p.is_dir()]
    report = []

    for pdb_dir in pdb_dirs:
        ligase_name = pdb_dir.parent.name
        report.append(f"\n=== 🧬 LIGASE: {ligase_name} ===")
        for pdb_file in sorted(pdb_dir.glob("*.pdb")):
            split_ligands(pdb_file, ligase_name, report)

    log_file = root / "Final_Split_Report.log"
    log_file.write_text("\n".join(report))
    print(f"📄 Split complete. Detailed log saved to {log_file}")


if __name__ == "__main__":
    main()
