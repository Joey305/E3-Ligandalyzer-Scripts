#!/usr/bin/env python3
"""
LigaseBatchFetch.py — batch fetcher for ligand-bound E3 ligases from RCSB PDB.

For each ligase entry in the LIGASES dict:
    - Searches RCSB using multiple aliases
    - Aligns all hits to reference chain from template PDB
    - Extracts most similar chain (≥ ~70% identity)
    - Saves filtered PDBs (ligase + bound organic ligand)
    - Excludes waters, ions, and small fragments
    - Saves SDFs, logs, and a compressed ZIP for each ligase

Output structure:
  LIGASE_NAME/
    CIF/
    PDB/
    SDF/
  LIGASE_NAME_log.csv
  LIGASE_NAME_YYYY-MM-DD.zip
"""

import os, re, csv, zipfile, datetime, time, requests
from pathlib import Path
from typing import Dict, List, Tuple, Set
from Bio.PDB import PDBIO, NeighborSearch, Select
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner, substitution_matrices
from requests.adapters import HTTPAdapter, Retry
import requests, gemmi
from typing import List
import subprocess, tempfile, shutil



# ────────────────────────────────────────────────
# CONFIG
DATE_TAG = datetime.date.today().isoformat()
MW_MIN, MW_MAX = 150, 700
DISTANCE_CUTOFF = 6.0  # Å
SEQUENCE_IDENTITY_THRESHOLD = 0.3
# EXCLUDE = {"HOH","NA","K","ZN","MG","CL","SO4","GOL","EDO","TRS","MES","HEP","CAC","ACT","PEG","DMS","PO4" "CIT"}


# ────────────────────────────────────────────────
# Big exclusion lists (names only). These cover:
# - Ions & metals
# - Buffers/solvents/cryoprotectants/PEGs/detergents
# - Nucleotides & common cofactors
# - Sugars / glycans
# - Common PTM residue codes / PLP adducts
# Add anything project-specific to EXTRA_EXCLUDE below.
IONS_METALS = {
    # alkali/alkaline/transition/post-transition/lanthanides (common)
    "H", "D", "NA", "K", "LI", "RB", "CS",
    "CA", "MG", "SR", "BA",
    "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG",
    "AL", "GA", "IN", "TL", "SN", "PB",
    "AG", "AU", "PT",
    "YB", "LA", "CE", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "LU",
    # halides & small anions
    "CL", "BR", "IOD", "I", "F",
    "SCN", "CN", "OCN",
    "NO3", "NO2", "SO3", "SO4", "PO4", "HPO", "CO3",
}

BUFFERS_SOLVENTS = {
    # water / heavy water
    "HOH", "H2O", "DOD",
    # common organics
    "GOL", "EDO", "PGO", "MPD", "IPA", "EOH", "FMT", "ACT", "ACE", "ACO",
    "CIT", "TLA", "MLE", "OXA", "GLY", "ALA", "DAL",
    # zwitterion buffers / good-bad buffers
    "TRS", "TES", "MES", "HEP", "HEPES", "MOP", "MOPS", "PIP", "PIPES", "TAPS", "BIC",
    # PEGs & fragments
    "PEG", "P33", "P40", "PG4", "PG5", "PG6", "EPE", "PE4", "1PE",
    # detergents & lipids (broad brush)
    "BOG", "B0G", "LDAO", "CYM", "CHAPS", "CHS", "CLR", "CHL",
    # reducing agents
    "BME", "DTT", "TCEP",
    # other lab reagents
    "EDTA", "ATP", "ADP", "AMP-PNP", "AMPPCP",
}

NUCLEOTIDES_COFAC = {
    # nucleotides & analogs
    "ATP", "ADP", "AMP", "GTP", "GDP", "GNP", "CTP", "UTP", "TTP", "ITP", "IDP", "IMP",
    "SAM", "SAH", "COA", "FAD", "FMN", "NAD", "NADH", "NAP",
    # hemes/porphyrins
    "HEM", "HEA", "HEB", "HEC", "PHE", "BCL", "BLA",
    # PLP & adducts
    "PLP", "PMP", "PLG", "LLP",  # PLG = PLP-glycine adduct; LLP = PLP-Lys adduct
}

SUGARS_GLYCANS = {
    "NAG", "NDG", "BMA", "MAN", "GAL", "GLC", "FUC", "SIA", "NANA", "NGC", "A2G", "BGC",
    "FCA", "LMT", "MNT", "XYS", "LAC", "MAL", "SUC", "TRE",
}

PTM_RESIDUES = {
    # frequent amino-acid modifications in HET form
    "MSE", "SEC", "PYL",
    "SEP", "TPO", "PTR",       # phospho-Ser/Thr/Tyr
    "KCX", "CME",              # NZ-carboxylysine; S-carboxymethylcys
    "MLY", "MLZ", "M3L", "M2L",# mono/di/tri-methyl-Lys
    "HYP", "HIC",              # hydroxyproline, 4-imidazolone-5-propionic acid (rare)
    "PCA", "PYR",              # pyroglutamate variants
    "CSO", "CSX",              # oxidized Cys
}

# your original minimal list for safety
LEGACY_EXCLUDE = {"GOL","EDO","TRS","MES","HEP","CAC","ACT","PEG","DMS","PO4"}


MODDED_RESIDUES = {
    # --- Cysteine variants / crosslinks ---
    "CSD",  # Cys-sulfinic acid
    "OCS",  # Cys-serine thioether
    "C5C",  # Cys–Cys disulfide mimic
    "C6C",  # Cys–Cys variant
    "CAS",  # S-(dimethylarsenic)cysteine
    "SAC",  # S-acetylcysteine
    "SMC",  # S-methylcysteine variant

    # --- Tyrosine / Serine / Lysine variants ---
    "TYS",  # Sulfotyrosine
    "ALY",  # N6-acetyllysine
    "M3L",  # (keep if you remove from PTM_RESIDUES)
    "NLE",  # Norleucine
    "DPR",  # Dehydroproline

    # --- Halogenated or labeled residues ---
    "FME",  # N-formylmethionine
    "YCM",  # Diiodotyrosine or crosslinked Tyr
    "YOF",  # 3-fluorotyrosine variant

    # --- Unnatural analogs / adducts ---
    "HIC",  # hydroxyisoleucine (if not already caught)
    "HSE",  # Homoserine
    "ABA",  # Alpha-aminobutyric acid
    "MLZ",  # methionine sulfone
    "MLE",  # oxidized or extended leucine variant
}


PEPTIDE_FRAGMENTS = {
    "ACE",  # Acetyl (N-terminal cap)
    "NH2", "NME",  # amidated/methylated C-terminal caps
    "FOR",  # Formyl cap
    "GLY", "ALA", "VAL", "LEU", "ILE",  # free amino acids (non-polymeric)
    "PRO", "SER", "THR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG",
    "TYR", "PHE", "TRP", "MET", "CYS", "HIS",  # free or capping residues
    "CYA", "CYX", "CYO",  # cysteine dimer/crosslink fragments
}

SMALL_ORGANICS_LAB_REAGENTS = {
    "DMS", "DMF", "DMSO", "MEO", "MPT", "BUT", "TBUT", "MCH", "HEX",
    "URE", "ACN", "BEN", "TOL", "MTB", "TRI", "XYL", "CRY", "HCL",
}


PROSTHETIC_GROUPS = {
    "HEM", "HEA", "HEB", "HEC",  # Heme groups
    "BCL", "BLA",                # Chlorophyll/bacteriochlorophyll
    "FAD", "FMN", "NAD", "NADH", "NAP",
    "COA", "SAM", "SAH", "TPQ", "DHE",
}


ANTIBIOTICS_STABILIZERS = {
    "KAN", "GEN", "NEO", "AMK", "STR", "CHL", "ERY",  # antibiotics
    "TRI", "TRS", "C12", "C14", "C16", "C18", "BCL",  # detergents/lipids
}

ENGINEERED_TAGS = {
    "HIS", "SER", "GLY", "MET", "GFP", "MBP", "SUMO", "GST",  # purification fusions
    "LBT", "CRO", "SYG", "CRH",  # lanthanide-binding tags, chaperones
}

ARTIFACT_LIGANDS = {
    "UNL",  # Unknown ligand
    "XXX",  # Placeholder ligand
    "XXX1", "XXX2",  # old-style dummy entries
    "STL",  # stain or artifact
    "UNX",  # unknown atom type
    "UNK",  # generic unknown residue
    "MOL",  # generic molecule
}


DRUGLIKE_CONTROLS = {
    # intentionally kept “known drugs” for positive controls
    "IMN", "IMD", "ATP", "NDP", "MTX", "ADN", "NAG", "SO4",
}

# project-specific catchall — add anything you dislike here
EXTRA_EXCLUDE = set()

# unified list used everywhere
EXCLUDE_ALL = (
    IONS_METALS
    | BUFFERS_SOLVENTS
    | NUCLEOTIDES_COFAC
    | SUGARS_GLYCANS
    | PTM_RESIDUES
    | MODDED_RESIDUES
    | PEPTIDE_FRAGMENTS
    | SMALL_ORGANICS_LAB_REAGENTS
    | PROSTHETIC_GROUPS
    | ANTIBIOTICS_STABILIZERS
    | ARTIFACT_LIGANDS
    | ENGINEERED_TAGS
    | LEGACY_EXCLUDE
    | EXTRA_EXCLUDE
)





# Ligase definitions — fill in as needed
LIGASES = {




    # "DCAF15": {
    #     "template_pdb": "6SJ7",
    #     "template_chain": "A",
    #     "search_terms": ["DCAF15"]
    # },
    # "KEAP1": {
    #     "template_pdb": "6ZEY",
    #     "template_chain": "A",
    #     "search_terms": ["KEAP1",]
    # },
    # "DCAF16": {
    #     "template_pdb": "8OV6",
    #     "template_chain": "B",
    #     "search_terms": ["DCAF16"]
    # },
    # "DCAF11": {
    #     "template_pdb": "8FZT",
    #     "template_chain": "A",
    #     "search_terms": ["DCAF11", "piperlongumine"]
    # },
    # "RNF114": {
    #     "template_pdb": "8E3M",
    #     "template_chain": "A",
    #     "search_terms": ["RNF114", "nimbolide"]
    # # },
    # "MDM2": {
    #     "template_pdb": "4LWV",
    #     "template_chain": "A",
    #     "search_terms": ["MDM2"]
    # },
    # "XIAP": {
    #     "template_pdb": "4KMP",
    #     "template_chain": "A",
    #     "search_terms": ["XIAP", "BIRC4", "SMAC mimetic"]
    # },
    # "cIAP1": {
    #     "template_pdb": "3MUP",
    #     "template_chain": "A",
    #     "search_terms": ["cIAP1", "BIRC2", "SMAC mimetic"]
    # },
    # "cIAP2": {
    #     "template_pdb": "3M0D",
    #     "template_chain": "D",
    #     "search_terms": ["cIAP2", "BIRC3", "SMAC mimetic"]
    # },
    # "RNF4": {
    #     "template_pdb": "4PPE",
    #     "template_chain": "A",
    #     "search_terms": ["RNF4"]
    # },
    # "TRIM21": {
    #     "template_pdb": "9Q9O",
    #     "template_chain": "A",
    #     "search_terms": ["TRIM21", "E3 UBIQUITIN-PROTEIN LIGASE TRIM21"]
    # },
    # "DCAF1": {
    #     "template_pdb": "7UFV",
    #     "template_chain": "A",
    #     "search_terms": ["DCAF1", "VPRBP", "Vpr-binding protein"]
    # },
    # "FBXO7": {
    #     "template_pdb": "7NUL",
    #     "template_chain": "A",
    #     "search_terms": ["FBXO7", "SCF complex"]
    # },
    # "FBXW7": {
    #     "template_pdb": "4A5M",
    #     "template_chain": "A",
    #     "search_terms": ["FBXW7", "SCF substrate adaptor"]
    # },
    # "SOCS2": {
    #     "template_pdb": "6I5J",
    #     "template_chain": "A",
    #     "search_terms": ["SOCS2", "CRL5"]
    # },
    # "CUL3_SPOP": {
    #     "template_pdb": "4HS2",
    #     "template_chain": "A",
    #     "search_terms": ["CUL3", "SPOP", "BTB complex"]
    # },
    # "RNF126": {
    #     "template_pdb": "6MEB",
    #     "template_chain": "A",
    #     "search_terms": ["RNF126"]
    # },
    # "RNF8": {
    #     "template_pdb": "4WHV",
    #     "template_chain": "C",
    #     "search_terms": ["RNF8", "E3 ubiquitin-protein ligase RNF8", "Ubiquitin ligase protein RNF8"]
    # },
    # "PARKIN": {
    #     "template_pdb": "5CAW",
    #     "template_chain": "A",
    #     "search_terms": ["PARKIN", "PRKN"]
    # },
    # "CHIP": {
    #     "template_pdb": "7TB1",
    #     "template_chain": "A",
    #     "search_terms": ["E3 ubiquitin-protein ligase CHIP", "STUB1", "Hsc70" ]
    # },
    # "HUWE1": {
    #     "template_pdb": "7NH3",
    #     "template_chain": "A",
    #     "search_terms": ["HUWE1", "MYC"]
    # },
    # "NEDD4L": {
    #     "template_pdb": "3JVZ",
    #     "template_chain": "C",
    #     "search_terms": ["NEDD4L"]
    # },
    # "RNF43": {
    #     "template_pdb": "4C9V",
    #     "template_chain": "A",
    #     "search_terms": ["RNF43"]
    # },


    # "CRBN": {
    #     "template_pdb": "4CI3",
    #     "template_chain": "B",
    #     "search_terms": ["Cereblon", "CRBN", "Protein cereblon"]
    # },

    "VHL": {
        "template_pdb": "4W9H",
        "template_chain": "C",
        "search_terms": ["VHL", "von Hippel-Lindau", "VHL protein"]
    },

    #  "Cbl-b": {
    #     "template_pdb": "9FQH",
    #     "template_chain": "A",
    #     "search_terms": ["E3 ubiquitin-protein ligase CBL-B", "CBL-B"]
    # },
    
    # "GID4": {
    #     "template_pdb": "7U3E",
    #     "template_chain": "A",
    #     "search_terms": ["GID4", "Glucose-induced degradation protein 4 homolog"]
    # },
}

# ────────────────────────────────────────────────
# HTTP utils
# ────────────────────────────────────────────────
def make_session():
    s = requests.Session()
    retries = Retry(total=5, backoff_factor=0.5,
                    status_forcelist=[429,500,502,503,504],
                    allowed_methods=["GET","POST"])
    s.mount("https://", HTTPAdapter(max_retries=retries))
    s.headers.update({"User-Agent":"LigaseBatchFetch/1.1"})
    return s
SESSION = make_session()


def http_get(url, timeout=10):
    return requests.get(url, timeout=timeout)

def http_post(url, **kw): kw.setdefault("timeout",60); return SESSION.post(url, **kw)

# ────────────────────────────────────────────────
def seq_identity_local(seq1, seq2):
    if not seq1 or not seq2: return 0.0
    aligner = PairwiseAligner(); aligner.mode="local"
    aligner.substitution_matrix=substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score=-10; aligner.extend_gap_score=-0.5
    alns=aligner.align(seq1,seq2)
    if not alns: return 0.0
    aln=alns[0]
    s1=seq1[aln.aligned[0][0][0]:aln.aligned[0][0][1]]
    s2=seq2[aln.aligned[1][0][0]:aln.aligned[1][0][1]]
    matches=sum(a==b for a,b in zip(s1,s2))
    return matches/max(len(s1),1)

def chain_sequence(structure, chain_id):
    builder=PPBuilder()
    for model in structure:
        for chain in model:
            if chain.id==chain_id:
                return "".join(str(pp.get_sequence()) for pp in builder.build_peptides(chain))
    return ""

# ────────────────────────────────────────────────
class LigaseSelect(Select):
    """
    Keep: protein atoms from keep_chains; and HET ligands explicitly in keep_ligs,
    excluding anything in EXCLUDE_ALL or that looks inorganic by structure.
    """
    def __init__(self, keep_chains, keep_ligs):
        self.keep_chains = set(keep_chains)
        self.keep_ligs = {x.upper() for x in keep_ligs}

    def accept_chain(self, chain):  # keep all chains; filter at residue level
        return True

    def accept_residue(self, res):
        chain_id = res.get_parent().id
        name = res.get_resname().strip().upper()
        het = res.id[0].strip()  # empty for polymer residues

        # Always keep the selected protein chain(s)
        if chain_id in self.keep_chains:
            return True

        # For other chains, only keep selected ligands that are not excluded and not ion-like
        if het:
            if name in self.keep_ligs and name not in EXCLUDE_ALL and not is_inorganic_like(res):
                return True

        return False

# ────────────────────────────────────────────────
def get_ligands_from_entry(pdb: str, pdb_file: Path = None) -> List[str]:
    """
    Extract ligand residue names directly from the *local PDB file*.
    Only returns unique HETATM residue names that are not in EXCLUDE_ALL.
    Captures full alphanumeric codes up to 8 characters (e.g. A1H17).
    """

    ligs = set()
    pdb = pdb.strip().upper()

    if not pdb_file or not pdb_file.exists():
        print(f"⚠️ No PDB file for {pdb} — cannot extract ligands.")
        return []

    try:
        with open(pdb_file, "r") as fh:
            for line in fh:
                if not line.startswith("HETATM"):
                    continue
                # extract residue name (columns 18–21 in standard PDB)
                code = line[17:27].strip().upper()
                # clean up extra chain / numbering info (e.g. "O6M C 9" -> "O6M")
                code = re.split(r"\s+", code)[0]

                if not code or len(code) < 2 or len(code) > 8:
                    continue
                if code in EXCLUDE_ALL:
                    continue

                ligs.add(code)

    except Exception as e:
        print(f"⚠️ Failed to read ligands from {pdb_file}: {e}")

    return sorted(ligs)




def ligand_info(code):
    url=f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}"
    try:
        r=http_get(url)
        if not r.ok: return {"mw":0,"formula":"","name":code}
        j=r.json().get("chem_comp",{})
        return {"mw":j.get("formula_weight",0),"formula":j.get("formula",""),"name":j.get("name",code)}
    except Exception:
        return {"mw":0,"formula":"","name":code}

def is_relevant_ligand(info):
    mw,form=info["mw"],info["formula"]
    if mw<MW_MIN: return False,"too_small"
    if mw>MW_MAX: return False,"too_large"
    if "C" not in form: return False,"inorganic"
    return True,"keep"

# ────────────────────────────────────────────────
def get_entry_list(terms:List[str])->List[str]:
    ids=set()
    for t in terms:
        payload={"query":{"type":"terminal","service":"full_text","parameters":{"value":t}},
                 "return_type":"entry",
                 "request_options":{"results_content_type":["experimental"],
                                    "paginate":{"start":0,"rows":200}}}
        r=http_post("https://search.rcsb.org/rcsbsearch/v2/query",json=payload)
        if not r.ok: continue
        try: j=r.json()
        except Exception: continue
        ids|={rec["identifier"] for rec in j.get("result_set",[])}
    return sorted(ids)

# ────────────────────────────────────────────────
def fetch_template_sequence(pdb_id,chain_id,base:Path):
    cache=base/f"{pdb_id}_{chain_id}_seq.txt"
    if cache.exists(): return cache.read_text().strip()
    url=f"https://files.rcsb.org/download/{pdb_id}.cif"
    r=http_get(url); r.raise_for_status()
    tmp=base/f"{pdb_id}.cif"; tmp.write_bytes(r.content)
    st=MMCIFParser(QUIET=True).get_structure(pdb_id,str(tmp))
    seq=chain_sequence(st,chain_id)
    cache.write_text(seq)
    return seq

def download_cif_structure(pdb: str, cif_dir: Path):
    """
    1. Download CIF (if not cached)
    2. Save as .tmp first
    3. Convert CIF → PDB via gemmi
    4. Copy cleaned PDB (ATOM + HETATM only) to /PDB/
    5. Return (Structure, pdb_path)
    """
    cif_dir.mkdir(parents=True, exist_ok=True)
    cif_path  = cif_dir / f"{pdb}.cif"
    pdb_tmp   = cif_dir / f"{pdb}.tmp.pdb"
    pdb_final = cif_dir / f"{pdb}.pdb"

    # ── Step 1: Download CIF if not cached ─────────────────────────────
    if not cif_path.exists():
        print(f"⬇️  Downloading CIF for {pdb}")
        url = f"https://files.rcsb.org/download/{pdb}.cif"
        r = http_get(url, timeout=45)
        r.raise_for_status()
        tmp_path = cif_dir / f"{pdb}.tmp"
        tmp_path.write_bytes(r.content)

        # Validate before renaming
        try:
            _ = gemmi.cif.read_file(str(tmp_path))
            shutil.move(str(tmp_path), str(cif_path))
        except Exception as e:
            tmp_path.unlink(missing_ok=True)
            raise RuntimeError(f"❌ Invalid CIF for {pdb}: {e}")

    # ── Step 2: Convert CIF → PDB ──────────────────────────────────────
    try:
        doc = gemmi.cif.read_file(str(cif_path))
        structure = gemmi.make_structure_from_block(doc.sole_block())
        structure.remove_hydrogens()
        structure.write_minimal_pdb(str(pdb_tmp))
        shutil.move(str(pdb_tmp), str(pdb_final))
        print(f"💾 Converted CIF → {pdb_final}")
    except Exception as e:
        print(f"⚠️ CIF→PDB conversion failed for {pdb}: {e}")
        pdb_final = None

    # ── Step 3: Copy a cleaned PDB version to /PDB/ ────────────────────
    if pdb_final and os.path.exists(pdb_final):
        pdb_clean = cif_dir.parent / "PDB" / f"{pdb}.pdb"
        os.makedirs(pdb_clean.parent, exist_ok=True)
        with open(pdb_final) as src, open(pdb_clean, "w") as dst:
            for line in src:
                if line.startswith(("ATOM", "HETATM")):
                    dst.write(line)
        print(f"💾 Copied clean PDB → {pdb_clean}")

    # ── Step 4: Parse with Biopython ───────────────────────────────────
    st = MMCIFParser(QUIET=True).get_structure(pdb, str(cif_path))
    return st, pdb_final

# ────────────────────────────────────────────────
def process_ligase(name, info):
    """
    Process a single ligase:
      - Fetch template sequence
      - Search RCSB entries
      - Download CIF and convert to PDB
      - Identify best chain
      - Detect ligands (robustly)
      - Save bound ligands and clean PDBs
    """
    base = Path(name)
    for subdir in ("CIF", "PDB", "SDF"):
        (base / subdir).mkdir(parents=True, exist_ok=True)

    log_path = base / f"{name}_log.csv"
    with open(log_path, "w", newline="") as log:
        writer = csv.DictWriter(log, fieldnames=["pdb", "ligand", "mw", "decision", "reason"])
        writer.writeheader()

        seq = fetch_template_sequence(info["template_pdb"], info["template_chain"], base)
        pdbs = get_entry_list(info["search_terms"])
        print(f"🔍 Found {len(pdbs)} structures for {name}")

        for pdb in pdbs:
            try:
                print(f"\n🔹 {name}: {pdb}")
                st, pdb_file = download_cif_structure(pdb, base / "CIF")

                # ───────────────────────────────────────────────
                # 🧩 Identify best chain
                # ───────────────────────────────────────────────
                scores = []
                for model in st:
                    for ch in model:
                        seq2 = "".join(str(pp.get_sequence()) for pp in PPBuilder().build_peptides(ch))
                        if not seq2:
                            continue
                        identity = seq_identity_local(seq, seq2)
                        scores.append((ch.id, identity))
                if not scores:
                    print("   ⚠️ No chains found with sequence.")
                    continue
                scores.sort(key=lambda x: x[1], reverse=True)
                best_id, best_sc = scores[0]
                print(f"   🧩 Best chain: {best_id} ({best_sc:.2f})")
                print("   📊 All chains:", ", ".join([f"{cid}:{sc:.2f}" for cid, sc in scores]))

                if best_sc < SEQUENCE_IDENTITY_THRESHOLD:
                    print(f"   ⚠️ Below threshold ({SEQUENCE_IDENTITY_THRESHOLD}) → skipping.")
                    writer.writerow({
                        "pdb": pdb, "ligand": "", "mw": "", 
                        "decision": "skip", "reason": "low_identity"
                    })
                    continue

                # ───────────────────────────────────────────────
                # 💊 Ligand detection and classification
                # ───────────────────────────────────────────────
                ligs = get_ligands_from_entry(pdb, pdb_file)
                print(f"💾 Converted CIF → {base}/CIF/{pdb}.pdb")
                print(f"💾 Copied clean PDB → {base}/PDB/{pdb}.pdb")
                print(f" 💊 Ligands found in FILE: {ligs}")

                if not ligs:
                    writer.writerow({
                        "pdb": pdb, "ligand": "", "mw": "", 
                        "decision": "skip", "reason": "no_ligands"
                    })
                    continue

                all_atoms = list(st.get_atoms())
                ns = NeighborSearch(all_atoms)
                main_atoms = [a for a in all_atoms if a.get_parent().get_parent().id == best_id]

                lig_in_best = set()
                lig_near_best = set()

                # 1️⃣ Primary Biopython scan
                for model in st:
                    for chain in model:
                        for res in chain:
                            rn = res.get_resname().strip().upper()
                            if rn not in ligs:
                                continue
                            if chain.id == best_id:
                                lig_in_best.add(rn)
                            else:
                                for a in main_atoms:
                                    for b in ns.search(a.coord, DISTANCE_CUTOFF):
                                        if b.get_parent().get_resname().strip().upper() == rn:
                                            lig_near_best.add(rn)
                                            break
                                    if rn in lig_near_best:
                                        break

                # 2️⃣ Fallback: parse PDB text to recover missing chain IDs
                with open(pdb_file) as fh:
                    for line in fh:
                        if line.startswith("HETATM"):
                            resname = line[17:20].strip().upper()
                            chain_id = line[21].strip()
                            if resname in ligs:
                                if chain_id == best_id:
                                    lig_in_best.add(resname)
                                elif chain_id and chain_id != best_id:
                                    lig_near_best.add(resname)

                print(f" 🔍 Ligands found in Best Chain: {sorted(lig_in_best) if lig_in_best else '[]'}")
                print(f" 💊 Ligands found Near (Bound to Best Chain): {sorted(lig_near_best) if lig_near_best else '[]'}")

                # ───────────────────────────────────────────────
                # 🧮 Ligand evaluation + SDF retrieval
                # ───────────────────────────────────────────────
                kept_ligands = []
                for l in ligs:
                    info_l = ligand_info(l)
                    keep, reason = is_relevant_ligand(info_l)
                    writer.writerow({
                        "pdb": pdb, "ligand": l, "mw": info_l["mw"],
                        "decision": "keep" if keep else "skip", "reason": reason
                    })
                    if keep:
                        kept_ligands.append(l)
                        print(f"   ✅ Kept ligand {l} ({info_l['mw']:.1f} Da)")
                        sdf_path = base / "SDF" / f"{name}_{l}.sdf"
                        for u in [
                            f"https://files.rcsb.org/ligands/download/{l}_ideal.sdf",
                            f"https://files.rcsb.org/ligands/download/{l}_model.sdf"
                        ]:
                            rr = http_get(u)
                            if rr.ok:
                                sdf_path.write_bytes(rr.content)
                                print(f"💾 Saved ligand SDF → {sdf_path}")
                                break
                        else:
                            print(f"⚠️ Failed to fetch SDF for {l}")

                if not kept_ligands:
                    continue

                # ───────────────────────────────────────────────
                # 💾 Save one filtered PDB per bound ligand
                #   - keep ONLY the best protein chain (best_id)
                #   - keep ONLY the target ligand HETATM (resname == lig_name)
                #   - drop all other HETATMs (waters/ions/etc.)
                # ───────────────────────────────────────────────
                # ───────────────────────────────────────────────
                # 💾 Save filtered PDBs for all bound ligands (in or near best chain)
                # ───────────────────────────────────────────────
                relevant_ligs = sorted(set(kept_ligands) & (lig_in_best | lig_near_best))
                if not relevant_ligs:
                    print(f"   ⚠️ No ligands to save for {pdb} after filtering.")
                else:
                    pdb_in = Path(pdb_file)
                    if not pdb_in.exists():
                        print(f"⚠️ Missing PDB file for {pdb}, skipping filtered saves.")
                    else:
                        for lig_name in relevant_ligs:
                            if lig_name in EXCLUDE_ALL:
                                print(f"   ⚠️ Skipping excluded ligand {lig_name} for {pdb}")
                                continue

                            out_pdb = base / "PDB" / f"{pdb}_{lig_name}.pdb"
                            prot_lines = lig_lines = 0

                            with open(pdb_in) as infile, open(out_pdb, "w") as out:
                                for line in infile:
                                    if not line.startswith(("ATOM", "HETATM")):
                                        continue
                                    resname = line[17:20].strip().upper()
                                    chain   = line[21].strip()

                                    # keep protein atoms from best_id, plus HETATM for this ligand
                                    if (line.startswith("ATOM") and chain == best_id) or (
                                        line.startswith("HETATM") and resname == lig_name
                                    ):
                                        out.write(line)
                                        if line.startswith("ATOM"):
                                            prot_lines += 1
                                        else:
                                            lig_lines += 1

                                out.write("END\n")

                            print(f"💾 Saved filtered PDB → {out_pdb} "
                                f"(protein lines: {prot_lines}, ligand {lig_name} lines: {lig_lines})")



            except Exception as e:
                print(f"   ⚠️ {pdb}: {e}")
                writer.writerow({
                    "pdb": pdb, "ligand": "", "mw": "", "decision": "error", "reason": str(e)
                })

    # ───────────────────────────────────────────────
    # 📦 Zip all results (once per ligase)
    # ───────────────────────────────────────────────
    zip_path = base / f"{name}_{DATE_TAG}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for subdir in ("CIF", "PDB", "SDF"):
            for f in (base / subdir).glob("*"):
                zf.write(f, f.relative_to(base))
    print(f"✅ {name} done → {zip_path.name}")


def is_inorganic_like(residue):
    """
    Heuristic: kick out classic ion-looking residues even if their names
    slipped past EXCLUDE_ALL. We treat as inorganic if:
      - <= 2 atoms total, and
      - every atom's element is in a typical ion set, OR
      - the residue's atoms contain no carbon at all.
    """
    atoms = list(residue.get_atoms())
    if not atoms:
        return True
    ion_elems = {
        "H","D","NA","K","LI","RB","CS",
        "CA","MG","SR","BA","MN","FE","CO","NI","CU","ZN","CD","HG",
        "AL","GA","IN","TL","SN","PB","AG","AU","PT",
        "CL","BR","I","F"
    }
    elems = {a.element.strip().upper() for a in atoms if getattr(a, "element", "")}
    no_carbon = "C" not in elems
    tiny = len(atoms) <= 2
    ionlike = elems.issubset(ion_elems)
    return (tiny and ionlike) or (tiny and no_carbon) or (no_carbon and ionlike)

# ────────────────────────────────────────────────
def main():
    print(f"📅 LigaseBatchFetch started — {len(LIGASES)} ligases to process")
    for name,info in LIGASES.items():
        try: process_ligase(name,info)
        except Exception as e: print(f"❌ Error processing {name}: {e}")
    print("\n🎉 All ligases processed successfully.")

# ────────────────────────────────────────────────
if __name__=="__main__":
    main()