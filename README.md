# 🧬 E3 Recruiter Ligandalyzer Scripts

<p align="center">
  <strong>Data-processing workflows used to build the curated datasets behind E3 Recruiter Ligandalyzer</strong>
</p>

<p align="center">
  <em>A structure-first cheminformatics and structural-biology pipeline for collecting, cleaning, annotating, standardizing, and packaging E3 ligase recruiter ligand data.</em>
</p>

<p align="center">
  <a href="https://E3Ligandalyzer.com">
    <img src="https://img.shields.io/badge/Launch-E3Ligandalyzer.com-success?style=for-the-badge&logo=googlechrome" alt="Launch E3 Recruiter Ligandalyzer">
  </a>
  <a href="https://github.com/Joey305/E3-Ligandalyzer-Scripts">
    <img src="https://img.shields.io/badge/View%20Code-GitHub-black?style=for-the-badge&logo=github" alt="View code on GitHub">
  </a>
  <a href="mailto:jmschulz@med.miami.edu?subject=E3%20Recruiter%20Ligandalyzer%20Question%20%2F%20Collaboration">
    <img src="https://img.shields.io/badge/Contact%20the%20Author-Joseph--Michael%20Schulz-blue?style=for-the-badge&logo=gmail" alt="Contact Joseph-Michael Schulz">
  </a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/E3%20Ligases-21-orange?style=for-the-badge" alt="21 curated E3 ligases">
  <img src="https://img.shields.io/badge/Unique%20Recruiter%20Ligands-494-blueviolet?style=for-the-badge" alt="494 unique recruiter ligands">
  <img src="https://img.shields.io/badge/LR%20Entries-602-lightgrey?style=for-the-badge" alt="602 LR-tagged recruiter-conformation entries">
  <img src="https://img.shields.io/badge/PDB%20Structures-568-9cf?style=for-the-badge" alt="568 recruiter-bound PDB structures">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Scaffolds-386-success?style=for-the-badge" alt="386 unique scaffolds">
  <img src="https://img.shields.io/badge/Scaffold%20Superclusters-419-informational?style=for-the-badge" alt="419 scaffold superclusters">
  <img src="https://img.shields.io/badge/RDKit-Cheminformatics-blue?style=for-the-badge" alt="RDKit cheminformatics">
  <img src="https://img.shields.io/badge/SASA-Structure%20Analysis-orange?style=for-the-badge" alt="SASA structural analysis">
</p>

---

## 🚀 Overview

**E3 Recruiter Ligandalyzer Scripts** contains the Python workflows used to build the curated datasets behind **E3 Recruiter Ligandalyzer**, a structure-first web platform for exploring E3 ligase recruiter ligands, scaffold relationships, physicochemical properties, solvent-accessible surface area (SASA), and downstream PROTAC Builder integration.

The public web platform is freely accessible at:

<p align="center">
  <a href="https://E3Ligandalyzer.com">
    <img src="https://img.shields.io/badge/Open%20the%20Platform-E3Ligandalyzer.com-success?style=for-the-badge&logo=googlechrome" alt="Open E3Ligandalyzer.com">
  </a>
</p>

This repository documents and preserves the workflows used to collect, clean, standardize, annotate, summarize, and package the recruiter ligand data used by the platform and associated manuscript.

---

## 🧭 Repository Navigation

<p align="center">
  <a href="#-quick-start">
    <img src="https://img.shields.io/badge/Get%20Started-Quick%20Start-orange?style=for-the-badge&logo=gnubash" alt="Quick start">
  </a>
  <a href="#-recommended-pipeline-order">
    <img src="https://img.shields.io/badge/View%20Pipeline-18%20Stages-blueviolet?style=for-the-badge&logo=python" alt="View pipeline order">
  </a>
  <a href="#-final-dataset-summary">
    <img src="https://img.shields.io/badge/View%20Dataset%20Summary-Final%20Counts-success?style=for-the-badge" alt="View dataset summary">
  </a>
  <a href="#-citation">
    <img src="https://img.shields.io/badge/Cite-E3%20Recruiter%20Ligandalyzer-lightgrey?style=for-the-badge&logo=googlescholar" alt="Citation">
  </a>
</p>

- [Project Summary](#-project-summary)
- [Repository Contents](#-repository-contents)
- [Expected Directory Structure](#-expected-directory-structure)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Recommended Pipeline Order](#-recommended-pipeline-order)
- [Key Output Tables](#-key-output-tables)
- [Important Data Definitions](#-important-data-definitions)
- [Final Dataset Summary](#-final-dataset-summary)
- [Database Build](#-database-build)
- [Citation](#-citation)
- [Data and Code Availability](#-data-and-code-availability)
- [Notes for Reproducibility](#-notes-for-reproducibility)
- [License](#-license)
- [Contact](#-contact)

---

## 🔬 Project Summary

**E3 Recruiter Ligandalyzer** was developed to support ligand-centric comparison of experimentally resolved E3 ligase recruiter chemotypes across curated E3 ligases.

The workflow processes recruiter-bound E3 ligase structures, extracts and standardizes ligand data, maps repeated structural instances to recruiter IDs, calculates ligand descriptors, assigns scaffold classes, computes SASA features, and builds final tables for downstream database integration.

The finalized manuscript-facing dataset contains:

| Dataset component | Count |
|---|---:|
| Curated E3 ligases | **21** |
| Unique recruiter ligand/component names | **494** |
| LR-tagged recruiter-conformation entries | **602** |
| Unique recruiter-bound PDB structures | **568** |
| Unique scaffolds | **386** |
| Scaffold superclusters | **419** |

In this workflow, a **unique recruiter ligand** refers to the original ligand/component identifier, while an **LR recruiter code** refers to a curated recruiter/conformation/instance tag. This distinction is important because a single ligand may appear in multiple crystal structures, binding poses, or conformational contexts.

---

## 📦 Repository Contents

The scripts are numbered approximately in the order they were used during dataset construction.

| Script | Purpose |
|---|---|
| `1_LIGASEfetch.py` | Searches/fetches E3 ligase structures and candidate ligand-containing complexes from structural sources such as the PDB. |
| `2_CleanPDBs.py` | Performs initial cleanup of ligase structure directories and PDB files. |
| `3_CLEANPDB_splitSINGLEhetatm.py` | Splits PDB files by single HETATM ligand entries to isolate recruiter-containing complexes. |
| `4_CLEAN_appendChainA.py` | Standardizes chain naming by appending or forcing chain identifiers where needed. |
| `5_SASA_LIGASES_AUDITED.py` | Computes audited SASA values for recruiter ligands in E3 ligase structures. |
| `5A_sasaclean.py` | Cleans and organizes SASA-related outputs and retired/mismatched structure files. |
| `6_CleanSASA.py` | Additional SASA table cleanup and normalization. |
| `7B_RenameDuplicateLigandCodes.py` | Handles duplicated ligand identifiers and assigns consistent recruiter/instance naming. |
| `8_Ligand_Metadata_From_SDF_HEAVYSAFE.py` | Extracts ligand metadata from SDF files with safeguards for heavy-atom handling. |
| `8B_Metadataclean.py` | Cleans and standardizes ligand metadata tables. |
| `9_LigaseTABLEclean.py` | Cleans ligase-level tables and prepares them for scaffold/descriptor integration. |
| `10_MCS_Matching.py` | Uses RDKit maximum common substructure matching to map ligand atoms and standardize ligand identities. |
| `10B_MapRecruiterSMILES.py` | Groups recruiter codes by SMILES and produces compact and wide recruiter–SMILES mapping tables. |
| `10C_CLean.py` | Rebuilds and synchronizes LR-code instance tables, metadata, SASA tables, and SMILES crosswalks. |
| `12_LigaseChemicalDescriptors.py` | Computes RDKit physicochemical and drug-likeness descriptors for recruiter ligands. |
| `13_TABLEclean.py` | Removes selected ligand/PDB entries from tables and related structure files when manual curation requires exclusion. |
| `14_ligasescaffolds.py` | Computes Bemis–Murcko scaffolds and initial scaffold diversity metrics. |
| `15_ADV_ScaffoldData.py` | Builds advanced scaffold-level summary tables. |
| `15B_ADV_ScafCluster.py` | Assigns scaffold superclusters and generates supercluster summary, frequency, and matrix tables. |
| `15C_ScaffoldUnifier.py` | Unifies scaffold identifiers across intermediate scaffold tables. |
| `16_TBLCreation.py` | Builds final crosswalk and master mapping tables, including recruiter metadata and duplicate-ligand summaries. |
| `16B_TableClean.py` | Performs final table cleaning and harmonization before database construction. |
| `17__BUILDdb.py` | Builds the SQLite database from finalized CSV tables in `Ligase_Table/`. |
| `18_Data_Summary.py` | Generates a manuscript-ready text summary of final dataset counts, descriptor summaries, scaffold counts, and SASA statistics. |

---

## 📁 Expected Directory Structure

The scripts assume a working project directory similar to:

```text
E3Ligandalyzer/
├── Ligases/
│   ├── CRBN/
│   ├── VHL/
│   ├── MDM2/
│   └── ...
├── Ligase_Table/
│   ├── Recruiter_Master_Map.csv
│   ├── Ligand_Instance_Recruiter_Codes.csv
│   ├── Ligase_Ligand_Metadata.csv
│   ├── Ligase_Ligand_SASA_summary.csv
│   ├── Ligase_Ligand_SASA_atoms.csv
│   ├── Ligase_Chemical_Descriptors.csv
│   ├── Ligase_Recruiters_Scaffold.csv
│   └── ...
├── *.py
└── README.md
```

The exact file set may vary depending on which stage of the pipeline is being reproduced.

---

## ⚙️ Installation

These scripts were developed for a Python/RDKit cheminformatics environment.

```bash
conda create -n e3ligandalyzer python=3.10 -y
conda activate e3ligandalyzer
conda install -c conda-forge rdkit pandas numpy biopython requests gemmi -y
pip install tqdm
```

Some scripts may use only standard-library modules. If additional imports are required, install them into the same environment.

---

## 🚀 Quick Start

Clone the repository:

```bash
git clone https://github.com/Joey305/E3-Ligandalyzer-Scripts.git
cd E3-Ligandalyzer-Scripts
```

Create and activate the environment:

```bash
conda create -n e3ligandalyzer python=3.10 -y
conda activate e3ligandalyzer
conda install -c conda-forge rdkit pandas numpy biopython requests gemmi -y
```

Run the final dataset summary after placing the finalized `Ligase_Table/` and `Ligases/` folders in the project root:

```bash
python 18_Data_Summary.py
```

This generates:

```text
E3Recruiter_Data_Summary.txt
Ligase_Table/Derived_Ligand_Conformation_Summary.csv
Ligase_Table/Derived_RecruiterCode_To_Ligand_Map.csv
```

---

## 🧬 Recommended Pipeline Order

The full historical workflow was iterative and included manual curation. The approximate processing order is:

```bash
python 1_LIGASEfetch.py
python 2_CleanPDBs.py
python 3_CLEANPDB_splitSINGLEhetatm.py
python 4_CLEAN_appendChainA.py
python 5_SASA_LIGASES_AUDITED.py
python 5A_sasaclean.py
python 6_CleanSASA.py
python 7B_RenameDuplicateLigandCodes.py
python 8_Ligand_Metadata_From_SDF_HEAVYSAFE.py
python 8B_Metadataclean.py
python 9_LigaseTABLEclean.py
python 10_MCS_Matching.py
python 10B_MapRecruiterSMILES.py
python 10C_CLean.py
python 12_LigaseChemicalDescriptors.py
python 13_TABLEclean.py
python 14_ligasescaffolds.py
python 15_ADV_ScaffoldData.py
python 15B_ADV_ScafCluster.py
python 15C_ScaffoldUnifier.py
python 16_TBLCreation.py
python 16B_TableClean.py
python 17__BUILDdb.py
python 18_Data_Summary.py
```

Some scripts are designed to be rerun after manual curation or table updates. Inspect generated logs, backup files, and output tables before proceeding to the next stage.

---

## 📊 Key Output Tables

The finalized workflow produces or expects several important tables in `Ligase_Table/`.

| Table | Description |
|---|---|
| `Recruiter_Master_Map.csv` | Master mapping between LR recruiter codes and original ligand/component names. |
| `Recruiter_SMILES_Map.csv` | Mapping between LR recruiter codes and SMILES strings. |
| `Recruiter_SMILES_Wide.csv` | Wide-format grouping of recruiter codes by SMILES. |
| `Ligand_Instance_Recruiter_Codes.csv` | Curated LR-tagged recruiter/conformation entries with ligase/PDB context. |
| `Recruiter_Code_Crosswalk.csv` | Crosswalk linking recruiter codes, ligand names, SMILES, structures, and metadata. |
| `Ligase_Ligand_Metadata.csv` | Ligand metadata table after cleanup and LR-code integration. |
| `Ligase_Ligand_SASA_summary.csv` | Ligand-level SASA summary table. |
| `Ligase_Ligand_SASA_atoms.csv` | Atom-level SASA table. |
| `Ligase_Chemical_Descriptors.csv` | RDKit descriptor table for recruiter entries. |
| `Ligase_Recruiters_Scaffold.csv` | Mapping between recruiters and scaffold IDs. |
| `Ligase_Scaffold_Data.csv` | Scaffold-level summary data. |
| `Ligase_Scaffold_Frequency.csv` | Scaffold frequency data by ligase. |
| `Ligase_Scaffold_Superclusters.csv` | Scaffold supercluster assignments. |
| `Ligase_Scaffold_Supercluster_Frequency.csv` | Supercluster frequency table. |
| `Ligase_Scaffold_Supercluster_Matrix.csv` | Matrix representation of scaffold supercluster coverage. |
| `Scaffold_Unified_Map.csv` | Unified scaffold identifier mapping. |
| `Derived_Ligand_Conformation_Summary.csv` | Derived table reporting how many LR codes/conformations map to each ligand. |
| `Derived_RecruiterCode_To_Ligand_Map.csv` | Derived expanded mapping between LR recruiter codes and ligand/component names. |

---

## 🧠 Important Data Definitions

### Unique recruiter ligand

A unique recruiter ligand is defined using the original ligand/component name from `Recruiter_Master_Map.csv`, such as a PDB chemical component code.

### LR recruiter code

An LR code, such as `LR00001`, is a curated recruiter/conformation/instance identifier. Multiple LR codes may map to the same underlying ligand when that ligand appears in multiple structures, poses, or conformational contexts.

### Recruiter-bound structure

A recruiter-bound structure refers to a unique PDB/structure identifier associated with a curated E3 ligase–recruiter complex.

### Scaffold

Scaffolds are generated using Bemis–Murcko framework decomposition and unified across scaffold tables using the scaffold unification workflow.

### Scaffold supercluster

Scaffold superclusters group related scaffold identifiers to support higher-level comparison of recruiter chemotypes across ligases.

---

## 📌 Final Dataset Summary

The final manuscript-facing summary can be regenerated with:

```bash
python 18_Data_Summary.py
```

Current finalized counts:

```text
494 unique ligand/component recruiters
602 curated LR-tagged recruiter-conformation entries
568 unique recruiter-bound PDB structures
21 curated E3 ligases
386 unique scaffolds
419 scaffold superclusters
```

These counts should be regenerated whenever upstream curation or table-generation scripts are rerun.

---

## 🗄️ Database Build

After the CSV tables in `Ligase_Table/` have been finalized, build the SQLite database with:

```bash
python 17__BUILDdb.py
```

This script imports CSV tables from `Ligase_Table/` into a SQLite database for use by the web application.

---

## 🧬 Citation

If you use this workflow, code, or generated dataset, please cite the associated E3 Recruiter Ligandalyzer manuscript when available.

Suggested citation text while the manuscript is in preparation:

```text
Schulz JM, Reynolds RC, Schürer SC. E3 Recruiter Ligandalyzer: A Structure-First, Ligand-Centric Platform for E3 Ligase Recruiter Selection. Manuscript in preparation.
```

---

## 🌐 Data and Code Availability

The E3 Recruiter Ligandalyzer web platform is freely accessible at:

<p align="center">
  <a href="https://E3Ligandalyzer.com">
    <img src="https://img.shields.io/badge/Web%20Platform-E3Ligandalyzer.com-success?style=for-the-badge&logo=googlechrome" alt="E3Ligandalyzer.com">
  </a>
</p>

The scripts used to build and summarize the dataset are available at:

<p align="center">
  <a href="https://github.com/Joey305/E3-Ligandalyzer-Scripts">
    <img src="https://img.shields.io/badge/Code-GitHub%20Repository-black?style=for-the-badge&logo=github" alt="GitHub repository">
  </a>
</p>

Processed tables and database files may be added to the repository or released separately depending on manuscript and repository size requirements.

---

## 🔁 Notes for Reproducibility

- Run scripts from the project root unless a script-specific path is provided.
- Inspect backup files and logs after each cleaning step.
- Some curation steps are intentionally manual or semi-manual to remove crystallographic additives, nonspecific binders, duplicate ligand names, and inappropriate entries.
- Treat `18_Data_Summary.py` as the final consistency check before reporting manuscript-level counts.
- If any upstream table changes, rerun the downstream scaffold, descriptor, database, and summary steps as needed.

---

## ⚖️ License

A license should be added before public release.

Recommended options:

- **MIT License** for permissive code reuse.
- **Apache 2.0 License** for permissive reuse with explicit patent language.
- **CC BY 4.0** or **CC BY-NC 4.0** for data tables, if released separately.

---

## 📬 Contact

For questions, collaboration, or dataset inquiries:

<p align="center">
  <a href="mailto:jmschulz@med.miami.edu?subject=E3%20Recruiter%20Ligandalyzer%20Question%20%2F%20Collaboration">
    <img src="https://img.shields.io/badge/Joseph--Michael%20Schulz-jmschulz%40med.miami.edu-blue?style=for-the-badge&logo=gmail" alt="Email Joseph-Michael Schulz">
  </a>
  <a href="https://github.com/Joey305">
    <img src="https://img.shields.io/badge/GitHub-Joey305-black?style=for-the-badge&logo=github" alt="Joey305 GitHub profile">
  </a>
</p>

**Joseph-Michael Schulz**  
University of Miami | Schürer Lab  
GitHub: [Joey305](https://github.com/Joey305)

---

## 🙌 Practical Takeaway

This repository preserves the reproducible data-processing backbone behind **E3 Recruiter Ligandalyzer**: from structure retrieval and ligand extraction to SASA analysis, scaffold clustering, descriptor calculation, master-table generation, SQLite database construction, and manuscript-ready dataset summaries.
