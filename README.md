# Util Repository – README

This repository contains a collection of utility scripts organized into categories for **automation**, **scraping**, and **visualization**. These scripts were developed for research computing, data extraction, and workflow acceleration. The sections below include **auto-generated blocks** that can be refreshed by running `update_readme.py`.

---
![Contributions](https://img.shields.io/badge/contributions-welcome-orange.svg) ![Python](https://img.shields.io/badge/language-Python-blue.svg?logo=python) ![Bash](https://img.shields.io/badge/language-Bash-green.svg?logo=gnu-bash)

## 📂 Repository Structure (high-level)
```
util/
├── Starting_Points/                # templates for .sh and .py files
├── Scraping/                      # Data extraction and processing
├── LAMMPS_Stand/                 # Automation Pipeline for large scale LAMMPS runs with OPLS forcefields
└── Visualization_Scripts/         # Visualization and chemistry helpers
```

---

### 📌 Directory Purposes
- ![Automation](https://img.shields.io/badge/category-automation-lightgrey)  
  Tools for repetitive workflows (e.g., batch runs on CRC cluster).  

- ![Scraping](https://img.shields.io/badge/category-scraping-blue)  
  Pipelines for **large-scale data collection**, especially from Semantic Scholar (PDF parsing, table extraction, CSV aggregation).  

- ![Visualization](https://img.shields.io/badge/category-visualization-green)  
  Helpers for **computational chemistry figure generation**, molecular rendering, and quick reaction-path visualizations.  

- ![LAMMPS MD](https://img.shields.io/badge/LAMMPS-MD-lightgrey?style=for-the-badge&logo=opensourceinitiative)
  Pipeline for **LAMMPS runs with OPLS files**, UGE/qsub automation intended, with config.yaml for specifics

👉 For machine-readable metadata, see `script_index.json` or `script_index.csv`.

---

## 📜 Script Overview

<!-- BEGIN AUTO-OVERVIEW -->
**Inventory summary** — 29 scripts (Python: 17, Bash: 12) • Other: 3 • LAMMPS: 16 • Scraping: 3 • Visualization: 7 • Last scan: 2025-09-09 12:40

<!-- END AUTO-OVERVIEW -->

---

## 🔧 Detailed Inventory

<!-- BEGIN AUTO-SCRIPTS -->
### 🔹 Other

| Script | Lang | Description | Modified |
|---|---|---|---|
| [`bash_template.sh`](Starting_Points/bash_template.sh) | Bash | Examples: | 2025-08-20 |
| [`python_template.py`](Starting_Points/python_template.py) | Python | Short, human-readable sentence about what this script does. | 2025-08-20 |
| [`update_readme.py`](update_readme.py) | Python | Auto-update README tables from repository scripts (help-only). | 2025-08-20 |


### 🔹 Scraping

| Script | Lang | Description | Modified |
|---|---|---|---|
| [`csv_to_excel.py`](Scraping/csv_to_excel.py) | Python | Combine CSV files into an Excel workbook. | 2025-08-20 |
| [`pdf_to_table.py`](Scraping/pdf_to_table.py) | Python | Extract numeric-ish tables from PDFs into per-PDF CSVs and a combined CSV. | 2025-08-20 |
| [`semantic_scraper.py`](Scraping/semantic_scraper.py) | Python | Scrape numeric-leaning tables from OA PDFs via Semantic Scholar. | 2025-08-20 |

### 🔹 LAMMPS Stand

| Script | Lang | Description | Modified |
|---|---|---|---|
| [`batch.sh`](LAMMPS_Stand/batch_automation/batch.sh) | Bash | Main batch script, customizable for replicate number and file locations | 2025-09-09 |
| ['config.yaml'] (LAMMPS_Stand/batch_automation/batch.sh) | YAML | Setpoints and preferences for LAMMPS runs, handled by batch.sh | 2025-09-09 |
| [`run_job.sh`](LAMMPS_Stand/batch_automation/gen_scripts/run_job.sh) | Bash | UGE/qsub job script for individual runs, handled by batch.sh | 2025-09-09 |
| [`gen_in.py`](LAMMPS_Stand/batch_automation/gen_scripts/gen_in.py) | Python | Source file for lammps input files. Uses config.yaml keys for parameters | 2025-09-09 |
| [`OPLS_box.py`](LAMMPS_Stand/batch_automation/gen_scripts/OPLS_box.py) | Python | OPLS forcefield -> lammps system file convertor, handled by batch.sh | 2025-09-09 |
| [`compute_diffusion.py`](LAMMPS_Stand/batch_automation/calc_scripts/compute_diffusion.py) | Python | Homebrewed diffusion calculator using the lammps-generated RDF files | 2025-09-09 |
| [`compute_diffusion_mdanalysis.py`](LAMMPS_Stand/batch_automation/calc_scripts/compute_diffusion_mdanalysis.py) | Python | 3rd party diffusion calculator using the lammps-generated trajectories | 2025-09-09 |
| [`data_comp.sh`](LAMMPS_Stand/batch_automation/data_handle/data_comp.sh) | Bash | Compile all diffusion data from the .flag and results.txt files | 2025-09-09 |
| [`diffusion_job.sh`](LAMMPS_Stand/batch_automation/data_handle/diffusion_job.sh) | Bash | Job script for redoing post-hoc diffusion analysis if needed | 2025-09-09 |
| [`submit_all_diffusion.sh`](LAMMPS_Stand/batch_automation/data_handle/submit_all_diffusion.sh) | Bash | Batching script for diffusion_job.sh that submits it recursively for completed jobs | 2025-09-09 |
| [`monitor_runs.py`](LAMMPS_Stand/batch_automation/data_handle/monitor_runs.py) | Python | UGE/qsub engineering monitor reporting run state and exit flags | 2025-09-09 |
| [`hall_monitor.sh`](LAMMPS_Stand/batch_automation/data_handle/hall_monitor.sh) | Bash | UGE/qsub job script for monitor_runs.py | 2025-09-09 |
| [`verify_opls_vs_runlist.py`](LAMMPS_Stand/batch_automation/data_handle/verify_opls_vs_runlist.py) | Python | Confirms all OPLS files needed for run_list.csv are available | 2025-09-09 |
| [`scrape_py.py`](LAMMPS_Stand/ligpargen_automation/scrape_py.py) | Python | SMILES string -> OPLS .lmp file for lammps inputs, using the ligpargen software on chrome browser | 2025-09-09 |
| [`ligpargen_batch.sh`](LAMMPS_Stand/ligpargen_automation/ligpargen_batch.sh) | Bash | batch run scrape_py.py for OPLS file generation from SMILES strings | 2025-09-09 |
| [`make_elec.py`](LAMMPS_Stand/ligpargen_automation/make_elec.py) | Python | Tailored .csv -> .xslx file convertor to prep files for scrape_py.py | 2025-09-09 |

### 🔹 Visualization

| Script | Lang | Description | Modified |
|---|---|---|---|
| [`bulk_mol_id.sh`](Visualization_Scripts/bulk_mol_id.sh) | Bash | Environment variables: | 2025-09-02 |
| [`callgraph.py`](Visualization_Scripts/callgraph.py) | Python | Build and render a cleaned, layered call-graph SVG. | 2025-08-20 |
| [`gif_from_mov.sh`](Visualization_Scripts/gif_from_mov.sh) | Bash | Options (edit variables at top of file): | 2025-08-20 |
| [`harness.py`](Visualization_Scripts/harness.py) | Python | No --help detected | 2025-08-20 |
| [`mol_id.py`](Visualization_Scripts/mol_id.py) | Python | Draw a 2D molecule image from a SMILES, InChI, or InChIKey. | 2025-09-02 |
| [`quick_react.py`](Visualization_Scripts/quick_react.py) | Python | Create a rocking (forward-then-back) GIF from a multi-frame XYZ trajectory. | 2025-08-20 |
| [`run_smiles.sh`](Visualization_Scripts/run_smiles.sh) | Bash | head: reaction_list.txt: No such file or directory | 2025-09-02 |
<!-- END AUTO-SCRIPTS -->

---

## 🚀 Quickstart

### Clone and set up environment
```bash
git clone https://github.com/T-Burton-ND/Util
cd Util
conda env create -f environment.yml
conda activate util_env
```

### Typical usage
```bash
# Convert CSV directory into a single Excel workbook
python scraping/csv_to_excel.py data/ -s single -o combined.xlsx

# Convert a .mov into a .gif
bash Visualization_Scripts/gif_from_mov.sh input.mov ./gifs

# Draw a molecule and fetch metadata
python Visualization_Scripts/mol_id.py "CCO"
```

> **CRC @ Notre Dame (UGE/qsub):** For long runs (e.g., large scrapes), submit a batch job with `qsub`. Most scripts here are lightweight and fine to run on a login node.

---

## 🛠 Dependencies
Install the tools used across scripts (conda-forge recommended):
```bash
conda install -c conda-forge rdkit requests pandas tqdm camelot-py ghostscript opencv tk pdfplumber ase matplotlib pillow xlsxwriter openpyxl
# Plus: ffmpeg (for GIFs)
conda install -c conda-forge ffmpeg
```

---

## 🔄 Auto-updating this README
Run the updater to refresh the **Overview** and **Inventory** tables:
```bash
python update_readme.py            # dry-run preview printed to console
python update_readme.py --write    # apply changes to README.md
python update_readme.py --write-index, # Write script_index.json and script_index.csv
```

The updater will:
- Scan directories recursively (auto-detect new folders and files)
- Parse Python **docstrings** and Bash **header comments**
- Infer **usage** and **dependencies** (best-effort)
- Compute **lines of code**, **last modified** timestamps, and **language stats**
- Replace only the content between the `AUTO-OVERVIEW` and `AUTO-SCRIPTS` markers
- Optionally write a machine-readable `script_index.json` and `.csv`

---

## 🤝 Contributing
Pull requests are welcome! Please ensure scripts:
- Include a short header docstring/help (`-h/--help`)
- Follow existing folder structure
- Add dependencies to `environment.yml`

---

## 🧪 Testing (optional)
If you add tests, place them in `tests/` and run:
```bash
pytest -q
```

---

## 👤 Author
Maintained by **Thomas J. Burton** – Savoie Research Group, University of Notre Dame

For questions or contributions, please open an issue or pull request.

## 📖 Citation
If you use these scripts in your research, please cite:

Thomas J. Burton, *Util Repository*, 2025.  
See the [CITATION.cff](./CITATION.cff) file for machine-readable citation info.

![License](https://img.shields.io/badge/license-MIT-green.svg)
