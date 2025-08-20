# Util Repository – README

This repository contains a collection of utility scripts organized into categories for **automation**, **scraping**, and **visualization**. These scripts were developed for research computing, data extraction, and workflow acceleration. The sections below include **auto-generated blocks** that can be refreshed by running `update_readme.py`.

---
![Python](https://img.shields.io/badge/python-3.8+-blue.svg
## 📂 Repository Structure (high-level)
```
util/
├── Starting_Points/                # templates for .sh and .py files
├── scraping/                      # Data extraction and processing
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

👉 For machine-readable metadata, see `script_index.json` or `script_index.csv`.

---

## 📜 Script Overview

<!-- BEGIN AUTO-OVERVIEW -->
**Inventory summary** — 10 scripts (Python: 7, Bash: 3) • Other: 3 • Scraping: 3 • Visualization: 4 • Last scan: 2025-08-20 11:50

<!-- END AUTO-OVERVIEW -->

---

## 🔧 Detailed Inventory

<!-- BEGIN AUTO-SCRIPTS -->
### 🔹 Other

| Script | Lang | Description | Usage | LOC | Modified |
|---|---|---|---|---|---|
| [`bash_template.sh`](Starting_Points/bash_template.sh) | Bash | <file_name>.sh — <one-line summary> |  | 31 | 2025-08-20 |
| [`python_template.py`](Starting_Points/python_template.py) | Python | Short, human-readable sentence about what this script does. | usage: python_template.py [-h] -i INPUT [-o OUTPUT] | 34 | 2025-08-20 |
| [`update_readme.py`](update_readme.py) | Python | [--include-ext [INCLUDE_EXT ...]] | usage: update_readme.py [-h] [--repo-root REPO_ROOT] [--readme README]

### 🔹 Scraping

| Script | Lang | Description | Usage | LOC | Modified |
|---|---|---|---|---|---|
| [`csv_to_excel.py`](Scraping/csv_to_excel.py) | Python | [--engine {xlsxwriter,openpyxl}] [--version] | usage: combine_csv_to_excel.py [-h] [-o OUTPUT] [-r] [-s {tabs,single}]
| [`pdf_to_table.py`](Scraping/pdf_to_table.py) | Python | Extract numeric-ish tables from PDFs into per-PDF CSVs and a combined CSV. | usage: pdf_to_table.py [-h] -d DIR [--version] | 121 | 2025-08-20 |
| [`semantic_scraper.py`](Scraping/semantic_scraper.py) | Python | Scrape numeric-leaning tables from OA PDFs via Semantic Scholar. | usage: semantic_scraper.py [-h] [-n NUM] -q QUERY [-o OUTDIR] [--version] | 455 | 2025-08-20 |

### 🔹 Visualization

| Script | Lang | Description | Usage | LOC | Modified |
|---|---|---|---|---|---|
| [`bulk_mol_id.sh`](Visualization_Scripts/bulk_mol_id.sh) | Bash | bulk_mol_id.sh — Batch molecule rendering + metadata rename |  | 112 | 2025-08-20 |
| [`gif_from_mov.sh`](Visualization_Scripts/gif_from_mov.sh) | Bash | gif_from_mov.sh — Convert .mov to .gif with tunable settings |  | 61 | 2025-08-20 |
| [`mol_id.py`](Visualization_Scripts/mol_id.py) | Python | Draw a 2D molecule image from a SMILES, InChI, or InChIKey. | usage: mol_id.py [-h] [-o OUTPUT] [--no-multicolor] [--version] input | 347 | 2025-08-20 |
| [`quick_react.py`](Visualization_Scripts/quick_react.py) | Python | [--dpi DPI] [--ms MS] [--frame-digits FRAME_DIGITS] | usage: quick_react.py [-h] [--bond-cutoff BOND_CUTOFF] [--size W H]

<!-- END AUTO-SCRIPTS -->

---

## 🚀 Quickstart

### Clone and set up environment
```bash
git clone <your-repo-url> util
cd util
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
![Contributions](https://img.shields.io/badge/contributions-welcome-orange.svg)
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
