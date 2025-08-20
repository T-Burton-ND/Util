# Util Repository – README

This repository contains a collection of utility scripts organized into categories for **automation**, **scraping**, and **visualization**. These scripts were developed for research computing, data extraction, and workflow acceleration. The sections below include **auto-generated blocks** that can be refreshed by running `update_readme.py`.

---

## 📂 Repository Structure (high-level)
```
util/
├── automation/                    # Batch automation scripts
├── scraping/                      # Data extraction and processing
└── Visualization_Scripts/         # Visualization and chemistry helpers
```

---

## 📜 Script Overview

<!-- BEGIN AUTO-OVERVIEW -->
**Inventory summary** — 10 scripts (Python: 7, Bash: 3) • Other: 3 • Scraping: 3 • Visualization: 4 • Last scan: 2025-08-20 11:08

<!-- END AUTO-OVERVIEW -->

---

## 🔧 Detailed Inventory

<!-- BEGIN AUTO-SCRIPTS -->
### 🔹 Other

| Script | Lang | Description | Usage | LOC | Modified |
|---|---|---|---|---|---|
| [`bash_template.sh`](Starting_Points/bash_template.sh) | Bash | <file_name>.sh — <one-line summary> |  | 31 | 2025-08-20 |
| [`python_template.py`](Starting_Points/python_template.py) | Python | No --help detected | usage: python_template.py [-h] -i INPUT [-o OUTPUT] | 34 | 2025-08-20 |
| [`update_readme.py`](update_readme.py) | Python | No --help detected | usage: update_readme.py [-h] [--repo-root REPO_ROOT] [--readme README]
[--include-ext [INCLUDE_EXT ...]]
[--exclude-dirs [EXCLUDE_DIRS ...]]
[--no-ignore-hidden] [--write] [--write-index] | 332 | 2025-08-20 |

### 🔹 Scraping

| Script | Lang | Description | Usage | LOC | Modified |
|---|---|---|---|---|---|
| [`csv_to_excel.py`](Scraping/csv_to_excel.py) | Python | No --help detected | usage: combine_csv_to_excel.py [-h] [-o OUTPUT] [-r] [-s {tabs,single}]
[--engine {xlsxwriter,openpyxl}] [--version]
directory | 152 | 2025-08-20 |
| [`pdf_to_table.py`](Scraping/pdf_to_table.py) | Python | No --help detected | usage: pdf_to_table.py [-h] -d DIR [--version] | 121 | 2025-08-20 |
| [`semantic_scraper.py`](Scraping/semantic_scraper.py) | Python | No --help detected | usage: semantic_scraper.py [-h] [-n NUM] -q QUERY [-o OUTDIR] [--version] | 455 | 2025-08-20 |

### 🔹 Visualization

| Script | Lang | Description | Usage | LOC | Modified |
|---|---|---|---|---|---|
| [`bulk_mol_id.sh`](Visualization_Scripts/bulk_mol_id.sh) | Bash | bulk_mol_id.sh — Batch molecule rendering + metadata rename |  | 112 | 2025-08-20 |
| [`gif_from_mov.sh`](Visualization_Scripts/gif_from_mov.sh) | Bash | gif_from_mov.sh — Convert .mov to .gif with tunable settings |  | 61 | 2025-08-20 |
| [`mol_id.py`](Visualization_Scripts/mol_id.py) | Python | No --help detected | usage: mol_id.py [-h] [-o OUTPUT] [--no-multicolor] [--version] input | 347 | 2025-08-20 |
| [`quick_react.py`](Visualization_Scripts/quick_react.py) | Python | No --help detected | usage: quick_react.py [-h] [--bond-cutoff BOND_CUTOFF] [--size W H]
[--dpi DPI] [--ms MS] [--frame-digits FRAME_DIGITS]
[--outdir OUTDIR] [--gif GIF] [--version]
xyz | 161 | 2025-08-20 |
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
```
Advanced options:
```bash
python update_readme.py --run-help          # try to capture -h/--help for scripts
python update_readme.py --include-ext .py .sh
python update_readme.py --exclude-dirs .git .venv build dist
python update_readme.py --write-index       # write script_index.json + .csv
```

The updater will:
- Scan directories recursively (auto-detect new folders and files)
- Parse Python **docstrings** and Bash **header comments**
- Infer **usage** and **dependencies** (best-effort)
- Compute **lines of code**, **last modified** timestamps, and **language stats**
- Replace only the content between the `AUTO-OVERVIEW` and `AUTO-SCRIPTS` markers
- Optionally write a machine-readable `script_index.json` and `.csv`

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
