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
**Inventory summary** — 11 scripts (Python: 7, Bash: 4) • Automation: 1 • Other: 3 • Scraping: 3 • Visualization: 4 • Last scan: 2025-08-19 16:05

<!-- END AUTO-OVERVIEW -->

---

## 🔧 Detailed Inventory

<!-- BEGIN AUTO-SCRIPTS -->
### 🔹 Automation

| Script | Lang | Description | Usage | Deps | LOC | Modified |
|---|---|---|---|---|---|---|
| [`batch_run_mol_id.sh`](Automation/batch_run_mol_id.sh) | Bash | Bash script | Usage  : ./batch_run_mol_id.sh
------------------------------------------------------------------------------ | bash, grep, head, python, sed, tr | 93 | 2025-08-19 |

### 🔹 Other

| Script | Lang | Description | Usage | Deps | LOC | Modified |
|---|---|---|---|---|---|---|
| [`bash_template.sh`](Templates/bash_template.sh) | Bash | <file_name>.sh — <one‑line summary> | Usage:
./<file_name>.sh [OPTIONS] <args> | awk, bash, ffmpeg, sed, sh | 35 | 2025-08-19 |
| [`python_template.py`](Templates/python_template.py) | Python | <file_name>.py — <one‑line summary> ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Short description (2–3 sentences max). | usage: python_template.py [-h] -i INPUT [-o OUTPUT] |  | 34 | 2025-08-19 |
| [`update_readme.py`](update_readme.py) | Python | update_readme.py ─────────────────────────────────────────────────────────────────────── Auto-refreshes README sections by scanning the repository for scripts. Features • Recursively discovers .py/.sh scripts, skipping junk (e.g., .git, venvs) • Extracts: – Title (filename), relative path – Language (Python/Bash) – Description (module docstring or header comment) – Usage (from docstring, header, or optional --help capture) – Guessed dependencies (imports / external commands) – Lines of code, las… | usage: update_readme.py [-h] [--repo-root REPO_ROOT] [--readme README]
[--include-ext [INCLUDE_EXT ...]]
[--exclude-dirs [EXCLUDE_DIRS ...]] |  | 391 | 2025-08-19 |

### 🔹 Scraping

| Script | Lang | Description | Usage | Deps | LOC | Modified |
|---|---|---|---|---|---|---|
| [`csv_to_excel.py`](Scraping/csv_to_excel.py) | Python | combine_csv_to_excel.py ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Author : Thomas J. Burton – Savoie Research Group, UND Updated: 2025-06-06 License: MIT Description | usage: csv_to_excel.py [-h] [-o OUTPUT] [-r] [-s {tabs,single}]
[--engine {xlsxwriter,openpyxl}]
directory | pandas | 141 | 2025-06-19 |
| [`pdf_to_table.py`](Scraping/pdf_to_table.py) | Python | -extract_all_tables.py - -Extracts tables from every PDF in a directory and saves all useful tables -into a single combined CSV file. - -Usage: -    python pdf_to_table.py -d path/to/pdf_dir | usage: pdf_to_table.py [-h] -d DIR | camelot, pandas, pdfplumber, tqdm | 113 | 2025-08-19 |
| [`semantic_scraper.py`](Scraping/semantic_scraper.py) | Python | semantic_scraper.py ──────────────────────────────────────────────────────────────────────── A beginner-friendly pipeline that 1. queries the **Semantic Scholar Graph API** for open-access (OA) papers, 2. downloads their PDFs (robustly and resumably), 3. extracts every table the PDF parser can detect, 4. keeps only tables that look numeric (→ e.g. diesel / elemental analysis), 5. writes a run summary + full log so you can audit failures. The script prints **only tqdm progress bars** to the termi… | Usage:
python semantic_scraper.py -q 'Keyword Search" -n 999 -o output_dir -s tabs | camelot, pandas, pdfplumber, requests, tqdm | 438 | 2025-06-18 |

### 🔹 Visualization

| Script | Lang | Description | Usage | Deps | LOC | Modified |
|---|---|---|---|---|---|---|
| [`bulk_mol_id.sh`](Visualization_Scripts/bulk_mol_id.sh) | Bash | Bash script |  | awk, bash, mktemp, python, python3, sed, tee | 97 | 2025-08-19 |
| [`gif_from_mov.sh`](Visualization_Scripts/gif_from_mov.sh) | Bash | gif_from_mov.sh — Convert .mov to .gif with tunable settings | Usage: /Users/tburton2/Desktop/Repos/Util/Visualization_Scripts/gif_from_mov.sh path/to/input.mov path/to/output_folder
Example: ./gif_from_mov.sh input.mov ./gifs | bash, ffmpeg | 48 | 2025-08-19 |
| [`mol_id.py`](Visualization_Scripts/mol_id.py) | Python | mol_id.py – Molecule visualizer with online InChIKey support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Author  : Thomas J. Burton – Savoie Research Group, UND Updated : 2025-07-22 License : MIT Description | usage: mol_id.py [-h] [-o OUTPUT] [--no-multicolor] input | requests | 319 | 2025-08-07 |
| [`quick_react.py`](Visualization_Scripts/quick_react.py) | Python | -quick_react.py -Usage: python quick_react.py your_trajectory.xyz -Description: Loads a multi-frame XYZ trajectory, infers bonds manually using covalent radii and a bond cutoff, and creates a rocking (forward-reverse) GIF animation using matplotlib. Works with older versions of ASE. | Usage:
python quick_react.py your_trajectory.xyz | matplotlib | 129 | 2025-08-19 |
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
