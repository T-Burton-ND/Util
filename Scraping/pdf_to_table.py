#!/usr/bin/env python3
"""
pdf_to_table.py — Extract numeric-ish tables from PDFs into CSV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Extracts tables from every PDF in a directory and saves (a) per‑PDF CSVs
and (b) a single combined CSV of useful tables (numeric‑leaning heuristic).

Usage
-----
python pdf_to_table.py -d path/to/pdf_dir

Examples
--------
python pdf_to_table.py -d ./papers

Dependencies
------------
camelot-py, ghostscript, opencv, tk, pdfplumber (optional), pandas, tqdm
"""

import argparse
import os
import re
from pathlib import Path
from typing import List

import pandas as pd
from tqdm import tqdm

__version__ = "0.1.0"

DIGIT_TOKEN = re.compile(r"\d")

# ── Filter thresholds ──────────────────────────────────────────────
NUMERIC_RATIO_MIN = 0.20
NUMERIC_CELLS_MIN = 9


# ── Heuristic junk filter ──────────────────────────────────────────
def table_is_useful(df: pd.DataFrame) -> bool:
    if df.shape[1] < 2:
        return False
    data_rows = df.iloc[1:] if len(df) > 1 else df
    total, numericish = 0, 0
    for cell in data_rows.values.flatten():
        s = str(cell)
        if len(s) > 20 or s.startswith("("):
            continue
        total += 1
        if DIGIT_TOKEN.search(s):
            numericish += 1
    if total == 0:
        return False
    ratio = numericish / total
    return ratio >= NUMERIC_RATIO_MIN and numericish >= NUMERIC_CELLS_MIN


# ── Table extractor for a single PDF ───────────────────────────────
def extract_tables_from_pdf(pdf_path: Path) -> List[pd.DataFrame]:
    import camelot
    tables = []

    for flavor in ("lattice", "stream"):
        try:
            result = camelot.read_pdf(str(pdf_path), flavor=flavor, pages="all")
            for tbl in result:
                df = tbl.df
                if table_is_useful(df):
                    tables.append(df)
        except Exception:
            continue

    # Fallback to pdfplumber (image-based)
    if not tables:
        try:
            import pdfplumber
            with pdfplumber.open(str(pdf_path)) as pdf:
                for page in pdf.pages:
                    for raw in page.extract_tables():
                        df = pd.DataFrame(raw)
                        if table_is_useful(df):
                            tables.append(df)
        except Exception:
            pass

    return tables


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pdf_to_table.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Extract numeric-ish tables from PDFs into per-PDF CSVs and a combined CSV.",
        epilog=(
            "Usage:\n"
            "  python pdf_to_table.py -d path/to/pdf_dir\n\n"
            "Examples:\n"
            "  python pdf_to_table.py -d ./papers\n"
        ),
    )
    p.add_argument(
        "-d", "--dir", type=str, required=True,
        help="Directory containing PDFs to process"
    )
    p.add_argument(
        "--version", action="version",
        version="%(prog)s " + __version__
    )
    return p


# ── Main ───────────────────────────────────────────────────────────
def main():
    parser = build_parser()
    args = parser.parse_args()
    pdf_dir = Path(args.dir)

    all_tables = []
    output_dir = Path("per_pdf_tables")
    output_dir.mkdir(exist_ok=True)

    for pdf_file in tqdm(sorted(pdf_dir.rglob("*.pdf")), desc="Processing PDFs"):
        tables = extract_tables_from_pdf(pdf_file)
        per_pdf_rows = []

        for df in tables:
            df.insert(0, "source_pdf", pdf_file.name)
            all_tables.append(df)
            per_pdf_rows.append(df)

        if per_pdf_rows:
            rel_path = pdf_file.relative_to(pdf_dir).with_suffix("")
            rel_csv_path = output_dir / f"{rel_path}.csv"
            rel_csv_path.parent.mkdir(parents=True, exist_ok=True)
            per_pdf_df = pd.concat(per_pdf_rows, ignore_index=True)
            per_pdf_df.to_csv(rel_csv_path, index=False)

    if all_tables:
        combined = pd.concat(all_tables, ignore_index=True)
        combined.to_csv("all_tables_combined.csv", index=False)
        print(f"✅ Combined CSV saved as all_tables_combined.csv with {len(combined)} rows.")
    else:
        print("⚠️ No valid tables found.")


if __name__ == "__main__":
    main()
