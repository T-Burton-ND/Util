#!/usr/bin/env python3
"""
combine_csv_to_excel.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Author : Thomas J. Burton – Savoie Research Group, UND
Updated: 2025-06-06
License: MIT

Description
-----------
Walk through a directory, find every ``*.csv`` file, and write them into an
Excel workbook. Two layout styles are supported:

1. ``tabs``   – each CSV becomes its own worksheet (default).
2. ``single`` – all CSVs are written sequentially in **one** worksheet,
                separated by two blank rows.

Usage
-----
python combine_csv_to_excel.py <directory> [-r] [-s {tabs,single}] [-o OUTPUT]

Examples
--------
# One worksheet per CSV (default)
python combine_csv_to_excel.py data/

# All tables in a single worksheet, separated by blank lines
python combine_csv_to_excel.py data/ -s single -o all_tables.xlsx

Dependencies
------------
conda install -c conda-forge pandas xlsxwriter   # (openpyxl optional)
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import textwrap

import pandas as pd


# --------------------------------------------------------------------------- #
# Helper functions
# --------------------------------------------------------------------------- #
def find_csv_files(directory: pathlib.Path, recursive: bool) -> list[pathlib.Path]:
    pattern = "**/*.csv" if recursive else "*.csv"
    return sorted(directory.glob(pattern))


def sanitize_sheet_name(name: str) -> str:
    """Excel sheet names ≤31 chars, without []:*?/\\."""
    for ch in "[]:*?/\\":  # forbidden chars
        name = name.replace(ch, "_")
    return name[:31] or "sheet"


def csvs_to_excel(
    csv_paths: list[pathlib.Path],
    output_path: pathlib.Path,
    engine: str,
    style: str,
) -> None:
    """Write CSVs to Excel according to *style*."""
    if not csv_paths:
        sys.exit("No CSV files found.")

    with pd.ExcelWriter(output_path, engine=engine) as writer:
        if style == "tabs":
            # One worksheet per CSV (original behaviour)
            for path in csv_paths:
                sheet = sanitize_sheet_name(path.stem)
                print(f"→ Writing {path.name}  ➜  sheet “{sheet}”")
                df = pd.read_csv(path)
                df.to_excel(writer, sheet_name=sheet, index=False)

        elif style == "single":
            # All tables in one worksheet, with title row + 1 blank row between tables
            sheet = "combined"
            current_row = 0
            for path in csv_paths:
                df = pd.read_csv(path)
                print(f"→ Appending {path.name}  ➜  sheet “{sheet}” at row {current_row+1}")

                # Write the title row (filename, no extension)
                title_df = pd.DataFrame([[path.stem]])
                title_df.to_excel(
                    writer,
                    sheet_name=sheet,
                    index=False,
                    header=False,
                    startrow=current_row,
                )
                current_row += 1

                # Write the actual table
                df.to_excel(
                    writer,
                    sheet_name=sheet,
                    index=False,
                    header=True,
                    startrow=current_row,
                )
                current_row += len(df) + 2  # one blank row after each table

        else:
            raise ValueError(f"Unknown style: {style!r}")

    print(f"\nDone! Workbook saved to: {output_path.resolve()}")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
__version__ = "0.1.0"


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="combine_csv_to_excel.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Combine CSV files into an Excel workbook.",
        epilog=textwrap.dedent(
            """\
            Usage:
              python combine_csv_to_excel.py <directory> [-r] [-s {tabs,single}] [-o OUTPUT]

            Examples:
              # One worksheet per CSV (default)
              python combine_csv_to_excel.py data/

              # All tables on a single worksheet, separated by blank lines
              python combine_csv_to_excel.py data/ -s single -o all_tables.xlsx
            """
        ),
    )
    p.add_argument("directory", type=pathlib.Path, help="Directory with CSV files")
    p.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("combined.xlsx"),
        help="Output Excel file (default: combined.xlsx)",
    )
    p.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        help="Recurse into sub-directories when searching for CSVs",
    )
    p.add_argument(
        "-s",
        "--style",
        choices=["tabs", "single"],
        default="tabs",
        help="Workbook layout style: tabs (default) or single sheet",
    )
    p.add_argument(
        "--engine",
        choices=["xlsxwriter", "openpyxl"],
        default="xlsxwriter",
        help="Excel writer engine (xlsxwriter faster; openpyxl pure-Python)",
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return p


def main() -> None:
    args = make_parser().parse_args()

    if not args.directory.is_dir():
        sys.exit(f"Error: {args.directory} is not a directory.")

    csv_paths = find_csv_files(args.directory, args.recursive)
    csvs_to_excel(csv_paths, args.output, args.engine, args.style)


if __name__ == "__main__":
    main()
