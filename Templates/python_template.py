#!/usr/bin/env python3
"""
<file_name>.py — <one‑line summary>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Short description (2–3 sentences max).

Usage
-----
python <file_name>.py [OPTIONS] <positional-args>

Examples
--------
python <file_name>.py -i input.csv -o out.xlsx

Dependencies
------------
pandas, numpy  # (only list non-stdlib libs you really use)
"""

from __future__ import annotations
import argparse
# ... (imports)

__author__  = "Thomas J. Burton"
__version__ = "0.1.0"

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Short, human-readable sentence about what this script does."
    )
    p.add_argument("-i", "--input", required=True, help="Path to input file.")
    p.add_argument("-o", "--output", default="output.ext", help="Where to write results.")
    # add more args...
    return p

def main(args=None) -> int:
    ns = build_parser().parse_args(args=args)
    # --- your logic ---
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
