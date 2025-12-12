#!/usr/bin/env python3
"""
Create Electrolytes.xlsx (unique SMILES) from a CSV with unknown/odd headers.

- Scans all cells for SMILES-like strings (conservative whitelist).
- Deduplicates.
- Writes an Excel file with a single column: 'rdkit_smiles' to feed scrape_py.py.

Usage:
  python make_electrolytes_xlsx.py run.csv  # writes Electrolytes.xlsx
"""
import sys
import re
from pathlib import Path
import pandas as pd

def extract_unique_smiles_from_csv(csv_path: Path) -> pd.Series:
    dfs = []
    # Try normal and headerless reads
    for kwargs in ({}, {"header": None}):
        try:
            dfs.append(pd.read_csv(csv_path, **kwargs))
        except Exception:
            pass
    if not dfs:
        raise RuntimeError(f"Could not read {csv_path} as CSV.")

    # Flatten all cells from all attempts
    cells = []
    for d in dfs:
        cells.extend(d.to_numpy().ravel().tolist())

    s = pd.Series(cells, dtype="string").str.strip()
    s = s[s.notna() & (s != "")]
    # Conservative SMILES-ish filter
    pattern = re.compile(r'^[A-Za-z0-9@\+\-\[\]\(\)\\\/=#%\.]+$')
    s = s[s.apply(lambda x: bool(pattern.fullmatch(x)))]
    return s.drop_duplicates().reset_index(drop=True)

def main():
    if len(sys.argv) < 2:
        print("Usage: python make_electrolytes_xlsx.py <input.csv>")
        sys.exit(1)
    inp = Path(sys.argv[1])
    out = Path("Electrolytes.xlsx")
    smiles = extract_unique_smiles_from_csv(inp)
    pd.DataFrame({"rdkit_smiles": smiles}).to_excel(out, index=False)
    print(f"Wrote {len(smiles)} unique SMILES to {out}")

if __name__ == "__main__":
    main()
