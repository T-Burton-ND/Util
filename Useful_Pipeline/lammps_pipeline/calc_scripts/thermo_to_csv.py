#!/usr/bin/env python3
"""
Convert a whitespace-separated LAMMPS thermo text file into a CSV.

Usage:
    python thermo_to_csv.py thermo.txt
Output:
    thermo.csv  (same base name)
"""

import sys
import pandas as pd

def main():
    if len(sys.argv) != 2:
        print("Usage: python thermo_to_csv.py input_file.txt")
        sys.exit(1)

    infile = sys.argv[1]
    outfile = infile.rsplit(".", 1)[0] + ".csv"

    # Load as whitespace-delimited
    df = pd.read_csv(infile, sep=r"\s+", engine="python", comment="#")

    # Drop empty columns (if any)
    df = df.dropna(axis=1, how="all")

    # Write CSV
    df.to_csv(outfile, index=False)

    print(f"✔ Converted '{infile}' → '{outfile}'")
    print(f"✔ Columns detected: {list(df.columns)}")

if __name__ == "__main__":
    main()
