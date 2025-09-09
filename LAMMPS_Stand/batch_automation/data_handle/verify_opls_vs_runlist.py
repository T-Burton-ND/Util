#!/usr/bin/env python3
"""
Verify that all molecules referenced in run_list.csv have matching .lmp files.

Assumptions:
- run_list.csv has a header; first column is RunID, remaining columns are molecule names.
- OPLS files are named <molecule>.lmp (exact match on the stem by default).

Usage:
  python verify_opls_vs_runlist.py --csv run_list.csv --opls-dir ./opls_files
  # or if you have a text list of files instead of a directory:
  python verify_opls_vs_runlist.py --csv run_list.csv --opls-list opls_files.txt

Outputs:
- A clear report to stdout.
- Optional CSV with missing items: --report missing_report.csv
- Exits with code 0 if all good; 1 if any missing.
"""

import argparse
import csv
import os
import sys
import glob
from difflib import get_close_matches
from collections import defaultdict

def load_opls_stems_from_dir(opls_dir: str) -> set:
    files = glob.glob(os.path.join(opls_dir, "*.lmp"))
    return {os.path.splitext(os.path.basename(p))[0] for p in files}

def load_opls_stems_from_list(list_path: str) -> set:
    stems = set()
    with open(list_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            base = os.path.basename(s)
            if base.endswith(".lmp"):
                stems.add(os.path.splitext(base)[0])
            else:
                # If user listed bare stems or other files, add smartly
                stem = os.path.splitext(base)[0]
                stems.add(stem)
    return stems

def normalize_name(s: str) -> str:
    # Minimal normalization matching your batch script expectations
    return s.strip()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", help="Path to run_list.csv", default="./run_list.csv")
    ap.add_argument("--opls-dir", help="Directory containing *.lmp files", default="./opls_files")
    ap.add_argument("--opls-list", help="Text file listing OPLS files (one per line)")
    ap.add_argument("--report", help="Optional path to write a CSV of missing items")
    ap.add_argument("--case-insensitive", action="store_true",
                    help="Try case-insensitive matching before declaring missing")
    args = ap.parse_args()

    # Load OPLS stems
    if args.opls_dir:
        if not os.path.isdir(args.opls_dir):
            print(f"[ERROR] OPLS directory not found: {args.opls_dir}", file=sys.stderr)
            sys.exit(2)
        opls_stems = load_opls_stems_from_dir(args.opls_dir)
    else:
        if not os.path.isfile(args.opls_list):
            print(f"[ERROR] OPLS list file not found: {args.opls_list}", file=sys.stderr)
            sys.exit(2)
        opls_stems = load_opls_stems_from_list(args.opls_list)

    opls_stems_ci_map = {s.lower(): s for s in opls_stems}  # for case-insensitive fallback

    # Parse run_list.csv
    if not os.path.isfile(args.csv):
        print(f"[ERROR] run list CSV not found: {args.csv}", file=sys.stderr)
        sys.exit(2)

    missing = []  # tuples: (RunID, molecule_name, suggestion)
    per_run_missing = defaultdict(list)
    all_molecules = set()

    with open(args.csv, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header:
            print("[ERROR] CSV is empty.", file=sys.stderr)
            sys.exit(2)

        for row in reader:
            if not row:
                continue
            run_id = normalize_name(row[0])
            # Remaining columns are molecules (may include blanks)
            molecules = [normalize_name(x) for x in row[1:] if normalize_name(x)]
            for mol in molecules:
                all_molecules.add(mol)
                if mol in opls_stems:
                    continue
                # Try case-insensitive fallback
                suggestion = ""
                if args.case_insensitive:
                    lower = mol.lower()
                    if lower in opls_stems_ci_map:
                        # treat as matched; warn, but not missing
                        continue
                # Fuzzy suggestion
                suggestion_candidates = get_close_matches(mol, opls_stems, n=3, cutoff=0.6)
                if suggestion_candidates:
                    suggestion = "; ".join(suggestion_candidates)
                missing.append((run_id, mol, suggestion))
                per_run_missing[run_id].append((mol, suggestion))

    # Summary
    print("====================================================")
    print("OPLS verification report")
    print("====================================================")
    print(f"Total unique OPLS stems found: {len(opls_stems)}")
    print(f"Total unique molecules referenced in CSV: {len(all_molecules)}")

    if missing:
        print("\nMISSING FILES:")
        for run_id, mol, sugg in missing:
            if sugg:
                print(f"  - Run {run_id}: '{mol}.lmp' NOT FOUND. Suggestions: {sugg}")
            else:
                print(f"  - Run {run_id}: '{mol}.lmp' NOT FOUND.")
        print("\nPer-run missing counts:")
        for run_id, items in per_run_missing.items():
            print(f"  - {run_id}: {len(items)} missing")
    else:
        print("\nAll molecules have corresponding .lmp files. ✅")

    # Optional CSV report
    if args.report:
        import csv as _csv
        with open(args.report, "w", newline="") as rf:
            w = _csv.writer(rf)
            w.writerow(["RunID", "Molecule", "Expected_File", "Suggestions"])
            for run_id, mol, sugg in missing:
                w.writerow([run_id, mol, f"{mol}.lmp", sugg])
        print(f"\nMissing report written to: {args.report}")

    # Exit code indicates success/failure
    sys.exit(0 if not missing else 1)

if __name__ == "__main__":
    main()
