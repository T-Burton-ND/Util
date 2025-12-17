#!/usr/bin/env bash
# gather_autocorr.sh
# Usage: ./gather_autocorr.sh /path/to/head_dir
# Collects LiLi_autocorr.csv files from each run.* directory
# and copies them into auto_corr_csvs/, renamed as <base>.csv
# Example: run.1_1.x2.r1/LiLi_autocorr.csv → auto_corr_csvs/1_1.r1.csv

set -euo pipefail
IFS=$'\n\t'

HEAD_DIR="${1:-run_outputs}"
OUT_DIR="${HEAD_DIR%/}/auto_corr_csvs"

mkdir -p "$OUT_DIR"

# Loop through all run.* directories directly under HEAD_DIR
find "$HEAD_DIR" -maxdepth 1 -mindepth 1 -type d -name 'run.*' | while read -r RUNDIR; do
  run_base="$(basename "$RUNDIR")"

  # Extract base name: run.1_1.x2.r1 → 1_1.r1
  base_name="$(echo "$run_base" | sed -E 's/^run\.//; s/\.x[0-9]+//g')"

  src_file="$RUNDIR/LiLi_autocorr.csv"
  dest_file="$OUT_DIR/${base_name}.csv"

  if [[ -f "$src_file" ]]; then
    cp -p "$src_file" "$dest_file"
    echo "Copied $src_file → $dest_file"
  else
    echo "⚠️  Missing LiLi_autocorr.csv in $RUNDIR (skipping)"
  fi
done
