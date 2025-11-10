#!/bin/bash
# collect_flags.sh — recursively collect enthalpy.flag files into one CSV

ROOT_DIR="${1:-da_orca_jobs}"   # pass directory as arg or use default
OUT_FILE="all_enthalpies.csv"

echo "SMILES,enthalpy_kcalmol" > "$OUT_FILE"

find "$ROOT_DIR" -type f -name "enthalpy.flag" | sort | while read -r f; do
  cat "$f"
done >> "$OUT_FILE"

echo "✅ Collected $(wc -l < "$OUT_FILE") lines into $OUT_FILE"
ß