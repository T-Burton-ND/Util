#!/bin/bash
# batch.sh — read SMILES (single-column CSV), make per-SMILES dirs, build ORCA jobs, submit

module load python

echo
echo "Starting batch ORCA job submission..."
echo

CSV="${1:-hp_smi.csv}"
echo "Using SMILES CSV: $CSV"
JOBS_DIR="hp_orca_jobs"
echo "Jobs will be created in directory: $JOBS_DIR"
mkdir -p "$JOBS_DIR"

echo
echo "==============================="
echo "Reading SMILES and preparing jobs..."
echo "================================"

submitted=0

# Skip header if present; adjust if your CSV has no header
tail -n +2 "$CSV" | tr -d '\r' | while IFS=, read -r SMI; do
  [ -z "$SMI" ] && continue
  ID="$(printf "%s" "$SMI" | sha1sum | awk '{print substr($1,1,12)}')"
  echo
  echo
  echo "-------------------------------"
  echo "Processing SMILES: $SMI (ID: $ID)"
  DIR="$JOBS_DIR/$ID"
  mkdir -p "$DIR"
  echo "$SMI" > "$DIR/smiles.txt"

  # Make XYZ and ORCA input (add Freq so enthalpy is printed)
  echo "Generating XYZ and ORCA input files..."
  python ./calc_scripts/smiles2xyz.py "$SMI" -o "$DIR/molecule.xyz" || { echo "✗ Embed failed for $ID — skipping"; echo "$SMI" > "$DIR/embed_failed.flag"; continue; }

  python ./calc_scripts/xyz2orca.py "$DIR/molecule.xyz" -o "$DIR/job.inp" 

  # Copy submit script and submit from inside the job dir
  echo "Submitting ORCA job..."
  cp ./templates/run_orca.sh "$DIR/"
  (cd "$DIR" && qsub run_orca.sh)
  ((submitted++))

done

echo
echo "==============================="
echo "All jobs processed."
echo "Submitted: $submitted"
echo "==============================="
echo