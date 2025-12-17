#!/bin/bash

#$ -N orca_batcher
#$ -cwd
#$ -pe smp 8

module load python

echo
echo "Starting ORCA job preparation..."
echo

CSV="${1:-smiles.csv}"          # default filename; change if needed
echo "Using SMILES CSV: $CSV"
JOBS_DIR="orca_jobs"
JOBLIST="job_dirs.txt"
mkdir -p "$JOBS_DIR"

# empty the list each time
: > "$JOBLIST"

echo
echo "==============================="
echo "Reading SMILES and preparing jobs..."
echo "==============================="

# Skip header if present; adjust if your CSV has no header
tail -n +2 "$CSV" | tr -d '\r' | while IFS=, read -r SMI; do
  [ -z "$SMI" ] && continue
  ID="$(printf "%s" "$SMI" | sha1sum | awk '{print substr($1,1,12)}')"

  echo
  echo "-------------------------------"
  echo "Processing SMILES: $SMI (ID: $ID)"
  DIR="$JOBS_DIR/$ID"
  mkdir -p "$DIR"
  echo "$SMI" > "$DIR/smiles.txt"

  # If we already have a completed enthalpy, skip
  if [ -f "$DIR/enthalpy.flag" ]; then
      echo "enthalpy.flag exists for $ID — skipping (already done)."
      continue
  fi

  echo "Generating XYZ and ORCA input files..."
  python ./calc_scripts/smiles2xyz.py "$SMI" -o "$DIR/molecule.xyz" || {
      echo "✗ Embed failed for $ID — skipping"
      echo "$SMI" > "$DIR/embed_failed.flag"
      continue
  }

  python ./calc_scripts/xyz2orca.py "$DIR/molecule.xyz" -o "$DIR/job.inp" || {
      echo "✗ ORCA input generation failed for $ID — skipping"
      echo "$SMI" > "$DIR/input_failed.flag"
      continue
  }

  # Add this directory to the job list
  echo "$DIR" >> "$JOBLIST"

done

cp $JOBLIST ./orca_jobs/ 

echo
echo "==============================="
echo "All jobs prepared."
echo "Job list written to: $JOBLIST"
echo "Number of jobs: $(wc -l < "$JOBLIST")"
echo "Next step: qsub a job array (see run_orca_array.sh)."
echo "==============================="
echo
