#!/bin/bash
#$ -N orca_array
#$ -cwd
#$ -pe smp 8
#$ -l h_rt=12:00:00
#$ -j y
#$ -V
#$ -t 1-4302 
#$ -tc 75
#$ -o logs/orca_array_$TASK_ID.log

module load orca
module load python

JOBS_DIR="run2_orca_jobs"
JOBLIST="$JOBS_DIR/job_dirs.txt"

# Get the directory corresponding to this array task
DIR=$(sed -n "${SGE_TASK_ID}p" "$JOBLIST")

if [ -z "$DIR" ]; then
    echo "No directory found for SGE_TASK_ID=$SGE_TASK_ID"
    exit 1
fi

echo "Task $SGE_TASK_ID using directory: $DIR"
cd "$DIR" || { echo "Failed to cd to $DIR"; exit 1; }

INPUT="job.inp"
OUTPUT="job.out"
SMIFILE="smiles.txt"
FLAG="enthalpy.flag"

# If already done, skip (lets you safely re-run arrays)
if [ -f "$FLAG" ]; then
    echo "Found $FLAG, skipping ORCA run."
    exit 0
fi

echo "Running ORCA on $INPUT with 8 cores..."
orca "$INPUT" > "$OUTPUT"

echo "ORCA job completed. Output written to $OUTPUT."
echo "-----------------------------------"

echo "Running Enthalpy Extraction Script..."
ENTHALPY=$(python /groups/bsavoie2/tburton2/yarp_networks/orca_enthalpies/calc_scripts/extract_enthalpy.py "$OUTPUT")

# grab SMILES string (first line of smiles.txt if it exists)
if [ -f "$SMIFILE" ]; then
    SMILES=$(head -n 1 "$SMIFILE")
else
    SMILES="UNKNOWN"
fi

# write "SMILES,enthalpy" to the flag file
echo "$SMILES,$ENTHALPY" > "$FLAG"
echo "Wrote $FLAG → $SMILES,$ENTHALPY"

# optional cleanup (you might want to keep .out/.inp during debug)
rm -f *.opt \
      *.gbw \
      *.hess \
      *.engrad \
      *.densities* \
      *.xyz \
      *.bibtex \
      *.inp \
      *.out

echo "Cleaned up intermediate files."
