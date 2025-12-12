#!/bin/bash
#$ -N orca_job
#$ -cwd
#$ -pe smp 8
#$ -l h_rt=12:00:00
#$ -j y
#$ -V

module load orca


INPUT="job.inp"
OUTPUT="job.out"
SMIFILE="smiles.txt"
FLAG="enthalpy.flag"

echo "Running ORCA on $INPUT with 8 cores..."
orca $INPUT > $OUTPUT

echo "ORCA job completed. Output written to $OUTPUT."
echo "-----------------------------------"

echo "Running Enthalpy Extraction Script..."
ENTHALPY=$(python /groups/bsavoie2/tburton2/yarp_again_in_depth/orca_enthalpies/calc_scripts/extract_enthalpy.py "$OUTPUT")

# grab SMILES string (first line of smiles.txt if it exists)
if [ -f "$SMIFILE" ]; then
    SMILES=$(head -n 1 "$SMIFILE")
else
    SMILES="UNKNOWN"
fi

# write "SMILES,enthalpy" to the flag file
echo "$SMILES,$ENTHALPY" > "$FLAG"
echo "Wrote $FLAG → $SMILES,$ENTHALPY"

rm *.opt
rm *.gbw
rm *.tmp
rm *.hess
rm *.engrad
rm *.densities
echo "Cleaned up intermediate files."