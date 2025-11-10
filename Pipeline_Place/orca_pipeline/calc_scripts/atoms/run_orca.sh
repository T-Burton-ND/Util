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

echo "Running ORCA on $INPUT with 8 cores..."
orca $INPUT > $OUTPUT
