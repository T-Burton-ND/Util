#!/bin/bash
#$ -N orca_job
#$ -cwd
#$ -pe smp 8
#$ -l h_rt=12:00:00
#$ -j y
#$ -V

module load orca

for f in *.inp; do
  dir="${f%.inp}"
  mkdir -p "$dir"
  cp "$f" "$dir/job.inp"
  cp run_orca.sh "$dir/"
  (cd "$dir" && qsub run_orca.sh)
done
