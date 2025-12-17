#!/bin/bash

# --- UGE Settings ---
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -j y
#$ -o job_output_$JOB_ID.txt
#$ -l h_rt=84:00:00

# --- Track start time ---
START_TIME=$(date +%s)
START_STRING=$(date '+%Y-%m-%d %H:%M:%S')

# --- Load Modules ---
source ~/.bashrc  
conda activate base
module load intel
module load intelmpi
module load lammps

# --- Get Run Name ---
RUN_NAME="$1"

# --- Run LAMMPS ---
: "${NSLOTS:?NSLOTS not set by UGE}"     # fail fast if something is off
export OMP_NUM_THREADS=1                 # pure-MPI: 1 thread per rank

echo "Starting LAMMPS run for $RUN_NAME"
mpirun -np $NSLOTS lmp_mpi -in run.in -log log.lammps
echo "LAMMPS run completed for $RUN_NAME"

# --- Track end time ---
END_TIME=$(date +%s)
END_STRING=$(date '+%Y-%m-%d %H:%M:%S')
ELAPSED=$((END_TIME - START_TIME))

# --- Write flag + D result ---
{
  echo "Job started at: $START_STRING"
  echo "Job ended at:   $END_STRING"
  echo "Total wall time: ${ELAPSED} seconds"
  echo "Hostname: $(hostname)"
  echo "UGE Job ID: $JOB_ID"
  echo "Run Name: $RUN_NAME"
  echo
} > run_finished.flag
