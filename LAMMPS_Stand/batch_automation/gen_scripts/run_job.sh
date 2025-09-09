#!/bin/bash

# --- UGE Settings ---
#$ -pe smp 24
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
module load python 

# --- Get Run Name ---
RUN_NAME="$1"

# --- Run LAMMPS ---
: "${NSLOTS:=1}"
echo "Starting LAMMPS run for $RUN_NAME"
mpirun -np $NSLOTS lmp_mpi -in run.in -log log.lammps
echo "LAMMPS run completed for $RUN_NAME"

# --- Track end time ---
END_TIME=$(date +%s)
END_STRING=$(date '+%Y-%m-%d %H:%M:%S')
ELAPSED=$((END_TIME - START_TIME))

# --- Compute Diffusion Coefficient (your script) ---
echo "Running diffusion coefficient analysis..."
DIFFUSION_OUTPUT=$(python /temp180/bsavoie2/tburton2/Electrolytes/automated/calc_scripts/compute_diffusion.py msd_li.txt msd_cl.txt --outdir self_diff_out)
echo "$DIFFUSION_OUTPUT" | tee diffusion_results.txt

# ---- Analysis: MDAnalysis ----
# Prefer gzip’d trajectory if present
TRAJ="production.lammpstrj"
if [[ -f "${TRAJ}.gz" ]]; then
  TRAJ="${TRAJ}.gz"
fi

echo
echo "Using MDAnalysis to calculate diffusion coefficients on ${TRAJ} ..."
MD_DIFF_OUTPUTS=$(
  python /temp180/bsavoie2/tburton2/Electrolytes/automated/calc_scripts/compute_diffusion_mdanalysis.py \
    "${TRAJ}" \
    --dt-fs 0.5 --out MD_Analysis_out --fit-start-step 5000000 --fit-end-step 20000000 \
    | tee -a diffusion_results.txt
)


# --- Compute Viscosity (if applicable) ---
# Uncomment and modify the following lines if viscosity calculation is needed
# echo "Running viscosity analysis..."
# VISCOSITY_OUTPUT=$(python ../gen_scripts/compute_viscosity.py)
# echo "$VISCOSITY_OUTPUT" | tee viscosity_results.txt

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
  echo "=== Self-diffusion (msd_li/msd_cl) ==="
  echo "$DIFFUSION_OUTPUT"
  echo
  echo "=== MDAnalysis diffusion ==="
  echo "$MD_DIFF_OUTPUTS"
} > run_finished.flag

# --- Gzip the output files ---
gzip log.lammps
gzip msd_*.txt
gzip *.data
gzip *.lammpstrj