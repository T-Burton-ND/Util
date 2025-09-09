# --- UGE Settings ---
#$ -pe smp 12
#$ -cwd
#$ -V
#$ -j y
#$ -o job_output_$JOB_ID.txt
#$ -l h_rt=84:00:00

set -euo pipefail

echo "[diffusion_job] Running in $(pwd)"

# --- Compute Diffusion Coefficient (your script) ---
#echo "Running diffusion coefficient analysis..."
#DIFFUSION_OUTPUT=$(python ../calc_scripts/compute_diffusion.py msd_li.txt msd_cl.txt --outdir self_diff_out)
#echo "$DIFFUSION_OUTPUT" | tee diffusion_results.txt

# --- Use MDAnalysis to get Diffusion Coefficients ---
echo
echo "Using MDAnalysis to calculate diffusion coefficients..."
if [[ -f production.lammpstrj.gz ]]; then
  traj="production.lammpstrj.gz"
else
  traj="production.lammpstrj"
fi

python ../calc_scripts/compute_diffusion_mdanalysis.py "$traj" \
    --dt-fs 0.5 --out MD_Analysis_out --fit-start-step 5000000 --fit-end-step 20000000 \
    | tee -a diffusion_results.txt

echo "[diffusion_job] Done."
