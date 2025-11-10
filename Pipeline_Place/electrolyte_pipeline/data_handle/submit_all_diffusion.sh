#!/usr/bin/env bash
# Submit diffusion jobs in every run dir

set -euo pipefail

ROOT="/temp180/bsavoie2/tburton2/Electrolytes/automated"
RUNS_DIR="${ROOT}/run_outputs"
JOB_SCRIPT="${ROOT}/data_handle/diffusion_job.sh"

if [[ ! -f "$JOB_SCRIPT" ]]; then
  echo "[ERROR] Missing job script: $JOB_SCRIPT"
  exit 1
fi

echo "[submit_all_diffusion] Scanning $RUNS_DIR ..."
mapfile -t run_dirs < <(find "$RUNS_DIR" -maxdepth 1 -type d -name "run.*" | sort)

for rd in "${run_dirs[@]}"; do
  echo "Submitting diffusion job in $rd"
  cp "$JOB_SCRIPT" "$rd/"
  (cd "$rd" && qsub diffusion_job.sh)
done
