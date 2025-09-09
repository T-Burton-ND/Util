#!/bin/bash
# ------ UGE settings ------
#$ -N HALL_MONITOR
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 1
#$ -l h_rt=72:00:00
#$ -M tburton2@nd.edu
#$ -m e   # email at job end (success OR fail). Use -m ae if you also want "at start".

# If you need a module or conda, do it here:
# module load python  # or:
# source ~/.bashrc
# conda activate rxn_graphs   # <- whatever env has Python 3

# ---- Configure via env vars (override per qsub -v ...) ----
export MON_ROOT_DIR="${MON_ROOT_DIR:-$PWD}"
export MON_RUN_GLOB="${MON_RUN_GLOB:-run.*}"
export MON_FLAG_NAME="${MON_FLAG_NAME:-run_finish.flag}"
export MON_CSV_PATH="${MON_CSV_PATH:-run_states.csv}"
export MON_INTERVAL_MIN="${MON_INTERVAL_MIN:-5}"
export MON_STOP_WHEN_ALL_DONE="${MON_STOP_WHEN_ALL_DONE:-1}"
export MON_MAX_HOURS="${MON_MAX_HOURS:-240}"
export MON_INITIAL_DELAY_MIN="${MON_INITIAL_DELAY_MIN:-1}"
export MON_LOG_GLOBS="${MON_LOG_GLOBS:-job_output_*.txt,*.out,log.lammps,*.log}"
export MON_STATE_PATH="${MON_STATE_PATH:-./data_handle/.monitor_state.json}"

# export MON_ERROR_FLAGS="run_error.flag,run_failed.flag"

# ---- Run ----
module load python
echo "[run_monitor] root=$MON_ROOT_DIR pattern=$MON_RUN_GLOB flag=$MON_FLAG_NAME interval=${MON_INTERVAL_MIN}min"
echo "[run_monitor] initial_delay=${MON_INITIAL_DELAY_MIN}min logs=${MON_LOG_GLOBS}"
python ./data_handle/monitor_runs.py   # adjust path if script lives in ./data_handle
rc=$?
echo "[run_monitor] exiting with code $rc"
exit $rc
