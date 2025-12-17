#!/bin/bash

# --- UGE Settings ---
#$ -pe smp 24
#$ -cwd
#$ -V
#$ -j y
#$ -l h_rt=84:00:00
#$ -N K100_Auto_Corr
#$ -m e
#$ -M tburton2@nd.edu

# --- Save current directory ---
START_DIR=$(pwd)

# --- Track start time ---
START_TIME=$(date +%s)
START_STRING=$(date '+%Y-%m-%d %H:%M:%S')

# --- Load Modules ---
source ~/.bashrc
conda activate base
module load python


# --- Enter the 'Run Outputs' directory ---
cd ./run_outputs || { echo "Failed to enter ./run_outputs"; exit 1; }

mkdir -p auto_correlations/

# --- Loop through subdirectories ---
for dir in run.* ; do
    if [[ -d "$dir" ]]; then
        echo "==== Processing directory: $dir ===="
        
        # Enter the directory
        cd "$dir" || { echo "Failed to enter directory $dir"; continue; }
        
        bname=${dir##*/}             # run.1_1.x2.r2
        core=${bname#run.}           # 1_1.x2.r2
        RUN_NAME=${core%%.x*}        # 1_1
        ELECTROLYTES=${RUN_NAME%%_*} # 1  (before _)
        RUN_NUMBER=${RUN_NAME##*_}   # 1  (after _)

        echo "==== Running Python Script ===="
        # Call Python and capture its final line; assume it prints: "<RUN_NAME> <SHELLS>"
        python ../calc_scripts/analysis_py.py "$RUN_NAME"

        echo "==== Copying output to auto_correlations directory ===="
        cp Li_all_analysis.png ../auto_correlations/$RUN_NAME.png

        echo "==== Completed processing for $dir ===="

        # Go back up
        cd ..
    fi
done

# --- Track end time ---
END_TIME=$(date +%s)
END_STRING=$(date '+%Y-%m-%d %H:%M:%S')
ELAPSED=$((END_TIME - START_TIME))

echo "==== Job started at: $START_STRING ===="
echo "==== Job ended at:   $END_STRING ===="
echo "==== Elapsed time (seconds): $ELAPSED ===="

# --- Return to starting directory ---
cd "$START_DIR" || { echo "Failed to return to starting directory"; exit 1; }
