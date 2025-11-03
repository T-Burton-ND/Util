#!/bin/bash

# --- UGE Settings ---
#$ -pe smp 24
#$ -cwd
#$ -V
#$ -j y
#$ -l h_rt=84:00:00
#$ -N Se_Find_Shells
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

# --- Where to save the CSV (outside 'trial' so gzipping *.txt won't catch it) ---
OUTPUT_CSV="$START_DIR/shell_results.csv"
echo "run_name,electrolytes,run_number,shells" > "$OUTPUT_CSV"

# --- Enter the 'Run Outputs' directory ---
cd ./run_outputs || { echo "Failed to enter ./run_outputs"; exit 1; }

# --- Loop through subdirectories ---
for dir in run.* ; do
    if [[ -d "$dir" ]]; then
        echo "Processing directory: $dir"
        
        # Enter the directory
        cd "$dir" || { echo "Failed to enter directory $dir"; continue; }
        
        bname=${dir##*/}             # run.1_1.x2.r2
        core=${bname#run.}           # 1_1.x2.r2
        RUN_NAME=${core%%.x*}        # 1_1
        ELECTROLYTES=${RUN_NAME%%_*} # 1  (before _)
        RUN_NUMBER=${RUN_NAME##*_}   # 1  (after _)

        echo "Run name: $RUN_NAME | electrolytes: $ELECTROLYTES | run number: $RUN_NUMBER"


        # Unzip necessary files
        #[[ -f equilibrated.data.gz ]] && gunzip -f equilibrated.data.gz || { echo "⚠️ Missing equilibrated.data(.gz), skipping $dir"; cd ..; continue; }
        #[[ -f production.lammpstrj.gz ]] && gunzip -f production.lammpstrj.gz || { echo "⚠️ Missing production.lammpstrj(.gz), skipping $dir"; cd ..; continue; }


        # Call Python and capture its final line; assume it prints: "<RUN_NAME> <SHELLS>"
        PY_OUT=$(python /groups/bsavoie2/tburton2/Electrolytes/FF_carb/calc_scripts/cluster_py.py "$RUN_NAME" | tail -n 1)

        # Shell count is last whitespace-separated field on that line
        SHELLS=$(echo "$PY_OUT" | awk '{print $NF}')

        # Append to CSV
        echo "${RUN_NAME},${ELECTROLYTES},${RUN_NUMBER},${SHELLS}" >> "$OUTPUT_CSV"

        # Rezip files to save space
#        gzip ${RUN_NAME}.xyz
        #[[ -f equilibrated.data ]] && gzip -f equilibrated.data || echo "⚠️ Missing equilibrated.data, skipping gzip in $dir"
        #[[ -f production.lammpstrj ]] && gzip -f production.lammpstrj || echo "⚠️ Missing production.lammpstrj, skipping gzip in $dir"


        echo "Completed processing for $dir"

        # Go back up
        cd ..
    fi
done

# --- Track end time ---
END_TIME=$(date +%s)
END_STRING=$(date '+%Y-%m-%d %H:%M:%S')
ELAPSED=$((END_TIME - START_TIME))

echo "Job started at: $START_STRING"
echo "Job ended at:   $END_STRING"
echo "Elapsed time (seconds): $ELAPSED"

# --- Return to starting directory ---
cd "$START_DIR" || { echo "Failed to return to starting directory"; exit 1; }
