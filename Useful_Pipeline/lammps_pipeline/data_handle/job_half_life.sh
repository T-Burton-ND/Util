#!/bin/bash

# --- UGE Settings ---
#$ -pe smp 24
#$ -cwd
#$ -V
#$ -j y
#$ -l h_rt=84:00:00
#$ -N K100_Half_Life
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
OUTPUT_CSV="$START_DIR/half_life_results.csv"
echo "run_name,electrolytes,run_number,Mean Half Life,Std Deviation" > "$OUTPUT_CSV"

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

        echo "Run name: $RUN_NAME | electrolytes: $ELECTROLYTES | run number: $RUN_NUMBER" |

        python ../calc_scripts/shell_hl.py > python_output.txt 2>&1

        # Extract mean lifetime (text after "Mean Lifetime:")
        MEAN=$(grep "Mean Lifetime" python_output.txt | awk -F':' '{print $2}' | tr -d '[:space:]')

        # Extract std deviation (text after "Std:")
        STDDEV=$(grep "^Std" python_output.txt | awk -F':' '{print $2}' | tr -d '[:space:]')

        # Debug print
        echo "Parsed values → MEAN=$MEAN, STDDEV=$STDDEV"

        # Append to CSV
        echo "${RUN_NAME},${ELECTROLYTES},${RUN_NUMBER},${MEAN},${STDDEV}" >> "$OUTPUT_CSV"

        
        echo "Completed processing for $dir"

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
