#!/bin/bash

# ==============================
# Script for generating MSD s with MDA.
# Single Job for All Runs
# ==============================
# ------ UGE settings ------
#$ -N H_MSD_GEN
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 4
#$ -l h_rt=240:00:00
#$ -M tburton2@nd.edu
#$ -m e   # email at job end (success OR fail). Use -m ae if you also want "at start".

set -euo pipefail
shopt -s nullglob  # so 'run.*' with no matches doesn't iterate the literal

cleanup() {
      [[ -f production.lammpstrj ]] && gzip -f production.lammpstrj
      [[ -f equilibrated.data   ]] && gzip -f equilibrated.data
    }

# cd into the Outputs Super-Directory
RUN_ROOT="${RUN_ROOT:-./run_outputs}"
cd "$RUN_ROOT"

for dir in run.*; do
    [[ -d "$dir" ]] || continue
    echo "Processing directory: $dir"
    cd "$dir"
    CWD_DIR=$(pwd)
    echo "Current working directory: $CWD_DIR"
    
    run="${dir##*/}"              # e.g., run.l_1_1.x2.r1
    base="${run#*.}"              # strip first dot -> l_1_1.x2.r1
    base="${base%%.*}"            # take up to next dot -> l_1_1

    echo "Calculating MSD for run: $run (base: $base)"    
    #[[ -f equilibrated.data.gz ]] && gunzip -f equilibrated.data.gz || { echo "⚠️ Missing equilibrated.data(.gz), skipping $dir"; cd ..; continue; }
    #[[ -f production.lammpstrj.gz ]] && gunzip -f production.lammpstrj.gz || { echo "⚠️ Missing production.lammpstrj(.gz), skipping $dir"; cd ..; continue; }

    #trap cleanup EXIT
    python /groups/bsavoie2/tburton2/Electrolytes/FF_carb/calc_scripts/msd_py.py "$run"
    #[[ -f equilibrated.data ]] && gzip -f equilibrated.data || echo "⚠️ Missing equilibrated.data, skipping gzip in $dir"
    #[[ -f production.lammpstrj ]] && gzip -f production.lammpstrj || echo "⚠️ Missing production.lammpstrj, skipping gzip in $dir"

    cd ..
done

cd ..
echo "MSD calculations completed for all runs."

# --- Collect all MSD output files ---
COLLECT_DIR="$RUN_ROOT/MSD_Results"

# Safely make the directory if it doesn’t exist
mkdir -p "$COLLECT_DIR"

# Copy results, prefixing with run directory name to avoid collisions
# Copy each matching file as-is (same name)
find "$RUN_ROOT" -type f -name "mda_msd_run.*" -exec cp -vn {} "$COLLECT_DIR"/ \;

echo "✅ Copied all MSD output files (names unchanged) into: $COLLECT_DIR"

echo "Renaming collected MSD files..."

for file in "$COLLECT_DIR"/*; do
    [[ -f "$file" ]] || continue

    # Extract run directory name
    fname=$(basename "$file")
    echo "Renaming file: $fname"
    core="${fname#mda_msd_run.}"   # -> 1_1.x2.r1.txt
    echo "Core before cleanup: $core"
    core="${core%.txt}"              # -> 1_1.x2.r1
    echo "Core after removing .txt: $core"
    core="${core/.x2/}"              # -> 1_1.r1
    echo "Core after removing .x2: $core"
    core="${core//.r/_}"              # -> 1_1_1
    echo "Core after replacing .r with _: $core"
    new_name="msd_${core}.txt"
    echo "New name: $new_name"

    mv -v "$file" "$COLLECT_DIR/$new_name"
done

echo "✅ Renamed all MSD files in $COLLECT_DIR to msd_{run}.txt format."
