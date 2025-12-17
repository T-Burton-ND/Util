#!/bin/bash

# ==============================
# Batch LAMMPS Automation Script
# ==============================

set -euo pipefail

# --------- Where to put all run.* directories ---------
RUN_ROOT="${RUN_ROOT:-run_outputs}"
mkdir -p "$RUN_ROOT"

# --------- Paths ---------
YAML_FILE="config.yaml"
GEN_DIR="$(cd ./gen_scripts && pwd)"
OPLS_DIR="$(cd ./opls_files && pwd)"
CSV_FILE="run_list.csv"
LAMMPS_JOB_TEMPLATE="run_job.sh"

# --------- Replicate settings (override with env vars) ---------
REPL_FACTORS=(${REPL_FACTORS:-2 2 2})   # e.g., REPL_FACTORS="2 4 5" ./batch.sh
ION_COUNT=${ION_COUNT:-1}              # each ion species gets this many molecules
BASE_SEED="${BASE_SEED:-4928459}"       # per-replicate seed = BASE_SEED + index
TOTAL_ELEC="${TOTAL_ELEC:-40}"         # electrolytes (excluding ions) sum to this

# NEW: how many ion species are in the CSV row, from the left
# 1 => CSV is: run_id, cation, electrolyte_1, electrolyte_2, ...
# 2 => CSV is: run_id, cation, anion, electrolyte_1, electrolyte_2, ...
ION_SPECIES="${ION_SPECIES:-2}"
if ! [[ "$ION_SPECIES" =~ ^[12]$ ]]; then
  echo "ERROR: ION_SPECIES must be 1 or 2 (got '$ION_SPECIES')."
  exit 1
fi

# --------- Tiny YAML edit helper (no yq needed) ---------
yaml_set_key () {
  local in_yaml="$1" out_yaml="$2" key="$3" value="$4"
  python - "$in_yaml" "$out_yaml" "$key" "$value" <<'PY'
import sys, yaml
inp, outp, key, val = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
with open(inp, 'r') as f:
    data = yaml.safe_load(f) or {}
def parse(v):
    for caster in (int, float):
        try: return caster(v)
        except: pass
    if isinstance(v, str) and v.lower() in ("true","false"):
        return v.lower()=="true"
    return v
data[key] = parse(val)
with open(outp, 'w') as f:
    yaml.safe_dump(data, f, sort_keys=False)
PY
}

trim() { awk '{$1=$1;print}'; }

# --------- CSV + summary ---------
echo
echo "Starting batch processing from $CSV_FILE..."
echo "RunID,RepFactor,RepIndex,Seed,Molecules,MoleculeCounts,SubmissionTime" > run_summary.csv.tmp
echo >> "$CSV_FILE"   # ensure trailing newline

# ---------- Start Hall Monitor -------
echo 
echo "Starting Hall Monitor"
qsub -v MON_ROOT_DIR="$PWD/$RUN_ROOT" ./data_handle/hall_monitor.sh || true
echo 

# --------- Process each row ---------
tail -n +2 "$CSV_FILE" | while IFS=',' read -r RUN_ID rest_of_line; do
    RUN_ID="$(echo "$RUN_ID" | tr -d '\r' | trim)"

    # Collect molecule names (order matters)
    IFS=',' read -ra RAW_FILE_NAMES <<< "$rest_of_line"
    MOLECULE_NAMES=()
    for name in "${RAW_FILE_NAMES[@]}"; do
        clean_name="$(echo "$name" | tr -d '\r' | trim)"
        [[ -n "$clean_name" ]] && MOLECULE_NAMES+=("$clean_name")
    done

    # Skip if none
    if [ ${#MOLECULE_NAMES[@]} -eq 0 ]; then
        echo
        echo " ----------------------------------------"
        echo "Run $RUN_ID has no molecules listed. Skipping."
        echo " ----------------------------------------"
        continue
    fi

    # Validate presence of at least one ion species + one electrolyte
    if (( ${#MOLECULE_NAMES[@]} < ION_SPECIES + 1 )); then
        echo
        echo "----------------------------------------"
        echo "Run $RUN_ID — not enough molecules for requested ION_SPECIES=$ION_SPECIES"
        echo "Need at least: 1 run_id + $ION_SPECIES ion(s) + 1 electrolyte."
        echo "Skipping."
        echo "----------------------------------------"
        continue
    fi

    # Partition into ions vs electrolytes
    ION_NAMES=("${MOLECULE_NAMES[@]:0:ION_SPECIES}")
    ELEC_NAMES=("${MOLECULE_NAMES[@]:ION_SPECIES}")

    echo
    echo "----------------------------------------"
    echo "Run $RUN_ID"
    echo "ION_SPECIES: $ION_SPECIES"
    echo "Ions: ${ION_NAMES[*]}"
    echo "Electrolytes: ${ELEC_NAMES[*]}"
    echo "Replicate factors: ${REPL_FACTORS[*]}"
    echo "----------------------------------------"
    echo

    # --------- Step 1: Prepare .lmp Files (copy into GEN_DIR as your script expects) ---------
    cd "$GEN_DIR"

    LMP_FILES_TO_PASS=""
    for mol in "${MOLECULE_NAMES[@]}"; do
        lmp_file="${mol}.lmp"
        if [ ! -f "$OPLS_DIR/$lmp_file" ]; then
            echo "ERROR: Missing OPLS file: $OPLS_DIR/$lmp_file"
            echo "Skipping run $RUN_ID."
            rm -f *.lmp || true
            cd ..
            continue 2
        fi
        cp "$OPLS_DIR/$lmp_file" .         # keep .lmp extension; OPLS_box.py expects it
        LMP_FILES_TO_PASS+="$lmp_file "
    done
    LMP_FILES_TO_PASS="$(echo "$LMP_FILES_TO_PASS" | xargs)"  # trim trailing space

    # --------- Step 2: Compute Molecule Counts (ion(s) fixed, electrolytes distributed) ---------
    MOLECULE_COUNTS=""

    # Ion species counts
    for _ in "${ION_NAMES[@]}"; do
        MOLECULE_COUNTS+="$ION_COUNT "
    done

    # Electrolyte distribution
    ELEC_COUNT=${#ELEC_NAMES[@]}
    if (( ELEC_COUNT > 0 )); then
        total_elec="$TOTAL_ELEC"
        if (( total_elec < ELEC_COUNT )); then
            echo "NOTE: TOTAL_ELEC ($total_elec) < num electrolytes ($ELEC_COUNT), bumping to $ELEC_COUNT to keep ≥1 each."
            total_elec="$ELEC_COUNT"
        fi
        base=$(( total_elec / ELEC_COUNT ))
        rem=$(( total_elec % ELEC_COUNT ))

        for i in $(seq 1 "$ELEC_COUNT"); do
            c="$base"
            if (( i <= rem )); then c=$((c+1)); fi
            MOLECULE_COUNTS+="$c "
        done
    fi

    MOLECULE_COUNTS="$(echo "$MOLECULE_COUNTS" | xargs)"  # trim

    echo "Final Molecule Counts: $MOLECULE_COUNTS"
    echo 

    # --------- Step 3: Run OPLS_box.py ONCE ---------
    OUTPUT_DATA="system.data"
    BOX_SIZE=50

    echo "Running OPLS_box.py with files: \"$LMP_FILES_TO_PASS\""
    python OPLS_box.py "$LMP_FILES_TO_PASS" -N "$MOLECULE_COUNTS" -O "$OUTPUT_DATA" -L "$BOX_SIZE"

    # Check outputs
    if [ ! -f "$OUTPUT_DATA" ] || [ ! -f "settings.lmp" ]; then
        echo "Run $RUN_ID failed during OPLS_box.py execution!"
        rm -f *.lmp
        cd ..
        continue
    fi

    # Move generated files back to project root
    mv "$OUTPUT_DATA" settings.lmp ../

    # Clean up copied .lmp files from GEN_DIR
    rm -f *.lmp
    cd ..

    # --------- Step 4: For each replicate size: make a folder, override YAML, gen run.in, submit ---------
    rep_idx=0
    declare -A RF_COUNT=()
    
    for rf in "${REPL_FACTORS[@]}"; do
        rep_idx=$((rep_idx + 1))
        vseed=$((BASE_SEED + rep_idx))

        RF_COUNT["$rf"]=$(( ${RF_COUNT["$rf"]:-0} + 1 ))
        rord=${RF_COUNT["$rf"]}
        RUN_DIR="${RUN_ROOT}/run.${RUN_ID}.x${rf}.r${rord}"
        mkdir -p "$RUN_DIR"

        # Copy inputs into replicate folder
        cp system.data settings.lmp "$RUN_DIR/"

        # Create a temp yaml per replicate with overridden replicate & vseed
        TMP_YAML="${RUN_DIR}/lammps_set_${RUN_ID}_x${rf}.yaml"
        cp "$YAML_FILE" "$TMP_YAML"
        yaml_set_key "$TMP_YAML" "$TMP_YAML" replicate "$rf"
        yaml_set_key "$TMP_YAML" "$TMP_YAML" vseed "$vseed"

        # Generate run.in inside the replicate folder
        python "$GEN_DIR/gen_in.py" "$TMP_YAML"
        if [ ! -f "${RUN_DIR}/run.in" ]; then
            # gen_in.py may have written run.in in CWD; move it if needed
            if [ -f "run.in" ]; then mv run.in "$RUN_DIR/"; fi
        fi
        if [ ! -f "${RUN_DIR}/run.in" ]; then
            echo "Run $RUN_ID (x${rf}) failed during LAMMPS input generation!"
            continue
        fi

        # Job script and submit
        cp "$GEN_DIR/$LAMMPS_JOB_TEMPLATE" "$RUN_DIR/"
        (
          cd "$RUN_DIR"
          qsub -N "run_${RUN_ID}_x${rf}_r${rord}" "$LAMMPS_JOB_TEMPLATE" "$RUN_ID"  
        )
        # Log submission
        TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
        MOLECULE_LIST=$(IFS=';'; echo "${MOLECULE_NAMES[*]}")
        echo "$RUN_ID,$rf,$rep_idx,$vseed,\"$MOLECULE_LIST\",\"$MOLECULE_COUNTS\",$TIMESTAMP" >> run_summary.csv.tmp
    done

    # Keep base files for reference or clean them—your call. If you want them gone:
    rm -f system.data settings.lmp

done

echo
echo "Batch processing completed."

# Count lines with a timestamp (skip header)
SUCCESS_COUNT=$(tail -n +2 run_summary.csv.tmp | awk -F',' 'NF>0 && $NF ~ /^[0-9]{4}-[0-9]{2}-[0-9]{2}/ {count++} END{print count+0}')

echo
echo "Total jobs submitted successfully: $SUCCESS_COUNT"
echo

END_TIME=$(date '+%Y-%m-%d %H:%M:%S')

mkdir -p summaries
mv run_summary.csv.tmp "summaries/$(date '+%Y-%m-%d_%H-%M-%S')_run_summary.csv"
