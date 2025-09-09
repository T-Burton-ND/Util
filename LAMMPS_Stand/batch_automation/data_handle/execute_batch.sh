#!/bin/bash

# ==============================
# Batch LAMMPS Automation Script
# ==============================

# --------- Paths ---------
YAML_FILE="lammps_set.yaml"
GEN_DIR="$(cd ./gen_scripts && pwd)"
OPLS_DIR="$(cd ./opls_files && pwd)"
CSV_FILE="run_list.csv"
LAMMPS_JOB_TEMPLATE="run_job.sh"


# --------- Read CSV and Loop Over Runs ---------
echo "Starting batch processing from $CSV_FILE..."
echo "RunID,Molecules,MoleculeCounts,SubmissionTime,JobID" > run_summary.csv
cp run_list.csv exec_list.csv  # Backup original CSV
echo >> run_list.csv


tail -n +2 $CSV_FILE | while IFS=',' read -r RUN_ID rest_of_line; do
    RUN_ID=$(echo $RUN_ID | tr -d '\r' | xargs)

    # Parse the rest of the line into an array of molecule names
    IFS=',' read -ra RAW_FILE_NAMES <<< "$rest_of_line"

    MOLECULE_NAMES=()
    for name in "${RAW_FILE_NAMES[@]}"; do
        clean_name=$(echo $name | xargs | tr -d '\r')
        if [ -n "$clean_name" ]; then
            MOLECULE_NAMES+=("$clean_name")
        fi
    done

    # Skip if no molecules found
    if [ ${#MOLECULE_NAMES[@]} -eq 0 ]; then
        echo " "
        echo " ----------------------------------------"
        continue
    fi
    echo " "
    echo "----------------------------------------"
    echo "Run $RUN_ID"
    echo "Molecule Names: ${MOLECULE_NAMES[@]}"
    echo "----------------------------------------"
    echo " "
    RUN_DIR="run.${RUN_ID}"
    mkdir -p $RUN_DIR

    # --------- Step 1: Prepare .lmp Files ---------
    cd $GEN_DIR

    LMP_FILES_TO_PASS=""
    for mol in "${MOLECULE_NAMES[@]}"; do
        lmp_file="${mol}.lmp"
        cp "$OPLS_DIR/$lmp_file" .
        LMP_FILES_TO_PASS+="$lmp_file "
    done
    LMP_FILES_TO_PASS=$(echo $LMP_FILES_TO_PASS | xargs)  # Trim trailing space

    # --------- Step 2: Compute Molecule Counts ---------
    MOLECULE_COUNT=${#MOLECULE_NAMES[@]}
    MOLECULE_COUNTS="3 3"

    if [ $MOLECULE_COUNT -gt 2 ]; then
        REMAINING_FILES=$((MOLECULE_COUNT - 2))
        BASE_COUNT=$((45 / REMAINING_FILES))
        REMAINDER=$((45 % REMAINING_FILES))

        for i in $(seq 1 $REMAINING_FILES); do
            COUNT=$BASE_COUNT
            if [ $i -eq $REMAINING_FILES ]; then
                COUNT=$((COUNT + REMAINDER))
            fi
            MOLECULE_COUNTS+=" $COUNT"
        done
    fi

    echo "Final Molecule Counts: $MOLECULE_COUNTS"

    # --------- Step 3: Run OPLS_box.py ---------
    OUTPUT_DATA="system.data"
    BOX_SIZE=50

    echo "Running OPLS_box.py with files: \"$LMP_FILES_TO_PASS\""
    #echo "Molecule counts: \"$MOLECULE_COUNTS\""
    python OPLS_box.py "$LMP_FILES_TO_PASS" -N "$MOLECULE_COUNTS" -O $OUTPUT_DATA -L $BOX_SIZE


    # Check outputs
    if [ ! -f "$OUTPUT_DATA" ] || [ ! -f "settings.lmp" ]; then
        echo "Run $RUN_ID failed during OPLS_box.py execution!"
        rm -f *.lmp  # Cleanup copied files
        cd ..
        continue
    fi

    # Move generated files back to project root
    mv $OUTPUT_DATA settings.lmp ../

    # Clean up copied .lmp files
    rm -f *.lmp
    cd ..

    # --------- Step 4: Generate run.in ---------
    python $GEN_DIR/gen_in.py $YAML_FILE
    if [ ! -f "run.in" ]; then
        echo "Run $RUN_ID failed during LAMMPS input generation!"
        continue
    fi

    # --------- Step 5: Move Inputs to Run Folder ---------
    mv system.data settings.lmp run.in $RUN_DIR/

    # --------- Step 6: Copy and Modify run_job.sh ---------
    cp ./gen_scripts/run_job.sh $RUN_DIR/

    cd $RUN_DIR
    qsub -N run_${RUN_ID}_lammps run_job.sh $RUN_ID
    cd ..
    # --------- Log Submission to run_summary.csv ---------
    TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
    MOLECULE_LIST=$(IFS=','; echo "${MOLECULE_NAMES[*]}")

    # Append to summary CSV
    echo "$RUN_ID,$MOLECULE_LIST,$MOLECULE_COUNTS,$TIMESTAMP" >> run_summary.csv

done

echo " "
echo "Batch processing completed."
rm run_list.csv
mv exec_list.csv run_list.csv  # Restore original CSV