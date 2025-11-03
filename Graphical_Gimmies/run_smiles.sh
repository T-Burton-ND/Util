#!/bin/bash
# run_smiles.sh - Extract SMILES under "P" column and run mol_id.py
# All outputs go into mol_id.log

INPUT_FILE="reaction_list.txt"
LOG_FILE="mol_id.log"

# Clear log file at the start
: > "$LOG_FILE"

# Find the column index of "P" (1-based)
col_index=$(head -1 "$INPUT_FILE" | tr -s ' ' | tr '\t' ' ' | tr -s ' ' | \
            tr ' ' '\n' | nl -v1 | grep -w "P" | awk '{print $1}')

# Loop over SMILES in that column (skip header)
tail -n +2 "$INPUT_FILE" | awk -v col=$col_index '{print $col}' | while read -r smiles; do
    if [[ -n "$smiles" ]]; then
        echo ">>> Running mol_id.py on: $smiles" | tee -a "$LOG_FILE"
        python mol_id.py "$smiles" >> "$LOG_FILE" 2>&1
        echo "" >> "$LOG_FILE"
    fi
done
