#!/bin/bash

INPUT_FILE="smiles_list.txt"
OUTPUT_DIR="images"
USE_COLOR=true  # Set to false for B&W

mkdir -p "$OUTPUT_DIR"

while IFS= read -r smiles || [[ -n "$smiles" ]]; do
    # Generate safe filename (hash or index if needed)
    filename=$(echo "$smiles" | sed 's/[^a-zA-Z0-9]/_/g')
    output_path="${OUTPUT_DIR}/${filename}.png"
    
    if [ "$USE_COLOR" = true ]; then
        ./rdkit_molecule_viz.py "$smiles" -o "$output_path" --multicolor
    else
        ./rdkit_molecule_viz.py "$smiles" -o "$output_path"
    fi
done < "$INPUT_FILE"
