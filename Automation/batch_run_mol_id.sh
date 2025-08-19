#!/bin/bash

# ------------------------------------------------------------------------------
# Batch Metadata Extraction Script for mol_id.py
# Author : Thomas J. Burton — Savoie Research Group
# Usage  : ./batch_run_mol_id.sh
# ------------------------------------------------------------------------------

# batch_run_mol_id.sh — Batch metadata extraction via mol_id.py
# -------------------------------------------------------------
# Feeds a predefined list of molecules (SMILES or InChIKeys) to Visualization_Scripts/mol_id.py,
# captures RDKit canonical SMILES / InChI / InChIKey / IUPAC / common names, and writes a CSV.
#
# Usage:
#   ./batch_run_mol_id.sh
#
# Output:
#   mol_metadata.csv  (columns: input_string, rdkit_smiles, inchi, inchikey, iupac_name, common_names, n_carbons)
#
# Dependencies:
#   python with rdkit+requests, grep, sed, wc, tr
#
# Notes:
#   Adjust the 'inputs' array below as needed. The script calls ../Visualization_Scripts/mol_id.py.
# Define input molecules (SMILES or InChIKey strings)
inputs=(
  'O=C1OCCO1'
  'CCC1OC(=O)OC1C'
  ' O=c1occo1'
  'CC1OC(=O)OC1(C)C'
  'CCCC1COC(=O)O1'
  'CCCC1OC(=O)OC1C'
  "CC(C)C1OC(=O)OC1C"
  "CCC1OC(=O)OC1CC"
  "CCC1OC(=O)OC1(C)C"
  "CCC1(C)OC(=O)OC1C"
  "O=c1occo1"
  "CCc1oc(=O)oc1C"
  "Cc1oc(=O)oc1C"
  "CCc1oc(=O)oc1CC"
  "Cc1oc(=O)oc1C(C)C"
  "CCCc1coc(=O)o1"
  "CCCc1oc(=O)oc1C"
  "COC(=O)OC"
  "CCOC(=O)OCC"
  "CCOC(=O)OC"
  "CCOC(=O)OC(C)C"
  "CCCOC(=O)OCC"
  "CC(C)OC(=O)OC(C)C"
  "CCCOC(=O)OC(C)C"
  "CCCOC(=O)OCCC"
  "CCOC(C)=O"
  "CC(=O)OC(C)C"
  "CC(=O)OCC(C)C"
  "CCC(C)COC(C)=O"
  "CCC(=O)OC(C)C"
  "CCCOC(=O)CC"
  "CC(C)OC(=O)C(C)C"
)

# Output CSV
output_csv="mol_metadata.csv"
echo "input_string,rdkit_smiles,inchi,inchikey,iupac_name,common_names,n_carbons" > "$output_csv"

# Loop through each input
for input_str in "${inputs[@]}"; do
  echo "[•] Processing: $input_str"

  # Run mol_id.py and capture output
  result=$(python ../Visualization_Scripts/mol_id.py "$input_str" 2>&1)

  # Check for errors
  if echo "$result" | grep -q "\[✗\] Error"; then
    echo "[!] Error processing: $input_str"
    echo "$input_str,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR" >> "$output_csv"
    continue
  fi

  # Extract fields
  rdkit_smiles=$(echo "$result" | grep -E 'SMILES' | head -n 1 | sed 's/.*SMILES[^:]*: //')
  inchi=$(echo "$result" | grep -E 'InChI:' | sed 's/.*InChI: //')
  inchikey=$(echo "$result" | grep -E 'InChIKey:' | sed 's/.*InChIKey: //')
  iupac=$(echo "$result" | grep -E 'IUPAC name:' | sed 's/.*IUPAC name: //')
  common=$(echo "$result" | grep -E 'Common names:' | sed 's/.*Common names: //; s/, /;/g')

  # --- Modular Carbon Count ---
  if [[ "$rdkit_smiles" != "" && "$rdkit_smiles" != "N/A" ]]; then
    n_carbons=$(echo "$rdkit_smiles" | grep -o "C" | wc -l | tr -d ' ')
  else
    n_carbons="N/A"
  fi

  # Fallbacks
  rdkit_smiles="${rdkit_smiles:-N/A}"
  inchi="${inchi:-N/A}"
  inchikey="${inchikey:-N/A}"
  iupac="${iupac:-N/A}"
  common="${common:-N/A}"

  # Write row to CSV
  echo "\"$input_str\",\"$rdkit_smiles\",\"$inchi\",\"$inchikey\",\"$iupac\",\"$common\",\"$n_carbons\"" >> "$output_csv"
done

echo "[✓] All done. Results saved to: $output_csv"
