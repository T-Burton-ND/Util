#!/usr/bin/env bash
set -euo pipefail
# bulk_mol_id.sh — Batch molecule rendering + metadata rename
# -----------------------------------------------------------
# Runs mol_id.py on a list of SMILES strings and, when available, renames the
# generated image file using the IUPAC name. Handles temporary logs, unique
# filenames, and safe sanitization of output names.
#
# Usage:
#   PYTHON_BIN=python3 bash bulk_mol_id.sh
#
# Environment variables:
#   PYTHON_BIN   Python executable to use (default: python)
#
# Dependencies:
#   python (with rdkit + requests installed), awk, sed, tee, mktemp
#
# Examples:
#   bash bulk_mol_id.sh
#   PYTHON_BIN=/users/tburton2/.conda/envs/util_env/bin/python bash bulk_mol_id.sh
# end help
# ---- Config ----
PYTHON_BIN="${PYTHON_BIN:-python}"   # override: PYTHON_BIN=python3 bash bulk_mol_id.sh

__VERSION__="0.1.0"

# ---- Help ----
usage() {
  # Print header comment block up to the Config marker
  sed -n '/^# bulk_mol_id.sh/,/^# end help/p' "$0" | sed -E 's/^# ?//'
  exit 0
}

# Handle -h/--help/--version without changing core behavior
for arg in "$@"; do
  case "$arg" in
    -h|--help) usage ;;
    --version) echo "bulk_mol_id.sh ${__VERSION__}"; exit 0 ;;
  esac
done

# ---- Inputs ----
smiles=(
reaction                                 R                                                            P                                                            type            barrier   
XSASRUDTFFBDDK_0_0                       OOCCC=O                                                      CC(=O)COO                                                    P_unintended     71.594168
XSASRUDTFFBDDK_1_0                       OOCCC=O                                                      OCC=O.C=O                                                    intended         78.262428
XSASRUDTFFBDDK_2_0                       C=C=O.C=O.O                                                  O=CCC=O.O                                                    unintended       50.893881
XSASRUDTFFBDDK_3_0                       OOCCC=O                                                      OOC=O.C=C                                                    intended         51.199809
XSASRUDTFFBDDK_4_0                       OOCCC=O                                                      OC[CH]C=O.[OH]                                               P_unintended     55.967973
XSASRUDTFFBDDK_6_0                       OO/C=C/C=O.[H][H]                                            OO/C=C\C=O.[H][H]                                            R_unintended     49.619981
XSASRUDTFFBDDK_7_0                       OOCCC=O                                                      O=CCC=O.O                                                    intended         49.746654
XSASRUDTFFBDDK_9_0                       OOCCC=O                                                      OOCOC=C                                                      intended         83.164436
XSASRUDTFFBDDK_10_0                      OOCCC=O                                                      CC(=O)COO                                                    intended         71.594168
XSASRUDTFFBDDK_11_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.249044
XSASRUDTFFBDDK_12_0                      OOCCC=O                                                      OCOCC=O                                                      intended         81.383843
XSASRUDTFFBDDK_13_0                      OOCCC=O                                                      OCOOC=C                                                      intended         78.795411
XSASRUDTFFBDDK_14_0                      OOCCC=O                                                      OOCCC=O                                                      intended          5.002390
XSASRUDTFFBDDK_15_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.251434
XSASRUDTFFBDDK_16_0                      OOC=C.C=O                                                    COOOC=C                                                      R_unintended     52.108031
XSASRUDTFFBDDK_17_0                      OOCCC=O                                                      [O-]C=CC.O=[OH+]                                             P_unintended     84.760994
XSASRUDTFFBDDK_18_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.249044
XSASRUDTFFBDDK_19_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.251434
XSASRUDTFFBDDK_20_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.251434
XSASRUDTFFBDDK_21_0                      OOCCC=O                                                      OOC/C=C/O                                                    intended         63.747610
XSASRUDTFFBDDK_22_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.251434
XSASRUDTFFBDDK_23_0                      OOCCC=O                                                      OOC=C.C=O                                                    P_unintended     55.100382
XSASRUDTFFBDDK_24_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.249044
XSASRUDTFFBDDK_25_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.249044
XSASRUDTFFBDDK_26_0                      OOCCC=O                                                      CCOOC=O                                                      intended         70.131453
XSASRUDTFFBDDK_27_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.251434
XSASRUDTFFBDDK_28_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.251434
XSASRUDTFFBDDK_29_0                      OOCCC=O                                                      [O-]/C=C\C[O+]=O.[H][H]                                      P_unintended    110.408700
XSASRUDTFFBDDK_30_0                      OOCCC=O                                                      C=C=O.C=O.O                                                  P_unintended     61.489006
XSASRUDTFFBDDK_31_0                      OCCC(=O)O                                                    C1COOC1=O.[H][H]                                             R_unintended    122.934990
XSASRUDTFFBDDK_32_0                      OOCOC=C                                                      OOCOC=C                                                      unintended        5.148184
XSASRUDTFFBDDK_33_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.249044
XSASRUDTFFBDDK_34_0                      OOCCC=O                                                      C1COOOC1                                                     intended         88.083174
XSASRUDTFFBDDK_35_0                      OOCCC=O                                                      OC=O.CC=O                                                    P_unintended     32.251434
XSASRUDTFFBDDK_36_0                      OOCCC=O                                                      OC=O.CC=O
OC=O.CC=O                                                    P_unintended     32.249044
[O-][O+]=C.OC=C                                              unintended        0.468451
OC1OCCO1                                                     R_unintended     54.839866

)

# ---- Helpers ----
sanitize() {
  # Keep only letters/numbers; turn others into underscores; collapse repeats; trim ends; cap length.
  local s="$1"
  s="${s//[^a-zA-Z0-9]/_}"
  s="$(printf '%s' "$s" | sed -E 's/_+/_/g; s/^_+//; s/_+$//')"
  printf '%s' "${s:0:200}"
}

unique_path_with_suffix() {
  local base="$1" ext="$2"
  local candidate="${base}${ext}"
  local i=1
  while [[ -e "$candidate" ]]; do
    candidate="${base}_$i${ext}"
    ((i++))
  done
  printf '%s' "$candidate"
}

process_one() {
  local smi="$1"

  local tmpdir
  tmpdir="$(mktemp -d)"
  local log="$tmpdir/mol_id.out"

  # Portable: capture both stdout+stderr
  if ! "$PYTHON_BIN" mol_id.py "$smi" 2>&1 | tee "$log" >/dev/null; then
    echo "[!] mol_id.py failed for: $smi" >&2
    rm -rf "$tmpdir"
    return 1
  fi

  # Extract image path
  local img_path
  img_path="$(awk -F': ' '/^\[✓] Image saved to:/{print $2; exit}' "$log")"

  if [[ -z "${img_path:-}" || ! -e "$img_path" ]]; then
    echo "[!] Could not find created image path for: $smi" >&2
    rm -rf "$tmpdir"
    return 1
  fi

  # Extract IUPAC name (if present)
  local iupac
  iupac="$(awk -F': ' '/^\[✓] IUPAC name:/{print $2; exit}' "$log" || true)"

  if [[ -z "${iupac:-}" ]]; then
    echo "[i] No IUPAC name found; keeping original file: $img_path"
    rm -rf "$tmpdir"
    return 0
  fi

  local base dir ext dst
  base="$(sanitize "$iupac")"
  if [[ -z "$base" ]]; then
    echo "[i] IUPAC sanitized to empty; keeping original file: $img_path"
    rm -rf "$tmpdir"
    return 0
  fi

  dir="$(dirname "$img_path")"
  ext=".png"
  dst="$(unique_path_with_suffix "$dir/$base" "$ext")"

  mv -- "$img_path" "$dst"
  echo "[✓] Renamed: $img_path -> $dst"

  rm -rf "$tmpdir"
}

# ---- Main ----
for smi in "${smiles[@]}"; do
  process_one "$smi"
done
