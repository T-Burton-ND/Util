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
