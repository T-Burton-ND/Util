#!/usr/bin/env bash
set -euo pipefail

# Usage: ./recursive_unzip.sh [target_dir] [--keep]
# Default target_dir is current directory. Use --keep to keep .gz files.

TARGET_DIR="${1:-.}"
KEEP=false
if [[ "${2:-}" == "--keep" ]]; then
  KEEP=true
fi

# Extract tarballs first (tar.gz / .tgz)
find "$TARGET_DIR" -type f \( -name '*.tar.gz' -o -name '*.tgz' \) -print0 |
while IFS= read -r -d '' file; do
  dir=$(dirname "$file")
  echo "Extracting tar archive: $file -> $dir"
  tar -xzf "$file" -C "$dir"
  if ! $KEEP; then
    echo "Removing archive: $file"
    rm -f "$file"
  fi
done

# Decompress plain .gz files (but skip any .tar.gz already handled)
find "$TARGET_DIR" -type f -name '*.gz' ! -name '*.tar.gz' -print0 |
while IFS= read -r -d '' file; do
  out="${file%.gz}"
  echo "Decompressing: $file -> $out"
  if $KEEP; then
    # keep original: write decompressed output alongside original
    if ! gzip -dc -- "$file" > "$out"; then
      echo "Warning: Skipping $file (not in gzip format)" >&2
      rm -f "$out" # remove partial output if created
    fi
  else
    # default: replace .gz with decompressed file
    if ! gzip -d -- "$file"; then
      echo "Warning: Skipping $file (not in gzip format)" >&2
    fi
  fi
done