#!/usr/bin/env bash
# <file_name>.sh — <one-line summary>
# -----------------------------------------------------------------
# Short description of what the script does.
#
# Usage:
#   ./<file_name>.sh [OPTIONS]
#
# Options:
#   -i, --input FILE     Input file
#   -o, --output DIR     Output directory (default: ./out)
#   -h, --help           Show this help and exit
#
# Examples:
#   ./<file_name>.sh -i input.mov -o ./gifs
#

set -euo pipefail

usage() {
  sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'
}

# Default values
out="./out"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)  in="$2"; shift 2;;
    -o|--output) out="$2"; shift 2;;
    -h|--help)   usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 1;;
  esac
done

# --- your code goes here ---
