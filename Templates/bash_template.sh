#!/usr/bin/env bash
# <file_name>.sh — <one‑line summary>
# ---------------------------------------------------------------
# Short description (1–2 lines).
#
# Usage:
#   ./<file_name>.sh [OPTIONS] <args>
#
# Options:
#   -i, --input FILE     Input file
#   -o, --output DIR     Output directory (default: ./out)
#   -h, --help           Show this help and exit
#
# Dependencies:
#   ffmpeg, awk  # list external commands actually required
#
# Examples:
#   ./<file_name>.sh input.mov ./gifs
#
set -euo pipefail

usage() {
  sed -n '1,50p' "$0" | sed 's/^# \{0,1\}//' | sed 's/^$//'
}

# default values
out="./out"

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   in="$2"; shift 2;;
    -o|--output)  out="$2"; shift 2;;
    -h|--help)    usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 2;;
  esac
done

# ... your script ...
