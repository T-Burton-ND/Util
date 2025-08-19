#!/bin/bash
# gif_from_mov.sh — Convert .mov to .gif with tunable settings
# ------------------------------------------------------------
# Short description:
#   Wraps ffmpeg to produce an animated GIF from a .mov file, with configurable
#   duration, output width, FPS, and loop count. Creates the output folder if
#   needed and preserves aspect ratio (even dimensions via pad).
#
# Usage:
#   ./gif_from_mov.sh path/to/input.mov path/to/output_folder
#
# Options (edit variables at top of file):
#   DURATION      Length of GIF in seconds            (default: 10)
#   SCALE_WIDTH   Output width in pixels              (default: 480)
#   FPS           Frames per second                   (default: 15)
#   LOOP          0=infinite, 1=once, etc.           (default: 0)
#
# Dependencies:
#   ffmpeg
#
# Examples:
#   ./gif_from_mov.sh video.mov ./gifs
#   DURATION=6 FPS=12 SCALE_WIDTH=600 ./gif_from_mov.sh input.mov outdir
# =================== User Configurable Parameters ===================
DURATION=10        # Length of GIF in seconds
SCALE_WIDTH=480    # Output width in pixels (height auto-adjusts)
FPS=15             # Frames per second for the GIF
LOOP=0             # 0 = infinite loop, 1 = play once, etc.
# ===================================================================

# Usage: gif_from_mov.sh path/to/input.mov path/to/output_folder#

# Check for required arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 path/to/input.mov path/to/output_folder"
    echo "Example: ./gif_from_mov.sh input.mov ./gifs"
    exit 1
fi

INPUT_MOV="$1"
OUTPUT_DIR="$2"

# Extract filename without extension
BASENAME=$(basename "$INPUT_MOV" .mov)
OUTPUT_GIF="${OUTPUT_DIR}/${BASENAME}.gif"

# Create output folder if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Build and run ffmpeg command
ffmpeg -t "$DURATION" -i "$INPUT_MOV" \
  -vf "fps=${FPS},scale=${SCALE_WIDTH}:-1:flags=lanczos,pad=ceil(iw/2)*2:ceil(ih/2)*2" \
  -loop "$LOOP" "$OUTPUT_GIF"

echo "✅ GIF created: $OUTPUT_GIF"
