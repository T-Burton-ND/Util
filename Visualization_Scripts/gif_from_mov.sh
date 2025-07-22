#!/bin/bash

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
