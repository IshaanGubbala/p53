#!/bin/bash
# Render all p53 rescue animations with PyMOL

set -e  # Exit on error

PYMOL="/Applications/PyMOL.app/Contents/MacOS/PyMOL"
SCRIPT_DIR="reports/pymol_scripts"
OUTPUT_DIR="reports/figures/animations"

# Check PyMOL exists
if [ ! -f "$PYMOL" ]; then
    echo "ERROR: PyMOL not found at $PYMOL"
    echo "Please install PyMOL or update the PYMOL variable in this script"
    exit 1
fi

echo "========================================"
echo "Rendering p53 Rescue Animations"
echo "========================================"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Render animations
for target in R175H R248Q R273H; do
    echo "Rendering $target..."

    # Animation
    anim_script="$SCRIPT_DIR/animate_rescue_${target}.pml"
    if [ -f "$anim_script" ]; then
        echo "  → Animation (this will take ~5-10 minutes)..."
        "$PYMOL" -c "$anim_script" 2>&1 | grep -v "Warning:" || true
        echo "  ✓ Animation complete"
    fi

    # Comparison
    comp_script="$SCRIPT_DIR/compare_rescue_${target}.pml"
    if [ -f "$comp_script" ]; then
        echo "  → Side-by-side comparison..."
        "$PYMOL" -c "$comp_script" 2>&1 | grep -v "Warning:" || true
        echo "  ✓ Comparison complete"
    fi

    echo ""
done

echo "========================================"
echo "All animations rendered!"
echo "========================================"
echo ""
echo "Output files:"
ls -lh "$OUTPUT_DIR"/*.mp4 "$OUTPUT_DIR"/*.png 2>/dev/null || echo "  (files will appear in $OUTPUT_DIR)"
echo ""
echo "View results:"
echo "  open $OUTPUT_DIR"
