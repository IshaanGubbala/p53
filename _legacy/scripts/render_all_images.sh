#!/bin/bash
# Render all DNA binding visualization images

cd /Users/ishaangubbala/Documents/p53

# Create output directories
mkdir -p reports/figures/dna_binding_closeups
mkdir -p reports/figures/overlay_comparisons

echo "================================================================================"
echo "Rendering DNA Binding Visualizations"
echo "================================================================================"
echo ""

# Render overlay comparisons (showing all 3 structures aligned)
echo "1. Rendering overlay comparisons (all 3 states superimposed)..."
for target in R175H R248Q R273H; do
    echo "   - ${target} overlay..."
    /Applications/PyMOL.app/Contents/MacOS/PyMOL -c reports/pymol_scripts/overlay_${target}.pml 2>&1 | grep "saved" || echo "     (rendering...)"
done

echo ""
echo "2. Rendering individual close-ups (optional - remove labels first)..."
echo "   Skipping for now - scripts need label removal"

echo ""
echo "================================================================================"
echo "Done! Check outputs:"
echo "================================================================================"
echo ""
echo "Overlay comparisons: reports/figures/overlay_comparisons/"
ls -lh reports/figures/overlay_comparisons/*.png 2>/dev/null || echo "  (not yet rendered)"
echo ""
echo "Close-ups: reports/figures/dna_binding_closeups/"
ls -lh reports/figures/dna_binding_closeups/*.png 2>/dev/null || echo "  (not yet rendered)"
