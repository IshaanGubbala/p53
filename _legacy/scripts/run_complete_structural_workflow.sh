#!/bin/bash
# Complete workflow to create REAL structural comparisons showing conformational changes

set -e  # Exit on error

echo "========================================================================"
echo "p53 Rescue Structural Modeling & Visualization - Complete Workflow"
echo "========================================================================"
echo ""

# Check prerequisites
echo "Checking prerequisites..."
echo ""

# Check EvoEF2
if ! command -v EvoEF2 &> /dev/null; then
    echo "✗ EvoEF2 not found!"
    echo "  Install from: https://github.com/tommyhuangthu/EvoEF2"
    echo "  Or add to PATH: export PATH=/path/to/EvoEF2/build:\$PATH"
    exit 1
else
    echo "✓ EvoEF2 found: $(which EvoEF2)"
fi

# Check PyMOL
if [ ! -f "/Applications/PyMOL.app/Contents/MacOS/PyMOL" ]; then
    echo "✗ PyMOL not found at /Applications/PyMOL.app/"
    echo "  Install from: https://pymol.org/"
    exit 1
else
    echo "✓ PyMOL found"
fi

# Check OpenMM (optional)
if python -c "import openmm" 2>/dev/null; then
    echo "✓ OpenMM found (MD relaxation will be performed)"
    USE_MD=true
else
    echo "⚠ OpenMM not found (MD relaxation will be skipped)"
    echo "  Install with: conda install -c conda-forge openmm"
    USE_MD=false
fi

echo ""
echo "========================================================================"
echo "STEP 1: Model Mutant Structures"
echo "========================================================================"
echo ""
echo "This will create actual mutant structures with real conformational changes."
echo "Expected time: 15-30 minutes"
echo ""

read -p "Continue with structure modeling? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 0
fi

python scripts/model_rescue_structures.py

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Structure modeling failed!"
    echo "  Check error messages above"
    exit 1
fi

echo ""
echo "✓ Structure modeling complete!"
echo ""

# Check if structures were created
STRUCT_DIR="Data/processed/modeled_structures"
if [ ! -d "$STRUCT_DIR/R175H" ]; then
    echo "✗ No modeled structures found!"
    exit 1
fi

echo "Modeled structures:"
ls -lh "$STRUCT_DIR"/R175H/*.pdb | awk '{print "  " $9 " (" $5 ")"}'
echo ""

echo "========================================================================"
echo "STEP 2: Generate Comparison Scripts"
echo "========================================================================"
echo ""

python scripts/create_structural_comparison.py

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Script generation failed!"
    exit 1
fi

echo ""
echo "✓ Comparison scripts generated!"
echo ""

echo "========================================================================"
echo "STEP 3: Render High-Quality Images"
echo "========================================================================"
echo ""
echo "This will render side-by-side comparisons showing REAL structural changes."
echo "Expected time: 5-10 minutes"
echo ""

read -p "Continue with rendering? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Skipping rendering. You can run manually:"
    echo "  for script in reports/pymol_scripts/compare_structures_*.pml; do"
    echo "      /Applications/PyMOL.app/Contents/MacOS/PyMOL -c \"\$script\""
    echo "  done"
    exit 0
fi

PYMOL="/Applications/PyMOL.app/Contents/MacOS/PyMOL"
SCRIPT_DIR="reports/pymol_scripts"

for script in "$SCRIPT_DIR"/compare_structures_*.pml; do
    if [ -f "$script" ]; then
        target=$(basename "$script" .pml | sed 's/compare_structures_//')
        echo "Rendering $target..."
        "$PYMOL" -c "$script" 2>&1 | grep -E "(RMSD|Image saved|Error)" || true
        echo ""
    fi
done

echo ""
echo "========================================================================"
echo "COMPLETE! Structural Comparisons Generated"
echo "========================================================================"
echo ""

OUTPUT_DIR="reports/figures/structural_comparisons"
if [ -d "$OUTPUT_DIR" ]; then
    echo "Output images:"
    ls -lh "$OUTPUT_DIR"/*.png | awk '{print "  " $9 " (" $5 ")"}'
    echo ""

    echo "View results:"
    echo "  open $OUTPUT_DIR"
    echo ""

    # Count successfully generated images
    IMAGE_COUNT=$(ls "$OUTPUT_DIR"/*.png 2>/dev/null | wc -l)
    echo "Generated $IMAGE_COUNT high-quality comparison images"
    echo ""
fi

echo "Next steps for science fair:"
echo "  1. Review images to see real structural changes"
echo "  2. Check RMSD values (printed above) - should be 1-2 Å"
echo "  3. Use images in poster/presentation"
echo ""

if [ "$USE_MD" = false ]; then
    echo "TIP: Install OpenMM for more dramatic conformational changes:"
    echo "  conda install -c conda-forge openmm"
    echo "  Then re-run this script"
    echo ""
fi

echo "Documentation: see REAL_STRUCTURE_WORKFLOW.md"
echo ""
echo "✓ All done!"
