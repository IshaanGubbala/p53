#!/bin/bash
# Render all DNA binding close-up images

# Create output directory
cd /Users/ishaangubbala/Documents/p53
test -d reports/figures/dna_binding_closeups || {
    mkdir -p reports/figures/dna_binding_closeups
}

echo "Rendering DNA binding close-up images..."
echo ""

# Render all closeup scripts
for script in reports/pymol_scripts/closeup_*.pml; do
    if [ -f "$script" ]; then
        name=$(basename "$script" .pml)
        echo "Rendering $name..."
        /Applications/PyMOL.app/Contents/MacOS/PyMOL -c "$script" 2>&1 | grep "Image saved"
    fi
done

echo ""
echo "Done! Images saved to: reports/figures/dna_binding_closeups/"
ls -lh reports/figures/dna_binding_closeups/*.png 2>/dev/null || echo "No images found"
