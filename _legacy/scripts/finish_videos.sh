#!/bin/bash
# Create MP4 videos from all completed frame sequences

cd /Users/ishaangubbala/Documents/p53/reports/figures/animations

for target in R175H R248Q R273H; do
    for type in cancer rescued; do
        vid="${target}_${type}"

        # Check if frames exist
        frame_count=$(ls frames/${vid}_binding_*.png 2>/dev/null | wc -l)

        if [ "$frame_count" -eq 60 ]; then
            if [ ! -f "${vid}_binding.mp4" ]; then
                echo "Creating ${vid}_binding.mp4..."
                ffmpeg -y -framerate 30 -i frames/${vid}_binding_%04d.png \
                    -c:v libx264 -pix_fmt yuv420p ${vid}_binding.mp4 2>&1 | tail -1
            fi
        fi
    done
done

echo ""
echo "Videos created:"
ls -lh *_binding.mp4 2>/dev/null

echo ""
echo "Cleaning up frame files..."
rm -f frames/*_binding_*.png

echo "Done!"
