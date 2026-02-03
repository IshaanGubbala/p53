# Quick p53 Rescue Animation (No Ray-Tracing for Speed)
# Target: R175H

# Settings
viewport 1280, 720
set antialias, 1
set bg_rgb, [1, 1, 1]
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# Load all three structures
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb, wt
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R175H/cancer_build/R175H_cancer.pdb, cancer
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R175H/rescue_build/R175H_rescued.pdb, rescued

# Align
align cancer, wt
align rescued, wt

# Hide all initially
hide everything

# Create 90-frame movie (3 seconds @ 30fps)
mset 1 x90

# === SCENE 1: Wild-Type (Frames 1-30) ===
frame 1
show cartoon, wt
color gray70, wt
select wt_site, wt and resi 175
show spheres, wt_site
color cyan, wt_site
set sphere_scale, 1.5, wt_site

# Create label
pseudoatom label1, pos=[30,30,0]
label label1, 'WILD-TYPE (Healthy)'
set label_size, 40
set label_color, black

# Orient
orient wt
zoom wt, 5

# Smooth rotation
mview store, 1
turn y, 60
mview store, 30

# === SCENE 2: Cancer Mutant (Frames 31-60) ===
frame 31

# Transition
hide cartoon, wt
hide spheres, wt_site
delete label1

# Show cancer mutant
show cartoon, cancer
color gray70, cancer
select cancer_site, cancer and resi 175
show spheres, cancer_site
color red, cancer_site
set sphere_scale, 1.8, cancer_site

# Show destabilized region
select destab, cancer and (cancer_site around 8)
show surface, destab
set transparency, 0.7, destab
color salmon, destab

# Label
pseudoatom label2, pos=[30,30,0]
label label2, 'CANCER MUTANT (R175H)'
set label_color, red

# Orient
orient cancer
zoom cancer, 5

# Smooth rotation
mview store, 31
turn y, 60
mview store, 60

# === SCENE 3: Rescued (Frames 61-90) ===
frame 61

# Transition
hide cartoon, cancer
hide spheres, cancer_site
hide surface, destab
delete label2

# Show rescued
show cartoon, rescued
color gray70, rescued
select rescued_site, rescued and resi 175
show spheres, rescued_site
color yellow, rescued_site
set sphere_scale, 1.6, rescued_site

# Show re-stabilized region
select stab, rescued and (rescued_site around 8)
show surface, stab
set transparency, 0.7, stab
color palegreen, stab

# Label
pseudoatom label3, pos=[30,30,0]
label label3, 'RESCUED (Stability Restored)'
set label_color, green

# Orient
orient rescued
zoom rescued, 5

# Smooth rotation
mview store, 61
turn y, 60
mview store, 90

# Interpolate all views
mview interpolate
mview reinterpolate

# === RENDER ===
# Render to PNG sequence (NO ray-tracing for speed)
set cache_frames, 0
mpng /Users/ishaangubbala/Documents/p53/reports/figures/animations/frames/quick_rescue_R175H_

print('='*60)
print('Animation frames saved!')
print('Create video with ffmpeg:')
print('cd /Users/ishaangubbala/Documents/p53/reports/figures/animations')
print('ffmpeg -framerate 30 -i frames/quick_rescue_R175H_%04d.png -vf scale=1280:720 -c:v libx264 -pix_fmt yuv420p quick_rescue_R175H.mp4')
print('='*60)