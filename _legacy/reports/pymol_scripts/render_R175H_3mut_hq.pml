# PyMOL script for rendering p53 rescue mutations
# Target: R175H, Tier: 3mut
# ΔΔG Gain: -13.9, Risk: 0.000

# Load structure
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb
remove solvent

# Set background and basic display
bg_color white
hide everything
show cartoon
color gray80, all
set cartoon_smooth_loops, 0
set cartoon_flat_sheets, 0
set cartoon_fancy_helices, 1
set ray_shadows, 1
set ray_trace_mode, 1
set antialias, 2
set ambient, 0.4
set specular, 0.5
set shininess, 10
set depth_cue, 1
set ray_trace_fog, 1

# Highlight target mutation in RED
select target_mut, resi 175
show spheres, target_mut
color red, target_mut
set sphere_scale, 1.5, target_mut
set sphere_transparency, 0.0, target_mut

# Rescue mutation 1: A189S
select rescue_1, resi 189
show spheres, rescue_1
color green, rescue_1
set sphere_scale, 1.3, rescue_1
set sphere_transparency, 0.0, rescue_1

# Rescue mutation 2: M133L
select rescue_2, resi 133
show spheres, rescue_2
color green, rescue_2
set sphere_scale, 1.3, rescue_2
set sphere_transparency, 0.0, rescue_2

# Rescue mutation 3: Y163F
select rescue_3, resi 163
show spheres, rescue_3
color green, rescue_3
set sphere_scale, 1.3, rescue_3
set sphere_transparency, 0.0, rescue_3

# Show nearby residues (within 10Å)
select mutation_sites, resi 175+189+133+163
select nearby, (mutation_sites around 10) and not mutation_sites
show sticks, nearby
color yellow, nearby
set stick_radius, 0.15, nearby
set stick_transparency, 0.5, nearby

# Set view
orient
zoom all, 5

# Add labels
label target_mut and name CA, 'R175H'
label rescue_1 and name CA, 'A189S'
label rescue_2 and name CA, 'M133L'
label rescue_3 and name CA, 'Y163F'

set label_size, 20
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.3

# Ray trace and save
viewport 2400, 2000
ray 2400, 2000
png /Users/ishaangubbala/Documents/p53/reports/figures/pymol_renders/protein_pymol_R175H_3mut_hq.png, dpi=300

# Deselect all
deselect

print('Render complete!')