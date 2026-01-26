# PyMOL script for multi-view rendering
# Target: R175H, Tier: 1mut

# Load structure
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb
remove solvent

# Set background and basic display
bg_color white
hide everything
show cartoon
color gray80, all
set cartoon_smooth_loops, 0
set cartoon_fancy_helices, 1
set ray_shadows, 1
set antialias, 2
set ambient, 0.4
set specular, 0.5

# Highlight target mutation
select target_mut, resi 175
show spheres, target_mut
color red, target_mut
set sphere_scale, 1.5, target_mut

select rescue_1, resi 133
show spheres, rescue_1
color green, rescue_1
set sphere_scale, 1.3, rescue_1
select mutation_sites, resi 175+133
select nearby, (mutation_sites around 10) and not mutation_sites
show sticks, nearby
color yellow, nearby
set stick_radius, 0.15, nearby
set stick_transparency, 0.5, nearby

# Orient
orient
zoom all, 5

viewport 1200, 1000

# Front View
ray 1200, 1000
png /Users/ishaangubbala/Documents/p53/reports/figures/pymol_renders/protein_pymol_R175H_1mut_front.png, dpi=300

orient
zoom all, 5

# Side View
turn y, 90
ray 1200, 1000
png /Users/ishaangubbala/Documents/p53/reports/figures/pymol_renders/protein_pymol_R175H_1mut_side.png, dpi=300

orient
zoom all, 5

# Top View
turn x, 90
ray 1200, 1000
png /Users/ishaangubbala/Documents/p53/reports/figures/pymol_renders/protein_pymol_R175H_1mut_top.png, dpi=300

orient
zoom all, 5

# Back View
turn y, 180
ray 1200, 1000
png /Users/ishaangubbala/Documents/p53/reports/figures/pymol_renders/protein_pymol_R175H_1mut_back.png, dpi=300

print('Multi-view render complete!')