# Overlay: R175H - All 3 States Showing Structural Differences

set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set ambient, 0.5
set specular, 0.6
set depth_cue, 1
set bg_rgb, [1, 1, 1]
set cartoon_fancy_helices, 1

# Load all structures
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb, wt
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R175H/cancer_build/R175H_cancer.pdb, cancer
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R175H/rescue_build/R175H_rescued.pdb, rescued
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna

# Align all to wild-type
align cancer, wt
align rescued, wt
align dna, wt

# Offset structures slightly so all are visible
translate [-2, 0, 0], wt
# cancer stays at origin (0, 0, 0)
translate [2, 0, 0], rescued

hide everything

# Show DNA as transparent reference
show cartoon, dna
color gray80, dna
set cartoon_transparency, 0.7, dna

# Show all 3 proteins as transparent cartoons
show cartoon, wt
show cartoon, cancer
show cartoon, rescued
set cartoon_transparency, 0.8, wt
set cartoon_transparency, 0.8, cancer
set cartoon_transparency, 0.8, rescued
color cyan, wt
color red, cancer
color green, rescued

# Highlight mutation sites
select wt_site, wt and resi 175
select cancer_site, cancer and resi 175
select rescued_site, rescued and resi 175

show spheres, wt_site
show spheres, cancer_site
show spheres, rescued_site

color cyan, wt_site
color red, cancer_site
color yellow, rescued_site

set sphere_scale, 1.5, wt_site
set sphere_scale, 1.5, cancer_site
set sphere_scale, 1.5, rescued_site

# Show binding interface
select wt_interface, wt and resi 175+241+248+273+280
select cancer_interface, cancer and resi 175+241+248+273+280
select rescued_interface, rescued and resi 175+241+248+273+280

show sticks, wt_interface
show sticks, cancer_interface
show sticks, rescued_interface

set stick_radius, 0.25, wt_interface
set stick_radius, 0.25, cancer_interface
set stick_radius, 0.25, rescued_interface

set stick_transparency, 0.3, wt_interface
set stick_transparency, 0.3, cancer_interface
set stick_transparency, 0.3, rescued_interface

color cyan, wt_interface and elem C
color red, cancer_interface and elem C
color green, rescued_interface and elem C

# Orient
center wt_interface
orient wt_interface
zoom wt_interface, 5
turn y, 15
turn x, 10

# Render
viewport 1920, 1080
ray 1920, 1080
png /Users/ishaangubbala/Documents/p53/reports/figures/overlay_comparisons/R175H_overlay_comparison.png, dpi=300
print('Overlay image saved!')
