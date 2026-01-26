# PyMOL Animation: Before/After Rescue with DNA Binding
# Target: R248Q
# Rescue: M133L

# Set movie parameters
viewport 1920, 1080
set ray_trace_frames, 1
set ray_trace_mode, 1
set antialias, 2
set ambient, 0.5
set specular, 0.6
set ray_shadows, 1
set depth_cue, 1
set bg_rgb, [1, 1, 1]

# Load structures
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb, wildtype
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R248Q/cancer_build/R248Q_cancer.pdb, cancer_mutant
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R248Q/rescue_build/R248Q_rescued.pdb, rescued
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna

# Align all structures to wild-type
align cancer_mutant, wildtype
align rescued, wildtype
align dna, wildtype

# Initially hide all
hide everything

# =====================
# SCENE 1: Wild-Type (Frames 1-60)
# =====================
mset 1 x180
frame 1

# Show wild-type only
show cartoon, wildtype
show cartoon, dna
color gray70, wildtype
color blue, dna
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 0.5

# DNA-binding interface
select dna_contact, wildtype and (resi 248 or resi 241 or resi 248 or resi 273 or resi 280)
show sticks, dna_contact
color cyan, dna_contact
set stick_radius, 0.2, dna_contact

# Show hydrogen bonds between DNA and protein
distance hbonds_wt, (dna_contact), (dna), 3.5
hide labels, hbonds_wt
color yellow, hbonds_wt
set dash_width, 3, hbonds_wt
set dash_gap, 0.2, hbonds_wt

# Add text label
set label_size, 40
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.3
pseudoatom label_wt, pos=[50,50,50]
label label_wt, 'WILD-TYPE: Stable DNA Binding'

# Orient view
orient dna_contact
zoom dna_contact, 8

# Gentle rotation for Scene 1
mview store, 1
turn y, 15
mview store, 30
turn y, 15
mview store, 60
mview interpolate

# =====================
# SCENE 2: Cancer Mutant (Frames 61-120)
# =====================
frame 61

# Transition: fade out wild-type, fade in cancer mutant
hide cartoon, wildtype
hide sticks, dna_contact
delete hbonds_wt
delete label_wt

# Show cancer mutant
show cartoon, cancer_mutant
color gray70, cancer_mutant

# Highlight mutation site in RED
select cancer_site, cancer_mutant and resi 248
show spheres, cancer_site
color red, cancer_site
set sphere_scale, 1.8, cancer_site

# Show disrupted DNA binding interface
select mutant_interface, cancer_mutant and (resi 248 or resi 241 or resi 248 or resi 273 or resi 280)
show sticks, mutant_interface
color orange, mutant_interface
set stick_transparency, 0.4, mutant_interface

# Show BROKEN hydrogen bonds (fewer/weaker)
distance hbonds_mutant, (mutant_interface), (dna), 3.5
hide labels, hbonds_mutant
color red, hbonds_mutant
set dash_width, 2, hbonds_mutant
set dash_gap, 0.4, hbonds_mutant

# Structural distortion visualization
select distorted_region, cancer_mutant and (cancer_site around 8)
show surface, distorted_region
set transparency, 0.7, distorted_region
color red, distorted_region

# Label
pseudoatom label_mutant, pos=[50,50,50]
label label_mutant, 'CANCER MUTANT (R248Q): Destabilized, Poor DNA Binding'

# Re-orient
orient cancer_site
zoom cancer_site, 8

# Wobble motion to show instability
mview store, 61
turn y, -10
turn x, 5
mview store, 75
turn y, 20
turn x, -10
mview store, 90
turn y, -10
turn x, 5
mview store, 105
turn y, 5
mview store, 120
mview interpolate

# =====================
# SCENE 3: Rescued Structure (Frames 121-180)
# =====================
frame 121

# Transition: fade out cancer mutant, fade in rescued
hide cartoon, cancer_mutant
hide spheres, cancer_site
hide sticks, mutant_interface
hide surface, distorted_region
delete hbonds_mutant
delete label_mutant

# Show rescued structure
show cartoon, rescued
color gray70, rescued

# Highlight rescue mutations in GREEN
select rescue_1, rescued and resi 133
show spheres, rescue_1
color green, rescue_1
set sphere_scale, 1.5, rescue_1
label rescue_1 and name CA, 'M133L'
select rescued_cancer_site, rescued and resi 248
show spheres, rescued_cancer_site
color salmon, rescued_cancer_site
set sphere_scale, 1.6, rescued_cancer_site
label rescued_cancer_site and name CA, 'R248Q'

# Show RESTORED DNA binding interface
select rescued_interface, rescued and (resi 248 or resi 241 or resi 248 or resi 273 or resi 280)
show sticks, rescued_interface
color cyan, rescued_interface
set stick_radius, 0.2, rescued_interface

# Show RESTORED hydrogen bonds
distance hbonds_rescued, (rescued_interface), (dna), 3.5
hide labels, hbonds_rescued
color lime, hbonds_rescued
set dash_width, 3, hbonds_rescued
set dash_gap, 0.2, hbonds_rescued

# Highlight stabilized region (green surface)
select stabilized_region, rescued and (rescued_cancer_site around 8)
show surface, stabilized_region
set transparency, 0.7, stabilized_region
color green, stabilized_region

# Label
pseudoatom label_rescued, pos=[50,50,50]
label label_rescued, 'RESCUED: Stability Restored, DNA Binding Recovered'

# Orient view
orient rescued_interface
zoom rescued_interface, 8

# Smooth rotation for Scene 3
mview store, 121
turn y, 30
mview store, 150
turn y, 30
mview store, 180
mview interpolate

# =====================
# RENDER ANIMATION
# =====================

# Set output format and render settings
set movie_fps, 30
set ray_trace_frames, 1
set cache_frames, 0

print('Animation rendering started (this takes ~10 min)...')
# Use mpng to render individual PNG frames, then ffmpeg will combine
mpng /Users/ishaangubbala/Documents/p53/reports/figures/animations/frames/rescue_animation_R248Q_

print('Frames rendered! Use ffmpeg to create video:')
print('cd /Users/ishaangubbala/Documents/p53/reports/figures/animations')
print('ffmpeg -framerate 30 -i frames/rescue_animation_R248Q_%04d.png -c:v libx264 -pix_fmt yuv420p rescue_animation_R248Q.mp4')
