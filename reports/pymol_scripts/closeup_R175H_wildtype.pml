# Close-up: R175H Wildtype DNA Binding

# High-quality settings
set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set ambient, 0.5
set specular, 0.6
set depth_cue, 1
set bg_rgb, [1, 1, 1]
set cartoon_fancy_helices, 1

# Load structures
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb, protein
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna

# Align
align dna, protein

# Hide all initially
hide everything

# Show protein cartoon
show cartoon, protein
color gray70, protein

# Show DNA as cartoon
show cartoon, dna
color marine, dna
set cartoon_ring_mode, 3, dna
set cartoon_ladder_mode, 1, dna

# Highlight mutation site
select mut_site, protein and resi 175
show spheres, mut_site
color cyan, mut_site
set sphere_scale, 1.8, mut_site

# Show DNA binding interface as sticks
select interface, protein and resi 175+241+248+273+280
show sticks, interface
color orange, interface
set stick_radius, 0.35, interface

# Show hydrogen bonds to DNA
distance hbonds, interface, dna, 3.5
hide labels, hbonds
color lime, hbonds
set dash_width, 4, hbonds
set dash_gap, 0.1, hbonds

# Show binding pocket surface
select pocket, protein and (interface around 10)
show surface, pocket
set transparency, 0.6, pocket
color palegreen, pocket

# Highlight DNA bases near binding site
select dna_contact, dna and (interface around 5)
show sticks, dna_contact
color yellow, dna_contact and elem C
set stick_radius, 0.25, dna_contact

# Label
set label_size, 60
set label_color, blue
set label_bg_color, white
set label_bg_transparency, 0.3
pseudoatom label_obj, pos=[0,0,0]
label label_obj, 'WILD-TYPE\nNormal DNA Binding'

# ZOOM IN on binding interface
center interface
orient interface
zoom interface, 4

# Rotate to show binding clearly
turn y, 20
turn x, 10

# Render high-quality image
viewport 1920, 1080
ray 1920, 1080
png /Users/ishaangubbala/Documents/p53/reports/figures/dna_binding_closeups/R175H_wildtype_binding.png, dpi=300

print('Image saved!')
