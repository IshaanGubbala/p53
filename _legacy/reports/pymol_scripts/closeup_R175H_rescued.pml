# Close-up: R175H Rescued DNA Binding

set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set ambient, 0.5
set specular, 0.6
set depth_cue, 1
set bg_rgb, [1, 1, 1]
set cartoon_fancy_helices, 1

load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R175H/rescue_build/R175H_rescued.pdb, protein
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna

align dna, protein
hide everything

show cartoon, protein
color gray70, protein

show cartoon, dna
color marine, dna
set cartoon_ring_mode, 3, dna
set cartoon_ladder_mode, 1, dna

select mut_site, protein and resi 175
show spheres, mut_site
color yellow, mut_site
set sphere_scale, 1.8, mut_site

select interface, protein and resi 175+241+248+273+280
show sticks, interface
color orange, interface
set stick_radius, 0.35, interface

distance hbonds, interface, dna, 3.5
hide labels, hbonds
color lime, hbonds
set dash_width, 4, hbonds
set dash_gap, 0.1, hbonds

select pocket, protein and (interface around 10)
show surface, pocket
set transparency, 0.6, pocket
color lightgreen, pocket

select dna_contact, dna and (interface around 5)
show sticks, dna_contact
color yellow, dna_contact and elem C
set stick_radius, 0.25, dna_contact

set label_size, 60
set label_color, green
set label_bg_color, white
set label_bg_transparency, 0.3
pseudoatom label_obj, pos=[0,0,0]
label label_obj, 'RESCUED\nBinding Restored'

center interface
orient interface
zoom interface, 4
turn y, 20
turn x, 10

viewport 1920, 1080
ray 1920, 1080
png /Users/ishaangubbala/Documents/p53/reports/figures/dna_binding_closeups/R175H_rescued_binding.png, dpi=300
print('Image saved!')
