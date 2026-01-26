# Close-up: R273H Cancer DNA Binding

set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set ambient, 0.5
set specular, 0.6
set depth_cue, 1
set bg_rgb, [1, 1, 1]
set cartoon_fancy_helices, 1

load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R273H/cancer_build/R273H_cancer.pdb, protein
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna

align dna, protein
hide everything

show cartoon, protein
color gray70, protein

show cartoon, dna
color marine, dna
set cartoon_ring_mode, 3, dna
set cartoon_ladder_mode, 1, dna

select mut_site, protein and resi 273
show spheres, mut_site
color red, mut_site
set sphere_scale, 1.8, mut_site

select interface, protein and resi 273+241+248+280
show sticks, interface
color orange, interface
set stick_radius, 0.35, interface

distance hbonds, interface, dna, 3.5
hide labels, hbonds
color red, hbonds
set dash_width, 2, hbonds
set dash_gap, 0.4, hbonds

select pocket, protein and (interface around 10)
show surface, pocket
set transparency, 0.6, pocket
color salmon, pocket

select dna_contact, dna and (interface around 5)
show sticks, dna_contact
color yellow, dna_contact and elem C
set stick_radius, 0.25, dna_contact

set label_size, 60
set label_color, red
set label_bg_color, white
set label_bg_transparency, 0.3
pseudoatom label_obj, pos=[0,0,0]
label label_obj, 'CANCER MUTANT\nR273H - Disrupted Binding'

center interface
orient interface
zoom interface, 4
turn y, 20
turn x, 10

viewport 1920, 1080
ray 1920, 1080
png /Users/ishaangubbala/Documents/p53/reports/figures/dna_binding_closeups/R273H_cancer_binding.png, dpi=300
print('Image saved!')
