# p53-DNA Interaction Video
# Cancer mutant: R175H

# Settings
viewport 1280, 720
set bg_rgb, [1, 1, 1]
set antialias, 1
set cartoon_fancy_helices, 1

# Load structures
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R175H/cancer_build/R175H_cancer.pdb, protein
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna

# Align DNA to protein
align dna, protein

# Style protein
hide everything
show cartoon, protein
color gray70, protein

# Show DNA as sticks
show cartoon, dna
color blue, dna

# Highlight mutation site
select mut_site, protein and resi 175
show spheres, mut_site
color red, mut_site
set sphere_scale, 1.6, mut_site

# Highlight DNA binding interface
select interface, protein and resi 241+248+273+280
show sticks, interface
color cyan, interface
set stick_radius, 0.3, interface

# Show hydrogen bonds (key difference!)
distance hbonds, interface, dna, 3.5
hide labels, hbonds
color red, hbonds
set dash_width, 2, hbonds
set dash_gap, 0.3, hbonds

# Add surface to show binding pocket
select pocket, protein and (interface around 8)
show surface, pocket
set transparency, 0.75, pocket
color salmon, pocket

# Label
pseudoatom label, pos=[0,40,0]
label label, 'CANCER MUTANT (Disrupted Binding)'
set label_size, 50
set label_color, red

# Orient to show binding interface
orient interface
zoom interface, 6

# Create 60-frame movie (2 seconds @ 30fps)
mset 1 x60

# Smooth 360-degree rotation
mview store, 1
turn y, 360
mview store, 60
mview interpolate

# Render frames (no ray-tracing for speed)
set cache_frames, 0
mpng /Users/ishaangubbala/Documents/p53/reports/figures/animations/frames/R175H_cancer_binding_

print('='*60)
print('Frames rendered!')
print('Create video:')
print('cd /Users/ishaangubbala/Documents/p53/reports/figures/animations')
print('ffmpeg -framerate 30 -i frames/R175H_cancer_binding_%04d.png -c:v libx264 -pix_fmt yuv420p R175H_cancer_binding.mp4')
print('='*60)