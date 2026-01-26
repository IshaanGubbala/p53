# PyMOL Side-by-Side Comparison
# Target: R273H, Rescue: A189S

# Load structures
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb, wt
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R273H/cancer_build/R273H_cancer.pdb, mutant
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R273H/rescue_build/R273H_rescued.pdb, rescued
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna_wt
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna_mut
load /Users/ishaangubbala/Documents/p53/Data/raw/experimental_pdbs/2OCJ_full.pdb, dna_res

# Align all
align mutant, wt
align rescued, wt
align dna_wt, wt
align dna_mut, mutant
align dna_res, rescued

# Position structures side-by-side
translate [-40, 0, 0], wt
translate [-40, 0, 0], dna_wt
translate [40, 0, 0], rescued
translate [40, 0, 0], dna_res
# mutant stays at origin

# Style all proteins
hide everything
show cartoon, wt or mutant or rescued
show cartoon, dna_*
color gray70, wt or mutant or rescued
color blue, dna_*
set cartoon_fancy_helices, 1

# Wild-type: Cyan highlights
select wt_interface, wt and resi 273+241+248+273+280
show sticks, wt_interface
color cyan, wt_interface
distance hb_wt, wt_interface, dna_wt, 3.5
hide labels, hb_wt
color lime, hb_wt
set dash_width, 3, hb_wt

# Cancer mutant: Red highlights
select mut_site, mutant and resi 273
show spheres, mut_site
color red, mut_site
set sphere_scale, 1.8, mut_site
distance hb_mut, mutant, dna_mut, 3.5
hide labels, hb_mut
color red, hb_mut
set dash_width, 2, hb_mut
set dash_gap, 0.5, hb_mut

# Rescued: Green highlights
select res_site, rescued and resi 273
show spheres, res_site
color salmon, res_site
set sphere_scale, 1.6, res_site
select rescue_1, rescued and resi 189
show spheres, rescue_1
color green, rescue_1
set sphere_scale, 1.4, rescue_1
distance hb_res, rescued, dna_res, 3.5
hide labels, hb_res
color lime, hb_res
set dash_width, 3, hb_res

# Labels
pseudoatom label1, pos=[-40, 25, 0]
label label1, 'WILD-TYPE\nStable DNA Binding'
pseudoatom label2, pos=[0, 25, 0]
label label2, 'CANCER (R273H)\nDestabilized'
pseudoatom label3, pos=[40, 25, 0]
label label3, 'RESCUED\nA189S'

set label_size, 30
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.3

# Orient and render
orient all
zoom all, 5
viewport 3600, 1200
bg_color white
set ray_shadows, 1
set antialias, 2
ray 3600, 1200
png /Users/ishaangubbala/Documents/p53/reports/figures/animations/rescue_comparison_R273H.png, dpi=300

print('Side-by-side comparison complete!')