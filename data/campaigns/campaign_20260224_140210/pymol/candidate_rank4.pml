# PyMOL script for rescue candidate: Rank4
# Rescue mutations : L111V, C124W, A138N, C141L, L145F, V157I, V203A, N235K, Y236F, N239Q, C242S, R248W, R249S, T253V, I255F, G262E, N268R, V272A, V274I, C275S, A276S
# Cancer mutations  : R248W, R249S
# Delivery method   : gene_therapy
# Evidence score    : 0.000

load data/raw/p53_wt.pdb, wt_p53
hide everything
show cartoon, wt_p53
color gray80, wt_p53
color lightblue, wt_p53 and resi 94-292
select cancer_muts, resi 248+249
show sticks, cancer_muts
color red, cancer_muts
set stick_radius, 0.25, cancer_muts
label resi 248 and name CA, "R248W"
label resi 249 and name CA, "R249S"
select rescue_muts, resi 111+124+138+141+145+157+203+235+236+239+242+248+249+253+255+262+268+272+274+275+276
show sticks, rescue_muts
color green, rescue_muts
set stick_radius, 0.25, rescue_muts
label resi 111 and name CA, "L111V"
label resi 124 and name CA, "C124W"
label resi 138 and name CA, "A138N"
label resi 141 and name CA, "C141L"
label resi 145 and name CA, "L145F"
label resi 157 and name CA, "V157I"
label resi 203 and name CA, "V203A"
label resi 235 and name CA, "N235K"
label resi 236 and name CA, "Y236F"
label resi 239 and name CA, "N239Q"
label resi 242 and name CA, "C242S"
label resi 248 and name CA, "R248W"
label resi 249 and name CA, "R249S"
label resi 253 and name CA, "T253V"
label resi 255 and name CA, "I255F"
label resi 262 and name CA, "G262E"
label resi 268 and name CA, "N268R"
label resi 272 and name CA, "V272A"
label resi 274 and name CA, "V274I"
label resi 275 and name CA, "C275S"
label resi 276 and name CA, "A276S"
set title, "Rank4 | rescue: L111V, C124W, A138N, C141L, L145F, V157I, V203A, N235K, Y236F, N239Q, C242S, R248W, R249S, T253V, I255F, G262E, N268R, V272A, V274I, C275S, A276S | score: 0.000"
bg_color white
set ray_shadows, on
set ray_opaque_background, on
set label_size, 14
set label_color, black
set label_font_id, 7
set antialias, 2
set ray_trace_mode, 1
zoom resi 248+249+111+124+138+141+145+157+203+235+236+239+242+248+249+253+255+262+268+272+274+275+276, 15
deselect
