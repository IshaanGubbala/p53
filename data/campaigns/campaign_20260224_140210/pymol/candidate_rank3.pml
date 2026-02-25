# PyMOL script for rescue candidate: Rank3
# Rescue mutations : L111V, C124W, C135I, A138N, C141L, L145F, A159S, A161V, V203A, C229Y, H233L, N235K, Y236F, R248W, T253V, I255L, D259S, N268Q, V274I, C275S
# Cancer mutations  : R248W
# Delivery method   : protein_therapy
# Evidence score    : 0.000

load data/raw/p53_wt.pdb, wt_p53
hide everything
show cartoon, wt_p53
color gray80, wt_p53
color lightblue, wt_p53 and resi 94-292
select cancer_muts, resi 248
show sticks, cancer_muts
color red, cancer_muts
set stick_radius, 0.25, cancer_muts
label resi 248 and name CA, "R248W"
select rescue_muts, resi 111+124+135+138+141+145+159+161+203+229+233+235+236+248+253+255+259+268+274+275
show sticks, rescue_muts
color green, rescue_muts
set stick_radius, 0.25, rescue_muts
label resi 111 and name CA, "L111V"
label resi 124 and name CA, "C124W"
label resi 135 and name CA, "C135I"
label resi 138 and name CA, "A138N"
label resi 141 and name CA, "C141L"
label resi 145 and name CA, "L145F"
label resi 159 and name CA, "A159S"
label resi 161 and name CA, "A161V"
label resi 203 and name CA, "V203A"
label resi 229 and name CA, "C229Y"
label resi 233 and name CA, "H233L"
label resi 235 and name CA, "N235K"
label resi 236 and name CA, "Y236F"
label resi 248 and name CA, "R248W"
label resi 253 and name CA, "T253V"
label resi 255 and name CA, "I255L"
label resi 259 and name CA, "D259S"
label resi 268 and name CA, "N268Q"
label resi 274 and name CA, "V274I"
label resi 275 and name CA, "C275S"
set title, "Rank3 | rescue: L111V, C124W, C135I, A138N, C141L, L145F, A159S, A161V, V203A, C229Y, H233L, N235K, Y236F, R248W, T253V, I255L, D259S, N268Q, V274I, C275S | score: 0.000"
bg_color white
set ray_shadows, on
set ray_opaque_background, on
set label_size, 14
set label_color, black
set label_font_id, 7
set antialias, 2
set ray_trace_mode, 1
zoom resi 248+111+124+135+138+141+145+159+161+203+229+233+235+236+248+253+255+259+268+274+275, 15
deselect
