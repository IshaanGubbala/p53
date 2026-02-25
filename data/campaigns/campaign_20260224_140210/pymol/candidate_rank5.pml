# PyMOL script for rescue candidate: Rank5
# Rescue mutations : E2M, L26S, P27S, A76T, A86M, S94F, V97L, T102N, L114K, A119W, V122A, C124W, Y126M, A129S, A138N, V143L, L145F, S149R, E171L, G187R, A189P, E204Q, F212G, Y220P, F270E, R273C, A276K, P278R, R282W, E294C, G334V, R337L, F338Y, F341L, A353I, E358V, A364Q, L369Q, Q375P, L383V
# Cancer mutations  : R273C, R282W
# Delivery method   : protein_therapy
# Evidence score    : 0.000

load data/raw/p53_wt.pdb, wt_p53
hide everything
show cartoon, wt_p53
color gray80, wt_p53
color lightblue, wt_p53 and resi 94-292
select cancer_muts, resi 273+282
show sticks, cancer_muts
color red, cancer_muts
set stick_radius, 0.25, cancer_muts
label resi 273 and name CA, "R273C"
label resi 282 and name CA, "R282W"
select rescue_muts, resi 2+26+27+76+86+94+97+102+114+119+122+124+126+129+138+143+145+149+171+187+189+204+212+220+270+273+276+278+282+294+334+337+338+341+353+358+364+369+375+383
show sticks, rescue_muts
color green, rescue_muts
set stick_radius, 0.25, rescue_muts
label resi 2 and name CA, "E2M"
label resi 26 and name CA, "L26S"
label resi 27 and name CA, "P27S"
label resi 76 and name CA, "A76T"
label resi 86 and name CA, "A86M"
label resi 94 and name CA, "S94F"
label resi 97 and name CA, "V97L"
label resi 102 and name CA, "T102N"
label resi 114 and name CA, "L114K"
label resi 119 and name CA, "A119W"
label resi 122 and name CA, "V122A"
label resi 124 and name CA, "C124W"
label resi 126 and name CA, "Y126M"
label resi 129 and name CA, "A129S"
label resi 138 and name CA, "A138N"
label resi 143 and name CA, "V143L"
label resi 145 and name CA, "L145F"
label resi 149 and name CA, "S149R"
label resi 171 and name CA, "E171L"
label resi 187 and name CA, "G187R"
label resi 189 and name CA, "A189P"
label resi 204 and name CA, "E204Q"
label resi 212 and name CA, "F212G"
label resi 220 and name CA, "Y220P"
label resi 270 and name CA, "F270E"
label resi 273 and name CA, "R273C"
label resi 276 and name CA, "A276K"
label resi 278 and name CA, "P278R"
label resi 282 and name CA, "R282W"
label resi 294 and name CA, "E294C"
label resi 334 and name CA, "G334V"
label resi 337 and name CA, "R337L"
label resi 338 and name CA, "F338Y"
label resi 341 and name CA, "F341L"
label resi 353 and name CA, "A353I"
label resi 358 and name CA, "E358V"
label resi 364 and name CA, "A364Q"
label resi 369 and name CA, "L369Q"
label resi 375 and name CA, "Q375P"
label resi 383 and name CA, "L383V"
set title, "Rank5 | rescue: E2M, L26S, P27S, A76T, A86M, S94F, V97L, T102N, L114K, A119W, V122A, C124W, Y126M, A129S, A138N, V143L, L145F, S149R, E171L, G187R, A189P, E204Q, F212G, Y220P, F270E, R273C, A276K, P278R, R282W, E294C, G334V, R337L, F338Y, F341L, A353I, E358V, A364Q, L369Q, Q375P, L383V | score: 0.000"
bg_color white
set ray_shadows, on
set ray_opaque_background, on
set label_size, 14
set label_color, black
set label_font_id, 7
set antialias, 2
set ray_trace_mode, 1
zoom resi 273+282+2+26+27+76+86+94+97+102+114+119+122+124+126+129+138+143+145+149+171+187+189+204+212+220+270+273+276+278+282+294+334+337+338+341+353+358+364+369+375+383, 15
deselect
