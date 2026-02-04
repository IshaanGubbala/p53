load data/raw/p53_wt.pdb, wt_p53
hide everything
show cartoon, wt_p53
color gray70, wt_p53
set transparency, 0.4
bg_color white
select muts, resi 1+2+66+81+91+393
show sticks, muts
color red, muts
set stick_radius, 0.3
zoom muts, 20
set ray_opaque_background, off
set ray_shadows, on
png data/processed/best_design.png, width=800, height=600, ray=1
quit