load data/raw/p53_wt.pdb, wt_p53
hide everything
show cartoon, wt_p53
color gray70, wt_p53
set transparency, 0.4
bg_color white
select muts, resi 24+115+120+133+141+156+160+174+175+178+179+185+186+191+205+209+214+217+229+230+236+237+238+241+243+245+246+248+253+256+258+267+272+273+274+275+277+279+280+281+286+287+313+315
show sticks, muts
color red, muts
set stick_radius, 0.3
zoom muts, 20
set ray_opaque_background, off
set ray_shadows, on
png data/processed/cand_2.png, width=800, height=600, ray=1
quit