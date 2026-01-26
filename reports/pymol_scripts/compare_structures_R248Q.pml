# PyMOL Structural Comparison: Real Conformational Changes
# Target: R248Q, Rescue: M133L

# Settings for high-quality render
set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set ambient, 0.5
set specular, 0.6
set depth_cue, 1
set bg_rgb, [1, 1, 1]
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 0.8

# Load all three structures
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R248Q/wild_type.pdb, wt
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R248Q/cancer_build/R248Q_cancer.pdb, cancer
load /Users/ishaangubbala/Documents/p53/Data/processed/modeled_structures/R248Q/rescue_build/R248Q_rescued.pdb, rescued

# Align all to wild-type for comparison
align cancer, wt
align rescued, wt

# Position structures side-by-side
translate [-50, 0, 0], wt
# cancer stays at origin
translate [50, 0, 0], rescued

# Hide everything initially
hide everything

# === WILD-TYPE (Left) ===
show cartoon, wt
color gray70, wt

# Highlight the position that WILL be mutated (for comparison)
select wt_future_mut, wt and resi 248
show sticks, wt_future_mut
color cyan, wt_future_mut
set stick_radius, 0.25, wt_future_mut

# Show stable structure indicator (green surface)
select wt_region, wt and (resi 248 around 8)
show surface, wt_region
set transparency, 0.7, wt_region
color palegreen, wt_region

# === CANCER MUTANT (Middle) ===
show cartoon, cancer
color gray70, cancer

# Highlight cancer mutation site (RED)
select cancer_mut, cancer and resi 248
show spheres, cancer_mut
color red, cancer_mut
set sphere_scale, 1.8, cancer_mut
label cancer_mut and name CA, 'R248Q (R→Q)'

# Show destabilized region (red surface)
select cancer_destab, cancer and (resi 248 around 8)
show surface, cancer_destab
set transparency, 0.6, cancer_destab
color salmon, cancer_destab

# Show any structural distortions as sticks
select cancer_distorted, cancer and (resi 248+243+244+245+246+247+248+249+250+251+252+253)
show sticks, cancer_distorted
color orange, cancer_distorted
set stick_radius, 0.2, cancer_distorted
set stick_transparency, 0.5, cancer_distorted

# === RESCUED (Right) ===
show cartoon, rescued
color gray70, rescued

# Highlight rescue mutations (GREEN)
select rescue_mut_1, rescued and resi 133
show spheres, rescue_mut_1
color green, rescue_mut_1
set sphere_scale, 1.5, rescue_mut_1
label rescue_mut_1 and name CA, 'M133L'

# Show cancer mutation still present (but stabilized)
select rescued_cancer, rescued and resi 248
show spheres, rescued_cancer
color yellow, rescued_cancer
set sphere_scale, 1.4, rescued_cancer
label rescued_cancer and name CA, 'R248Q'

# Show re-stabilized region (green surface)
select rescued_stable, rescued and (resi 248 around 8)
show surface, rescued_stable
set transparency, 0.7, rescued_stable
color lightgreen, rescued_stable

# === CALCULATE RMSD (shows actual structural differences) ===
align cancer, wt, cycles=0
print('='*60)
print('STRUCTURAL DIFFERENCES (RMSD):')
print('='*60)

# Reset positions for RMSD calc
translate [50, 0, 0], wt
translate [0, 0, 0], cancer
translate [-50, 0, 0], rescued

# RMSD: wild-type vs cancer (shows destabilization)
align cancer, wt
print('Wild-type → Cancer mutant RMSD (destabilization):')

# RMSD: cancer vs rescued (shows rescue effect)
align rescued, cancer
print('Cancer mutant → Rescued RMSD (re-stabilization):')

# RMSD: wild-type vs rescued (how close to original)
align rescued, wt
print('Wild-type → Rescued RMSD (recovery):')
print('='*60)

# Reposition for visualization
translate [-50, 0, 0], wt
translate [0, 0, 0], cancer
translate [50, 0, 0], rescued

# === LABELS ===
set label_size, 35
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.3

# Title labels
pseudoatom label_wt, pos=[-50, 30, 0]
label label_wt, 'WILD-TYPE\nStable'

pseudoatom label_cancer, pos=[0, 30, 0]
label label_cancer, 'CANCER\nR248Q Destabilized'

pseudoatom label_rescued, pos=[50, 30, 0]
label label_rescued, 'RESCUED\nM133L'

# === ORIENT AND RENDER ===
orient all
zoom all, 6

# Set viewport for widescreen
viewport 3600, 1200

# Ray trace at high quality
print('Rendering high-quality image (this takes ~2-3 min)...')
ray 3600, 1200
png /Users/ishaangubbala/Documents/p53/reports/figures/structural_comparisons/R248Q_structural_comparison.png, dpi=300

print('='*60)
print('Image saved: /Users/ishaangubbala/Documents/p53/reports/figures/structural_comparisons/R248Q_structural_comparison.png')
print('='*60)