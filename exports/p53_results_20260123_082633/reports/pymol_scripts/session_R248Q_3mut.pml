# PyMOL Interactive Session for Rescue Visualization
# Target: R248Q, Tier: 3mut
# Rescue mutations: M133L, R196Q, R213Q

# Load structure
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb

# Initial setup
bg_color white
hide everything
show cartoon, all
color gray80, all
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# Color by secondary structure
color salmon, ss h
color lightblue, ss s
color gray80, ss l+''

# Highlight target mutation in RED
select target_mut, resi 248
show sticks, target_mut
show spheres, target_mut
color red, target_mut
set sphere_scale, 0.3, target_mut

# Highlight rescue mutations in GREEN
select rescue_muts, resi 133+196+213
show sticks, rescue_muts
show spheres, rescue_muts
color green, rescue_muts
set sphere_scale, 0.3, rescue_muts

# Show nearby context (within 8Å)
select context, (resi 248 around 8) or (resi 133+196+213 around 8)
show lines, context and not (target_mut or rescue_muts)
color gray70, context and not (target_mut or rescue_muts)

# Add labels
label target_mut and name CA, '"Target: R248Q"'
label resi 133 and name CA, '"Rescue 1: M133L"'
label resi 196 and name CA, '"Rescue 2: R196Q"'
label resi 213 and name CA, '"Rescue 3: R213Q"'

# Set viewing angle
zoom resi 248+133+196+213
orient

# Enable better rendering
set ray_trace_mode, 1
set ray_shadows, 0
set antialias, 2
set cartoon_side_chain_helper, on

# Create selections for easy toggling
disable all
enable all

# Print instructions
print '======================================'
print 'PyMOL Interactive Session Loaded'
print '======================================'
print 'Target: R248Q'
print 'Rescue: M133L, R196Q, R213Q'
print ''
print 'Controls:'
print '  - Mouse drag: rotate'
print '  - Mouse scroll: zoom'
print '  - Right-click drag: translate'
print ''
print 'Useful commands:'
print '  hide everything; show cartoon'
print '  show surface, all'
print '  color_by_element'
print '  ray 1200,1200'
print '======================================'