# PyMOL script to visualize rescue mutations
# Target: R175H, Tier: 3mut
# Rescue mutations: A189S, M133L, Y163F

# Load structure
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb
hide everything

# Show protein as cartoon
show cartoon, all
color gray80, all

# Highlight target mutation in red
select target, resi 175
show sticks, target
color red, target
label target and name CA, 'Target: R175H'

# Highlight rescue mutations in green
select rescue, resi 189+133+163
show sticks, rescue
color green, rescue
label resi 189 and name CA, 'Rescue 1: A189S'
label resi 133 and name CA, 'Rescue 2: M133L'
label resi 163 and name CA, 'Rescue 3: Y163F'

# Show nearby residues for context
select nearby, (resi 175 around 8) or (resi 189+133+163 around 8)
show lines, nearby and not (target or rescue)

# Center view on target and rescues
zoom resi 175+189+133+163

# Set ray tracing for high-quality image
set ray_trace_mode, 1
set ray_shadows, 0

# Save image
png /Users/ishaangubbala/Documents/p53/reports/figures/structure_R175H_3mut.png, width=1200, height=1200, dpi=300, ray=1

print 'Visualization complete!'