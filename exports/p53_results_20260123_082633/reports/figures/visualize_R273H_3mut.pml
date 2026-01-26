# PyMOL script to visualize rescue mutations
# Target: R273H, Tier: 3mut
# Rescue mutations: R196Q, S215A, Y163F

# Load structure
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb
hide everything

# Show protein as cartoon
show cartoon, all
color gray80, all

# Highlight target mutation in red
select target, resi 273
show sticks, target
color red, target
label target and name CA, 'Target: R273H'

# Highlight rescue mutations in green
select rescue, resi 196+215+163
show sticks, rescue
color green, rescue
label resi 196 and name CA, 'Rescue 1: R196Q'
label resi 215 and name CA, 'Rescue 2: S215A'
label resi 163 and name CA, 'Rescue 3: Y163F'

# Show nearby residues for context
select nearby, (resi 273 around 8) or (resi 196+215+163 around 8)
show lines, nearby and not (target or rescue)

# Center view on target and rescues
zoom resi 273+196+215+163

# Set ray tracing for high-quality image
set ray_trace_mode, 1
set ray_shadows, 0

# Save image
png /Users/ishaangubbala/Documents/p53/reports/figures/structure_R273H_3mut.png, width=1200, height=1200, dpi=300, ray=1

print 'Visualization complete!'