# PyMOL script to visualize rescue mutations
# Target: R248Q, Tier: 3mut
# Rescue mutations: M133L, R196Q, R213Q

# Load structure
load /Users/ishaangubbala/Documents/p53/Data/raw/alphafold/AF-P04637-F1-model_v6.pdb
hide everything

# Show protein as cartoon
show cartoon, all
color gray80, all

# Highlight target mutation in red
select target, resi 248
show sticks, target
color red, target
label target and name CA, 'Target: R248Q'

# Highlight rescue mutations in green
select rescue, resi 133+196+213
show sticks, rescue
color green, rescue
label resi 133 and name CA, 'Rescue 1: M133L'
label resi 196 and name CA, 'Rescue 2: R196Q'
label resi 213 and name CA, 'Rescue 3: R213Q'

# Show nearby residues for context
select nearby, (resi 248 around 8) or (resi 133+196+213 around 8)
show lines, nearby and not (target or rescue)

# Center view on target and rescues
zoom resi 248+133+196+213

# Set ray tracing for high-quality image
set ray_trace_mode, 1
set ray_shadows, 0

# Save image
png /Users/ishaangubbala/Documents/p53/reports/figures/structure_R248Q_3mut.png, width=1200, height=1200, dpi=300, ray=1

print 'Visualization complete!'