#!/bin/bash
# p53-proteoMgCAD MD Validation - GROMACS
# Generated: 2026-02-25T07:03:09.441820
# Target: R273C_R282W

NAME="R273C_R282W"
SEQUENCE="MMEPQSDPSVEPPLSQETFSDLWKLSSENNVLSPLPSQAMDDLMLSPDDI..."  # Truncated

echo "=== p53-proteoMgCAD MD Validation ==="
echo "Target: $NAME"

# 1. Generate topology
gmx pdb2gmx -f ${NAME}.pdb -o ${NAME}_processed.gro -water tip3p -ff amber14

# 2. Define box
gmx editconf -f ${NAME}_processed.gro -o ${NAME}_box.gro -c -d 1.0 -bt cubic

# 3. Add solvent
gmx solvate -cp ${NAME}_box.gro -cs spc216.gro -o ${NAME}_solv.gro -p topol.top

# 4. Add ions
gmx grompp -f ions.mdp -c ${NAME}_solv.gro -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o ${NAME}_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# 5. Energy minimization
gmx grompp -f em.mdp -c ${NAME}_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# 6. NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# 7. NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# 8. Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md

# 9. Analysis
echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o ${NAME}_noPBC.xtc -pbc mol -center
echo "Backbone Backbone" | gmx rms -s md.tpr -f ${NAME}_noPBC.xtc -o rmsd.xvg -tu ns
echo "Backbone" | gmx rmsf -s md.tpr -f ${NAME}_noPBC.xtc -o rmsf.xvg -res

echo "=== Simulation Complete ==="
echo "RMSD: rmsd.xvg"
echo "RMSF: rmsf.xvg"
