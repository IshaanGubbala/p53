# Initiative 2: Setup Guide

**Purpose:** Install required tools and prepare infrastructure for MD simulations + drug docking
**Target:** macOS (Darwin 25.3.0)
**Date:** January 26, 2026

---

## Required Software Installation

### 1. GROMACS (Molecular Dynamics Engine)

**Option A: Conda Installation (Recommended - Easy)**
```bash
# Create dedicated conda environment
conda create -n md_env python=3.10
conda activate md_env

# Install GROMACS
conda install -c conda-forge gromacs

# Verify installation
gmx --version
# Expected: GROMACS version 2023.x or later
```

**Option B: Homebrew (macOS)**
```bash
brew install gromacs

# Verify
gmx --version
```

**Option C: From Source (For GPU Support)**
```bash
# Prerequisites
brew install cmake fftw

# Download GROMACS
wget ftp://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz
tar xfz gromacs-2023.3.tar.gz
cd gromacs-2023.3

# Configure (CPU-only)
mkdir build && cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs

# Build and install
make -j 4
sudo make install

# Add to PATH
echo 'export PATH=/usr/local/gromacs/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

**GPU Acceleration (Optional but Recommended):**
- Requires NVIDIA GPU + CUDA Toolkit
- On macOS (M1/M2): Use conda version (Metal support via OpenCL)
- Speedup: ~10-20× faster than CPU

---

### 2. AutoDock Vina (Molecular Docking)

**Installation:**
```bash
# Activate MD environment
conda activate md_env

# Install AutoDock Vina
conda install -c conda-forge autodock-vina

# Verify
vina --version
# Expected: AutoDock Vina 1.2.x
```

---

### 3. Open Babel (Chemical Structure Manipulation)

**Installation:**
```bash
conda activate md_env
conda install -c conda-forge openbabel

# Verify
obabel --version
# Expected: Open Babel 3.x
```

---

### 4. Python Packages (Analysis & Automation)

**Installation:**
```bash
conda activate md_env

# Core packages (likely already installed)
pip install numpy pandas matplotlib seaborn

# MD analysis
pip install MDAnalysis nglview

# Chemical informatics
pip install rdkit biopython pubchempy

# Visualization
pip install pymol-open-source

# Verify
python -c "import MDAnalysis; print(MDAnalysis.__version__)"
```

---

## Cloud Computing Alternative (No Local Installation)

If local installation is challenging or GPU is unavailable, use cloud computing:

### Option 1: Google Colab (Free GPU)
```python
# In Colab notebook:
!apt-get install -y gromacs
!conda install -c conda-forge autodock-vina openbabel

# Upload input files, run simulations, download results
```

**Pros:** Free GPU (Tesla T4), no setup
**Cons:** 12-hour session limit, need to manage interruptions

### Option 2: AWS EC2 (Pay-as-you-go)
```bash
# Launch p3.2xlarge instance (Tesla V100 GPU)
# ~$3/hour, 100 GPU hours = $300 total

# SSH into instance
ssh -i key.pem ubuntu@ec2-instance

# Install GROMACS with GPU support (already compiled)
# Run all simulations in parallel
```

**Pros:** Reliable, scalable, professional
**Cons:** Costs money (~$300 for full Initiative 2)

### Option 3: University HPC Cluster (Free if Available)
- Contact university IT department
- Most universities have HPC clusters with GROMACS pre-installed
- Submit batch jobs via SLURM/PBS

---

## Directory Structure Setup

Run this script to create the required directory structure:

```bash
# Create directories
mkdir -p Data/raw/drugs/{structures,prepared}
mkdir -p Data/processed/md_simulations/{structures,topologies,trajectories,analysis,ensembles}
mkdir -p Data/processed/docking/{results,consensus}
mkdir -p Data/processed/synergy/{scores}
mkdir -p reports/initiative2/{figures,supplementary}
mkdir -p src/md
mkdir -p src/docking
mkdir -p configs/md

echo "✅ Directory structure created"
```

---

## Validation Tests

After installation, run these tests to verify everything works:

### Test 1: GROMACS Installation
```bash
# Generate a simple water box
gmx solvate -cs spc216.gro -o water.gro -box 3 3 3

# Expected: Creates water.gro with ~700 water molecules
# If successful, GROMACS is working correctly
```

### Test 2: AutoDock Vina Installation
```bash
# Test with example files
vina --help

# Should print usage information
```

### Test 3: Open Babel Installation
```bash
# Convert a sample molecule
echo "C1=CC=CC=C1" > benzene.smi
obabel benzene.smi -O benzene.pdb --gen3d

# Should create benzene.pdb with 3D coordinates
```

### Test 4: Python Packages
```bash
python << EOF
import MDAnalysis
import rdkit
import pubchempy
print("✅ All Python packages installed correctly")
EOF
```

---

## Quick Start: Pilot Run

Once tools are installed, test the full pipeline on one rescue:

```bash
# 1. Prepare structure (manual - use PyMOL or existing mutant from EvoEF2)
cp Data/processed/rescues/R175H/mutant_models/A189S_M133L_S95T.pdb \
   Data/processed/md_simulations/structures/

# 2. Run short MD (10 ns test)
cd Data/processed/md_simulations
gmx pdb2gmx -f structures/A189S_M133L_S95T.pdb -o processed.gro -water tip3p
gmx editconf -f processed.gro -o box.gro -c -d 1.0 -bt cubic
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top
gmx grompp -f ../../../configs/md/mdp/em.mdp -c solvated.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# Continue with equilibration and production...
# (Full scripts provided in src/md/)

# 3. Test docking
vina --receptor em.pdbqt --ligand ../../raw/drugs/prepared/COTI-2.pdbqt \
     --center_x 30 --center_y 25 --center_z 40 \
     --size_x 20 --size_y 20 --size_z 20 \
     --out docked.pdbqt

# 4. Analyze results
python src/md/analyze_trajectory.py
python src/docking/calculate_synergy.py
```

---

## Troubleshooting

### Issue: GROMACS "Cannot find force field"
**Solution:**
```bash
# Download AMBER force field
gmx pdb2gmx -h  # Check available force fields
# If AMBER99SB-ILDN not listed, install via conda or download manually
```

### Issue: Vina "Cannot find config file"
**Solution:**
```bash
# Create config.txt with all parameters
# Don't rely on command-line arguments alone
```

### Issue: Out of Memory during MD
**Solution:**
```bash
# Reduce system size: smaller box, fewer ions
gmx editconf -f box.gro -o smaller_box.gro -c -d 0.8  # Smaller padding
```

### Issue: Slow MD on CPU
**Solution:**
- Use cloud GPU (AWS p3.2xlarge)
- Or reduce simulation time: 25 ns instead of 50 ns
- Or reduce number of rescues: focus on top 5 instead of 10

---

## Status Checklist

Before starting Initiative 2 production runs, verify:

- [ ] GROMACS installed and tested (`gmx --version`)
- [ ] AutoDock Vina installed (`vina --version`)
- [ ] Open Babel installed (`obabel --version`)
- [ ] Python packages installed (MDAnalysis, RDKit, PubChemPy)
- [ ] Directory structure created
- [ ] GPU available (check with `nvidia-smi` or use cloud)
- [ ] Pilot run successful (1 rescue, 1 drug, short MD)

**Once all checkboxes are ✅, proceed to Phase 1: MD Simulations**

---

## Estimated Time & Cost

**Installation:** 1-2 hours (one-time setup)
**Pilot run:** 2-3 hours (validation)
**Full Initiative 2:** 2-3 weeks (as per roadmap)

**Cost (if using cloud GPU):**
- AWS p3.2xlarge: ~$300 for 100 GPU hours
- Google Colab Pro: $10/month (can complete Initiative 2 in 1 month)
- University HPC: Free (if available)

**Recommendation:** Start with local CPU for pilot, move to cloud GPU for production

---

*Setup guide prepared: January 26, 2026*
*Next: Install tools and run validation tests*
