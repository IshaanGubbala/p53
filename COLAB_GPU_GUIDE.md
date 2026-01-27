# Running GROMACS MD on Google Colab (FREE GPU)

**Goal:** Complete 10 ns MD simulation on free GPU in ~6-10 hours (vs 5 days on local CPU)

---

## ⚡ Quick Summary

**Cost:** FREE
**GPU:** NVIDIA T4 (~10-15× faster than your CPU)
**Session limit:** 12 hours
**Setup time:** 30-60 minutes
**10 ns runtime:** 6-10 hours
**Babysitting:** Need to restart every ~10-12 hours

---

## 📋 Step-by-Step Setup

### Step 1: Upload Files to Google Drive (5 min)

```bash
# On your local machine, create a zip of current checkpoint
cd /Users/ishaangubbala/Documents/p53/Data/processed/md_simulations/A189S_M133L_S95T

# Package checkpoint files
zip -r colab_restart.zip md.tpr md.cpt system.gro topol.top *.mdp

# Size should be ~50-100 MB
ls -lh colab_restart.zip
```

Then upload `colab_restart.zip` to your Google Drive.

---

### Step 2: Open Google Colab (2 min)

1. Go to: https://colab.research.google.com/
2. Click "New Notebook"
3. Change runtime to GPU:
   - Runtime → Change runtime type → GPU (T4) → Save

---

### Step 3: Install GROMACS with GPU Support (15-20 min)

Run this in a Colab cell:

```python
# Install GROMACS with CUDA support
!apt-get update -qq
!apt-get install -y gromacs

# Verify GPU is available
!nvidia-smi

# Check GROMACS version
!gmx --version | grep -E "VERSION|GPU"
```

**Note:** This installs a pre-compiled GROMACS with CUDA support from Ubuntu repos.

---

### Step 4: Mount Google Drive and Extract Files (5 min)

```python
from google.colab import drive
drive.mount('/content/drive')

# Extract checkpoint files
!cd /content && unzip /content/drive/MyDrive/colab_restart.zip
!ls -lh *.tpr *.cpt *.gro *.top
```

---

### Step 5: Run MD with GPU (6-10 hours for 10 ns)

```python
# Run GROMACS with GPU acceleration
!gmx mdrun -v -deffnm md -cpi md.cpt -ntomp 2 -nb gpu -pme gpu -update gpu

# The flags:
# -ntomp 2: Use 2 CPU threads (Colab gives 2 vCPUs)
# -nb gpu: Non-bonded on GPU
# -pme gpu: Electrostatics on GPU
# -update gpu: Coordinate updates on GPU
```

**Monitor progress:** This cell will run for hours and show live output. You can see:
```
Step 100000, time 200.0 ps
Performance: 28.5 ns/day
```

---

### Step 6: Download Results (10 min)

After completion:

```python
# Compress trajectory
!zip md_results.zip md.xtc md.log md.edr md.cpt

# Copy to Google Drive
!cp md_results.zip /content/drive/MyDrive/
```

Then download from Google Drive to your local machine.

---

## ⚠️ Handling 12-Hour Limit

If 10 ns takes >12 hours, Colab will disconnect. **Solution:**

### Before Disconnect (every 10-11 hours):

```python
# Manually stop simulation (Interrupt execution)
# Checkpoint is auto-saved every 15 min

# Upload latest checkpoint to Drive
!cp md.cpt /content/drive/MyDrive/md_restart.cpt
!cp md.log /content/drive/MyDrive/
!cp md.xtc /content/drive/MyDrive/

print("Checkpoint saved! Can safely disconnect.")
```

### After Reconnect:

```python
# Download checkpoint
!cp /content/drive/MyDrive/md_restart.cpt ./md.cpt

# Resume from checkpoint
!gmx mdrun -v -deffnm md -cpi md.cpt -ntomp 2 -nb gpu -pme gpu -update gpu
```

---

## 📊 Expected Performance

| System | Hardware | Performance | 10 ns Time |
|--------|----------|-------------|------------|
| **Current (local)** | Apple Silicon CPU | 2.0 ns/day | 5 days |
| **Colab Free (T4)** | Tesla T4 GPU | 20-30 ns/day | 8-12 hours |
| **Colab Pro (V100)** | Tesla V100 GPU | 40-60 ns/day | 4-6 hours |

---

## 💰 Cost Comparison

| Option | Cost | 10 ns Time | Full Pipeline (10×50ns) |
|--------|------|------------|-------------------------|
| **Local CPU** | Free | 5 days | 250 days 😱 |
| **Colab Free** | Free | 8-12 hours | ~5-7 days (if babysit) |
| **Colab Pro** | $10/month | 4-6 hours | ~2-3 days |
| **AWS GPU** | ~$27 | 9 hours | $135 (if parallel) |

---

## 🎯 My Recommendation for You

### For Pilot (10 ns):
**Use Colab Free** - Worth the setup time to save 5 days

### For Full Pipeline (10 rescues × 50 ns):
**Options ranked:**
1. **Colab Pro ($10)** - Best value if you can babysit
2. **Student credits (free)** - If you have .edu email
3. **AWS GPU ($135)** - If you want hands-off automation

---

## 🚀 Quick Start Commands

### On Your Mac (Prepare Files):
```bash
cd /Users/ishaangubbala/Documents/p53/Data/processed/md_simulations/A189S_M133L_S95T

# Stop current simulation
ps aux | grep gmx
kill [PID]  # Replace with actual PID

# Wait for checkpoint to save
sleep 30

# Package files
zip colab_restart.zip md.tpr md.cpt topol.top system.gro

# Upload to Google Drive manually or use:
# (requires rclone setup)
```

### On Colab:
```python
# Cell 1: Setup
!apt-get update -qq && apt-get install -y gromacs
!nvidia-smi

# Cell 2: Mount Drive
from google.colab import drive
drive.mount('/content/drive')

# Cell 3: Extract
!unzip /content/drive/MyDrive/colab_restart.zip

# Cell 4: Run
!gmx mdrun -v -deffnm md -cpi md.cpt -ntomp 2 -nb gpu -pme gpu -update gpu

# Cell 5: Save Results
!zip md_results.zip md.xtc md.log md.edr md.cpt
!cp md_results.zip /content/drive/MyDrive/
```

---

## 🐛 Troubleshooting

### "GPU not available"
- Check Runtime → Change runtime type → GPU → Save
- Restart runtime

### "GROMACS not found"
- Run: `!apt-get install -y gromacs`

### "Checkpoint file corrupt"
- Use `md_prev.cpt` instead of `md.cpt`

### "Session disconnected"
- Colab disconnects after 12 hours idle
- Save checkpoints to Drive every few hours
- Use Colab Pro for 24-hour sessions

---

## ⏱️ Time Estimate

**Total time to complete 10 ns on Colab:**
- Setup: 30 min
- Upload files: 5 min
- Install GROMACS: 15 min
- Run simulation: 8-12 hours
- Download results: 10 min
**Total: ~9-13 hours** (vs 5 days on local)

**Worth it?** Absolutely! ✅

---

## 📚 Alternative: Student Credits

If you have a **.edu email**, get **$300 free credits**:

1. Go to: https://cloud.google.com/edu/students
2. Sign up with student email
3. Get $300 credit (lasts 1 year)
4. Launch GPU instance (same as AWS p3.2xlarge)
5. Run everything with no session limits

**$300 = 100 hours of GPU = entire PhD project** 🎓

---

*Next steps: Let me know if you want to try Colab, and I can help with the setup!*
