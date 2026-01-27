# 🚀 Colab Setup - Quick Start Guide

**Goal:** Run your MD simulation on free GPU (8-12 hours vs 5 days on CPU)

---

## ✅ What's Ready

1. **✅ Simulation stopped** (saved at 170 ps)
2. **✅ Files packaged** (16 MB): `colab_restart.zip`
3. **✅ Notebook created**: `GROMACS_MD_Colab.ipynb`

**Location:** `/Users/ishaangubbala/Documents/p53/Data/processed/md_simulations/A189S_M133L_S95T/colab_restart.zip`

---

## 📋 5-Minute Setup Steps

### Step 1: Upload Files to Google Drive (2 min)

**Option A: Manual Upload (Easiest)**
1. Open Google Drive in browser: https://drive.google.com
2. Click "New" → "File upload"
3. Navigate to: `/Users/ishaangubbala/Documents/p53/Data/processed/md_simulations/A189S_M133L_S95T/`
4. Upload `colab_restart.zip` (16 MB)
5. Wait for upload to complete

**✅ Done when:** You see `colab_restart.zip` in your Google Drive root folder

---

### Step 2: Open Google Colab (1 min)

1. Go to: https://colab.research.google.com/
2. Click **"Upload"** tab
3. Click **"Choose File"**
4. Select: `/Users/ishaangubbala/Documents/p53/GROMACS_MD_Colab.ipynb`

**✅ Done when:** Notebook opens in Colab

---

### Step 3: Enable GPU (30 seconds)

1. Click **"Runtime"** (top menu)
2. Click **"Change runtime type"**
3. Hardware accelerator → Select **"GPU"** (T4)
4. Click **"Save"**

**✅ Done when:** You see "GPU" in the top-right corner

---

### Step 4: Run the Notebook (10-15 min setup, then 8-12 hours simulation)

**Just click through each cell from top to bottom:**

1. **Cell 1:** Verify GPU (should see NVIDIA GPU info)
2. **Cell 2:** Install GROMACS (~15-20 min) ☕
3. **Cell 3:** Mount Google Drive (authorize when prompted)
4. **Cell 4:** Extract files (should see your files listed)
5. **Cell 5:** Run simulation (~8-12 hours) ⏰
6. **Cell 6:** Check results
7. **Cell 7:** Save to Google Drive

**✅ Done when:** Cell 7 completes and you see "Results saved to Google Drive!"

---

### Step 5: Download Results (5 min)

1. Go to Google Drive: https://drive.google.com
2. Find `md_results.zip` (~500-1000 MB)
3. Right-click → Download
4. Save to your Mac

---

## ⏰ Timeline

| Step | Time | What You Do |
|------|------|-------------|
| 1. Upload to Drive | 2 min | Upload file |
| 2. Open Colab | 1 min | Open notebook |
| 3. Enable GPU | 30 sec | Change settings |
| 4. Setup (Cells 1-4) | 20 min | Run cells, authorize |
| 5. Simulation (Cell 5) | **8-12 hours** | **Wait (can close browser)** |
| 6. Save (Cell 7) | 5 min | Run cell |
| 7. Download | 5 min | Download from Drive |
| **Total** | **~9-13 hours** | **(vs 5 days on CPU!)** |

---

## 🎯 What to Expect During Simulation

**Cell 5 will show live output like this:**

```
Step 100000, time 200.0 ps
Performance: 28.5 ns/day, 0.842 hours/ns
Temperature 300.0 K, Pressure -12.3 bar

Step 150000, time 300.0 ps
Performance: 29.2 ns/day, 0.822 hours/ns
Temperature 299.8 K, Pressure 8.5 bar

... (continues for 8-12 hours)

Step 5000000, time 10000.0 ps
Finished mdrun on rank 0
```

**You can:**
- ✅ Close browser (simulation continues)
- ✅ Put computer to sleep
- ✅ Come back later and check progress
- ⚠️ Don't close the Colab tab completely (minimize is OK)

---

## 🔄 If Colab Disconnects (>12 hours)

Colab free tier has 12-hour limit. If simulation isn't done:

**Before timeout (around 11 hours):**
1. Stop Cell 5 (click stop button)
2. Run checkpoint save cell (in notebook)
3. Checkpoint saved to Drive ✅

**After reconnecting:**
1. Re-run Cells 1-4 (setup again)
2. Run resume cell (in notebook)
3. Simulation continues from checkpoint ✅

**Tip:** Colab Pro ($10/month) has 24-hour sessions = no interruption

---

## 🐛 Common Issues

### "No GPU available"
→ Check Runtime → Change runtime type → GPU is selected

### "File not found: colab_restart.zip"
→ Make sure file is in Google Drive root folder (not in subfolder)

### "gmx: command not found"
→ Re-run Cell 2 (install GROMACS)

### "Session disconnected"
→ Save checkpoint (see "If Colab Disconnects" above)

---

## 📊 What You'll Get

After completion, `md_results.zip` will contain:

```
md.xtc      ~500-1000 MB    Trajectory (10 ns of protein dynamics)
md.log      ~1-5 MB         Simulation log (energy, temperature, etc.)
md.edr      ~50-100 MB      Energy data
md.cpt      ~8 MB           Final checkpoint
```

**Next steps after download:**
1. Extract on your Mac
2. Run trajectory analysis: `python src/md/analyze_trajectory.py`
3. Continue with docking pipeline

---

## 💰 Cost Breakdown

| Option | Cost | Time | Effort |
|--------|------|------|--------|
| **Local CPU (current)** | Free | 5 days | Low |
| **Colab Free (this)** | **Free** | **12 hours** | **Medium** |
| **Colab Pro** | $10/month | 12 hours | Low |
| **AWS GPU** | ~$27 | 9 hours | High |

**Best choice:** Try Colab Free first! If it works, use it for full pipeline too.

---

## ✅ Checklist

Before you start:
- [ ] `colab_restart.zip` uploaded to Google Drive root
- [ ] Notebook opened in Colab
- [ ] GPU enabled in Runtime settings
- [ ] Google account authorized

During simulation:
- [ ] Cell 5 running (shows live progress)
- [ ] Browser tab open (or minimized)
- [ ] Internet connected (or will auto-resume)

After completion:
- [ ] Cell 7 executed (saved to Drive)
- [ ] `md_results.zip` downloaded
- [ ] Files extracted on local machine

---

## 🎉 You're Ready!

**Start here:** https://colab.research.google.com/

**Upload this file:** `/Users/ishaangubbala/Documents/p53/GROMACS_MD_Colab.ipynb`

**Upload this to Drive first:** `/Users/ishaangubbala/Documents/p53/Data/processed/md_simulations/A189S_M133L_S95T/colab_restart.zip`

---

**Questions? Issues?** Check the troubleshooting section in the notebook or ask!

**Good luck!** 🚀 You're about to save 4+ days of waiting time!
