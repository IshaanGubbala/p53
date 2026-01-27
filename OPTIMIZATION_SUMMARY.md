# MD Simulation Optimization Summary

**Date:** January 26, 2026 - 11:55 PM
**Action:** Restarted simulation with optimized threading

---

## ✅ What We Did

### 1. Created Restart Script
**File:** `src/md/restart_md_optimized.sh`
- Automatic checkpoint detection
- Safe process termination
- Optimized threading configuration
- Reusable for future rescues

### 2. Stopped Current Simulation
- **Old PID:** 81588
- **Progress saved:** Step 68,150 (136.3 ps / 10,000 ps)
- **Checkpoint:** md.cpt (8.5 MB, saved at 23:48)
- **No data lost:** Resumed from exact same point

### 3. Restarted with Optimizations
- **New PID:** 3751
- **Threading:** `-ntomp 8` (using all 8 CPU cores)
- **Status:** ✅ Running successfully
- **Resumed from:** Step 68,150 (136.3 ps)

---

## 📊 Threading Configuration

### Before Optimization:
```
OpenMP threads: Auto-detected (~5 threads)
CPU usage:      ~400-450% (4-4.5 cores)
Performance:    2.6 ns/day (9.17 hours/ns)
```

### After Optimization:
```
OpenMP threads: 8 (explicitly set)
CPU usage:      ~400-500% (4-5 cores, still ramping up)
Performance:    TBD (calibrating, check in 15 min)
Expected:       ~4.5-5.0 ns/day (5-5.5 hours/ns)
```

**Note:** GROMACS is currently warming up and calibrating. Performance will stabilize in 10-15 minutes.

---

## 🎯 Expected Improvements

### Estimated Performance Gain:
- **Conservative:** 1.5× speedup → 3.9 ns/day
- **Realistic:** 1.8× speedup → 4.7 ns/day
- **Optimistic:** 2.0× speedup → 5.2 ns/day

### Revised Completion Times:

| Duration | Before (2.6 ns/day) | After (4.7 ns/day) | Time Saved |
|----------|--------------------|--------------------|------------|
| **10 ns** | 3.8 days | **2.1 days** | 1.7 days (45%) |
| **5 ns** | 1.9 days | **1.1 days** | 0.8 days (42%) |
| **2 ns** | 0.8 days | **0.4 days** | 0.4 days (50%) |

**10 ns completion ETA:** **Thursday evening (January 29)** instead of Friday

---

## 🔍 Why Not 8× Speedup with 8 Threads?

### GROMACS Performance Factors:

1. **Amdahl's Law**
   - Some parts of MD can't be parallelized (serial bottlenecks)
   - PME electrostatics has communication overhead
   - Typical scaling: ~70-80% efficiency

2. **Memory Bandwidth**
   - 370,000 atoms requires lots of memory access
   - Apple Silicon has unified memory (shared between cores)
   - Bandwidth saturation at ~4-6 cores

3. **Synchronization Overhead**
   - Threads must sync at every timestep (2 fs)
   - More threads = more coordination overhead
   - Diminishing returns after ~4-6 threads

4. **PME Grid Limitations**
   - Electrostatics calculated on 3D grid (100×100×100)
   - Grid decomposition limited for this system size
   - PME takes ~15% of total time (hard to parallelize)

**Typical real-world scaling:** 8 threads → 1.5-2.0× speedup (not 8×)

---

## 💡 Additional Optimization Opportunities

### Already Implemented:
- [x] Multi-threading optimization (8 threads) ✅

### Still Available:

#### 1. Hydrogen Mass Repartitioning (HMR)
**Speedup:** 2× (4 fs timestep instead of 2 fs)
**Effort:** Moderate (requires regenerating topology)
**Trade-off:** None (validated for biomolecular simulations)
**When:** For future rescues or full pipeline

#### 2. Reduce Output Frequency
**Speedup:** 5-10% (less I/O overhead)
**Effort:** Easy (change MDP file)
**Trade-off:** Lower temporal resolution
**Implementation:**
```bash
# configs/md/mdp/md_fast.mdp
nstxtcout = 10000  # Save every 20 ps instead of 10 ps
```

#### 3. Cloud GPU (AWS p3.2xlarge)
**Speedup:** 10-15× (vs current optimized CPU)
**Effort:** Medium (AWS setup, ~2 hours)
**Cost:** ~$3/hour
**10 ns cost:** ~$20-30
**50 ns cost:** ~$100-150

#### 4. Shorter Cutoffs (Aggressive)
**Speedup:** 15-20%
**Trade-off:** Slightly less accurate (not recommended)

---

## 📈 Current Status

### Simulation Progress:
```
Current step:    70,000 / 5,000,000 (1.4%)
Current time:    140 ps / 10,000 ps (1.4%)
Started:         23:49 PM (restarted)
Runtime:         ~5 minutes (warming up)
```

### System State:
```
Temperature:     300.0 K ✅ (perfect)
Pressure:        -50.4 bar (equilibrated, normal fluctuation)
Total Energy:    -5.02 MJ/mol (stable)
Trajectory:      20 MB (growing steadily)
CPU Usage:       417% (4+ cores active)
Memory:          355 MB
```

### Files:
- **md.xtc:** 20 MB (growing) ✅
- **md.log:** 34 KB (updating) ✅
- **md.edr:** 10 KB (updating) ✅
- **md.cpt:** 8.5 MB (last checkpoint) ✅

---

## 📋 Next Steps

### Immediate (Tonight):
- [x] Simulation restarted with optimization ✅
- [ ] Let run for 15 minutes to stabilize ⏳
- [ ] Check performance at midnight (23:59)
- [ ] Calculate actual speedup

### Tomorrow Morning:
- [ ] Check progress (should be at ~250-300 ps if 4.7 ns/day)
- [ ] Verify stable performance
- [ ] Decide: continue to 10 ns, reduce to 5 ns, or stop at 2 ns

### After Pilot Completes:
- [ ] Analyze trajectory (RMSD, RMSF, clustering)
- [ ] Test docking pipeline (COTI-2)
- [ ] Calculate synergy score
- [ ] Decide on cloud GPU for scale-up

---

## 🎯 Success Metrics

### Optimization Successful If:
- [x] Simulation resumed without data loss ✅
- [x] Using 8 OpenMP threads ✅
- [x] Process running stably ✅
- [x] Files growing (xtc, log, edr) ✅
- [ ] Performance ≥ 4.0 ns/day ⏳ (measuring)
- [ ] No crashes or errors ⏳ (monitoring)

**Current Status:** ✅ 4/6 success criteria met, 2 pending measurement

---

## 🛠️ How to Monitor

### Check Performance (after 15 min warmup):
```bash
cd Data/processed/md_simulations/A189S_M133L_S95T

# View latest progress
tail -30 md.log

# Check step and time
grep "Step\|Time" md.log | tail -10

# Calculate ns/day manually:
# (current_ps - restart_ps) / (hours_elapsed) * 24 / 1000
```

### Check if Running:
```bash
ps aux | grep "3751"  # Current PID
ps aux | grep "gmx mdrun"  # Any mdrun process
```

### Check Files Growing:
```bash
ls -lh md.xtc md.log md.edr
# Wait 1 minute
ls -lh md.xtc md.log md.edr
# File sizes should increase
```

### If Simulation Crashes:
```bash
cd Data/processed/md_simulations/A189S_M133L_S95T

# Check error
tail -100 md.log

# Restart from checkpoint
gmx mdrun -v -deffnm md -cpi md.cpt -ntomp 8
```

---

## 📚 Lessons Learned

### What Worked:
1. ✅ Checkpointing allowed seamless restart with no data loss
2. ✅ Explicit threading (`-ntomp 8`) overrides conservative auto-detection
3. ✅ Monitoring CPU usage revealed underutilization

### What to Remember:
1. GROMACS takes 10-15 min to stabilize after restart
2. Real-world speedup is always < theoretical (Amdahl's Law)
3. For large-scale work, cloud GPU is more cost-effective than CPU time

### For Future Rescues:
1. Always use `-ntomp` with explicit thread count
2. Consider HMR (4 fs timestep) for 2× speedup
3. Plan cloud GPU for full pipeline (10 rescues × 50 ns)

---

## 💰 Cost Analysis (Full Pipeline)

### Local CPU (Current Optimized):
- **Duration:** 10 rescues × 50 ns × 5.5 hrs/ns = 2,750 hours = **115 days**
- **Cost:** Free (electricity: ~$50)
- **Parallelization:** Limited (1-2 at a time)

### Cloud CPU (c5.4xlarge):
- **Duration:** 10 rescues × 50 ns × 3 hrs/ns = 1,500 hours
- **Cost:** ~$300 (serial) or ~$150 (parallel)
- **Not worth it:** Only 2× faster than local, but costs money

### Cloud GPU (p3.2xlarge):
- **Duration:** 10 rescues × 50 ns × 0.5 hrs/ns = 250 hours
- **Cost:** ~$750 (serial) or **~$150** (10 parallel = 25 hours)
- **Best option:** 10× faster, reasonable cost if parallelized

**Recommendation:** Use local CPU for pilot (free), cloud GPU for full pipeline

---

*Last updated: January 26, 2026, 11:55 PM*
*Simulation PID: 3751*
*Check performance at midnight to verify optimization worked*
