#!/bin/bash
# Restart MD simulation with optimized settings
# Usage: ./restart_md_optimized.sh <rescue_name>

set -e

RESCUE=$1
MD_DIR="Data/processed/md_simulations/${RESCUE}"

if [ -z "$RESCUE" ]; then
    echo "❌ Usage: $0 <rescue_name>"
    echo "   Example: $0 A189S_M133L_S95T"
    exit 1
fi

if [ ! -d "$MD_DIR" ]; then
    echo "❌ Directory not found: $MD_DIR"
    exit 1
fi

if [ ! -f "$MD_DIR/md.cpt" ]; then
    echo "❌ Checkpoint file not found: $MD_DIR/md.cpt"
    echo "   Cannot resume simulation without checkpoint."
    exit 1
fi

echo "🔄 Restarting MD simulation for: $RESCUE"
echo "   Directory: $MD_DIR"
echo "   Resuming from checkpoint: md.cpt"
echo ""

# Check for running GROMACS processes
RUNNING_PID=$(ps aux | grep "gmx mdrun" | grep -v grep | grep "$RESCUE" | awk '{print $2}' || echo "")

if [ ! -z "$RUNNING_PID" ]; then
    echo "⚠️  Found running mdrun process (PID: $RUNNING_PID)"
    echo "   Stopping it now..."
    kill $RUNNING_PID
    sleep 2
    echo "✅ Stopped previous run"
    echo ""
fi

# Get number of available CPU cores
NCORES=$(sysctl -n hw.ncpu 2>/dev/null || echo "8")
echo "💻 Detected $NCORES CPU cores"
echo ""

# Change to MD directory
cd "$MD_DIR"

# Start optimized mdrun
echo "🚀 Starting optimized mdrun with $NCORES threads..."
echo "   Command: gmx mdrun -v -deffnm md -cpi md.cpt -ntomp $NCORES"
echo ""

# Run in background and capture PID
nohup gmx mdrun -v -deffnm md -cpi md.cpt -ntomp $NCORES > mdrun_restart.log 2>&1 &
NEW_PID=$!

sleep 2

# Check if process started successfully
if ps -p $NEW_PID > /dev/null; then
    echo "✅ MD simulation restarted successfully!"
    echo "   PID: $NEW_PID"
    echo "   Using: $NCORES CPU threads"
    echo ""
    echo "📊 Monitor progress:"
    echo "   tail -f $MD_DIR/md.log"
    echo "   tail -f $MD_DIR/mdrun_restart.log"
    echo ""
    echo "⏱️  Expected performance: ~4.5-5.0 ns/day"
    echo "   10 ns completion: ~2 days from now"
else
    echo "❌ Failed to start mdrun"
    echo "   Check log: $MD_DIR/mdrun_restart.log"
    exit 1
fi
