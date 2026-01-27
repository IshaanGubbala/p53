#!/bin/bash
# Check MD simulation performance
# Usage: ./check_performance.sh <rescue_name>

RESCUE=${1:-"A189S_M133L_S95T"}
MD_DIR="Data/processed/md_simulations/${RESCUE}"

if [ ! -d "$MD_DIR" ]; then
    echo "❌ Directory not found: $MD_DIR"
    exit 1
fi

echo "📊 MD Simulation Performance Check"
echo "==================================="
echo ""

# Check if running
PID=$(ps aux | grep "gmx mdrun" | grep -v grep | grep "$RESCUE" | awk '{print $2}' || echo "")
if [ -z "$PID" ]; then
    echo "⚠️  No mdrun process found for $RESCUE"
    echo "   Simulation may be complete or crashed"
    echo ""
else
    CPU=$(ps aux | grep "$PID" | grep -v grep | awk '{print $3}')
    echo "✅ Process running"
    echo "   PID:       $PID"
    echo "   CPU Usage: ${CPU}%"
    echo ""
fi

# Get latest progress from log
cd "$MD_DIR"

if [ ! -f "md.log" ]; then
    echo "❌ Log file not found: md.log"
    exit 1
fi

# Extract latest step and time
LATEST_STEP=$(grep "Step" md.log | grep -v "step" | tail -1 | awk '{print $2}')
LATEST_TIME=$(grep "Time" md.log | grep -v "Time:" | tail -1 | awk '{print $2}')

# Extract earlier checkpoint for comparison (10 entries back)
EARLIER_STEP=$(grep "Step" md.log | grep -v "step" | tail -11 | head -1 | awk '{print $2}')
EARLIER_TIME=$(grep "Time" md.log | grep -v "Time:" | tail -11 | head -1 | awk '{print $2}')

echo "📈 Current Progress:"
echo "   Step: $(printf "%'d" $LATEST_STEP) / 5,000,000"
echo "   Time: ${LATEST_TIME} ps / 10,000 ps"
echo "   $(echo "scale=1; $LATEST_TIME / 100" | bc)% complete"
echo ""

# Get file sizes
XTC_SIZE=$(ls -lh md.xtc 2>/dev/null | awk '{print $5}')
LOG_SIZE=$(ls -lh md.log 2>/dev/null | awk '{print $5}')

echo "📁 Output Files:"
echo "   Trajectory: $XTC_SIZE (md.xtc)"
echo "   Log:        $LOG_SIZE (md.log)"
echo ""

# Calculate performance (if we have data)
if [ ! -z "$EARLIER_STEP" ] && [ ! -z "$EARLIER_TIME" ]; then
    python3 << EOF
import sys
from datetime import datetime, timedelta

latest_step = $LATEST_STEP
latest_time = $LATEST_TIME
earlier_step = $EARLIER_STEP if $EARLIER_STEP > 0 else 68150
earlier_time = $EARLIER_TIME if $EARLIER_TIME > 0 else 136.3

step_diff = latest_step - earlier_step
time_diff_ps = latest_time - earlier_time

if time_diff_ps > 0:
    # Assume ~5 seconds per log entry (5000 steps = 10 ps, logged every 10 ps)
    # So 10 entries back = ~50 seconds
    wall_time_seconds = 50  # Rough estimate

    ps_per_second = time_diff_ps / wall_time_seconds
    ns_per_day = (ps_per_second * 86400) / 1000
    hours_per_ns = 1 / (ps_per_second / 1000) if ps_per_second > 0 else 999

    print("⚡ Performance (rough estimate):")
    print(f"   {ns_per_day:.2f} ns/day")
    print(f"   {hours_per_ns:.2f} hours/ns")
    print("")

    # Estimate completion
    remaining_ps = 10000 - latest_time
    hours_remaining = remaining_ps / (ps_per_second / 3600) if ps_per_second > 0 else 9999
    days_remaining = hours_remaining / 24

    completion = datetime.now() + timedelta(hours=hours_remaining)

    print("⏰ Estimated Completion:")
    print(f"   {days_remaining:.1f} days remaining")
    print(f"   {completion.strftime('%A, %B %d at %I:%M %p')}")
    print("")

    print("ℹ️  Note: This is a rough estimate.")
    print("   For accurate measurement, wait 15+ min after restart")
else:
    print("⚠️  Not enough data for performance calculation")
    print("   Check again in a few minutes")
EOF
else
    echo "⚠️  No comparison data available"
    echo "   Log may be too short or simulation just started"
fi

echo ""
echo "🔄 Run this script again in 10-15 minutes for stable performance"
