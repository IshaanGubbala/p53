#!/usr/bin/env python3
"""Monitor campaign until completion, then notify."""

import time
import os
import sys
from datetime import datetime
import re

LOG_FILE = r"c:\Users\rajee\OneDrive\Documents\p53\logs\p53cad_workflow.log"
REFRESH_INTERVAL = 30  # seconds

def count_energy_minimizations():
    """Count completed energy minimizations."""
    if not os.path.exists(LOG_FILE):
        return 0

    with open(LOG_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
        return content.count("Energy minimization:")

def check_completion():
    """Check if campaign is complete."""
    if not os.path.exists(LOG_FILE):
        return False

    with open(LOG_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
        return "CAMPAIGN COMPLETE" in content or "Campaign finished" in content

def get_latest_timestamp():
    """Get timestamp of last log entry."""
    if not os.path.exists(LOG_FILE):
        return None

    with open(LOG_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        for line in reversed(lines):
            match = re.match(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', line)
            if match:
                return match.group(1)
    return None

def format_time(seconds):
    """Format seconds into readable time."""
    hours, remainder = divmod(int(seconds), 3600)
    minutes, secs = divmod(remainder, 60)
    if hours > 0:
        return f"{hours}h {minutes}m"
    else:
        return f"{minutes}m {secs}s"

def main():
    print("=" * 80)
    print(" p53CAD Campaign Monitor - Auto-Complete")
    print("=" * 80)
    print()
    print("Monitoring campaign until completion...")
    print(f"Log file: {LOG_FILE}")
    print(f"Refresh interval: {REFRESH_INTERVAL}s")
    print()
    print("Press Ctrl+C to stop monitoring")
    print("=" * 80)
    print()

    start_time = time.time()
    last_count = 0

    try:
        while True:
            # Check completion
            if check_completion():
                print("\n" + "=" * 80)
                print(" 🎉 CAMPAIGN COMPLETE! 🎉")
                print("=" * 80)
                elapsed = time.time() - start_time
                print(f"\nMonitoring time: {format_time(elapsed)}")
                print("\nResults available at:")
                print("  - Streamlit: http://localhost:8501")
                print("  - Campaign dir: data/campaigns/campaign_20260216_090733/")
                print("\nTo enable CUDA for future runs:")
                print("  Run: setup_openmm_cuda.bat")
                print("=" * 80)
                break

            # Count progress
            count = count_energy_minimizations()
            latest = get_latest_timestamp()
            elapsed = time.time() - start_time

            # Display status
            now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            progress_pct = (count / 30) * 100 if count > 0 else 0

            print(f"[{now}] Physics validation: {count}/30 structures ({progress_pct:.1f}%)", end="")

            if count > last_count:
                print(" ✓ NEW!", end="")

            if latest:
                print(f" | Last update: {latest}", end="")

            print(f" | Monitoring: {format_time(elapsed)}")

            last_count = count

            # Wait
            time.sleep(REFRESH_INTERVAL)

    except KeyboardInterrupt:
        print("\n\nMonitoring stopped by user.")
        print(f"Progress: {last_count}/30 structures completed")
        print("Campaign still running in background.")

    except Exception as e:
        print(f"\n\nError: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
