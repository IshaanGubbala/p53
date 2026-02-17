#!/usr/bin/env python3
"""
Real-time campaign progress monitor with beautiful table display.

Usage:
    python scripts/monitor_campaign.py

Or run in a separate terminal while the campaign is running.
"""

import os
import sys
import time
from pathlib import Path
from datetime import datetime, timedelta
import re

# System monitoring
try:
    import psutil
    import pynvml
    MONITORING_AVAILABLE = True
except ImportError:
    MONITORING_AVAILABLE = False
    psutil = None
    pynvml = None

# ANSI color codes
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def clear_screen():
    """Clear terminal screen."""
    os.system('cls' if os.name == 'nt' else 'clear')

def get_system_resources():
    """Get current system resource usage."""
    resources = {
        'cpu_percent': 0,
        'ram_used_gb': 0,
        'ram_total_gb': 0,
        'gpu_util': 0,
        'gpu_mem_used_gb': 0,
        'gpu_mem_total_gb': 0,
        'gpu_mem_shared_gb': 0,
        'available': False
    }

    if not MONITORING_AVAILABLE:
        return resources

    try:
        # CPU and RAM
        resources['cpu_percent'] = psutil.cpu_percent(interval=0.1)
        mem = psutil.virtual_memory()
        resources['ram_used_gb'] = mem.used / (1024**3)
        resources['ram_total_gb'] = mem.total / (1024**3)

        # GPU (NVIDIA only)
        try:
            pynvml.nvmlInit()
            handle = pynvml.nvmlDeviceGetHandleByIndex(0)

            # GPU utilization
            util = pynvml.nvmlDeviceGetUtilizationRates(handle)
            resources['gpu_util'] = util.gpu

            # GPU memory
            mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
            resources['gpu_mem_used_gb'] = mem_info.used / (1024**3)
            resources['gpu_mem_total_gb'] = mem_info.total / (1024**3)

            pynvml.nvmlShutdown()
            resources['available'] = True
        except:
            pass  # GPU monitoring failed, skip

    except:
        pass  # Monitoring failed, return zeros

    return resources

def parse_log_file(log_path):
    """Parse the log file for progress information."""
    if not log_path.exists():
        return None

    try:
        with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
    except:
        return None

    # Only process lines from the most recent campaign start
    # Find the last "Pass A [1/" marker to ignore old campaign data
    last_start_idx = 0
    for i, line in enumerate(lines):
        if 'Pass A [1/' in line or ('engine.campaign' in line and 'Pass A [1/' in line):
            last_start_idx = i
    lines = lines[last_start_idx:]

    info = {
        'status': 'Initializing',
        'current_step': 'Starting',
        'variants_tested': 0,
        'best_score': None,
        'total_candidates': 0,
        'scenarios_completed': 0,
        'total_scenarios': 0,
        'validation_rank': 0,
        'validation_target': None,
        'validation_timestamps': [],  # List of (rank, timestamp) tuples
        'validation_start_time': None,
        'in_validation': False,
        'physics_validation_count': 0,  # Count of completed energy minimizations
        'physics_timestamps': [],  # List of (count, timestamp) tuples
        'in_physics_validation': False,
        'current_receptor': None,
        'elapsed_time': None,
        'elapsed_seconds': 0,
        'last_update': None,
        'errors': [],
        'warnings': [],
    }

    start_time = None
    campaign_started = False

    for line in lines:
        # Ignore doctor command output - only look at actual campaign runs
        if 'doctor' in line.lower() or 'diagnostic' in line.lower():
            continue

        # Mark when campaign actually starts
        if 'campaign' in line.lower() and 'run' in line.lower():
            campaign_started = True
        # Extract timestamp (only from campaign activity)
        timestamp_match = re.match(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', line)
        if timestamp_match and campaign_started:
            if not start_time:
                start_time = datetime.strptime(timestamp_match.group(1), '%Y-%m-%d %H:%M:%S')
            info['last_update'] = timestamp_match.group(1)

        # Extract status information (only from campaign activity)
        # Pass A/B indicate campaign has started
        if 'Pass A' in line or ('screen' in line.lower() and 'scenario' in line.lower()):
            info['current_step'] = 'Pass A: Screening variants'
            info['status'] = 'Running'
            campaign_started = True
            # Extract total scenarios from "Pass A [1/108]" format
            scenario_match = re.search(r'Pass A \[(\d+)/(\d+)\]', line)
            if scenario_match:
                info['total_scenarios'] = int(scenario_match.group(2))
        elif 'Pass B' in line or ('optimize' in line.lower() and 'scenario' in line.lower()):
            info['current_step'] = 'Pass B: Optimizing candidates'
            info['status'] = 'Running'
            campaign_started = True
        elif 'Energy minimization:' in line:
            # Physics validation phase - count actual completions
            info['current_step'] = 'Physics Validation (Energy Minimization)'
            info['status'] = 'Running'
            info['in_physics_validation'] = True
            info['in_validation'] = True
            campaign_started = True
            timestamp_match = re.match(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', line)

            # Increment count for each completed energy minimization
            info['physics_validation_count'] += 1
            if timestamp_match:
                ts = timestamp_match.group(1)
                info['physics_timestamps'].append((info['physics_validation_count'], ts))
                if not info['validation_start_time']:
                    info['validation_start_time'] = ts
        elif 'Tier 1 validation' in line and not info['in_physics_validation']:
            # ESMFold validation phase (before physics)
            info['current_step'] = 'ESMFold Structural Validation'
            info['status'] = 'Running'
            info['in_validation'] = True
            campaign_started = True
            # Extract rank and target: "Tier 1 validation: rank=2 target=R273C+R282W"
            rank_match = re.search(r'rank=(\d+)', line)
            target_match = re.search(r'target=(\S+)', line)
            timestamp_match = re.match(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', line)

            if rank_match:
                rank = int(rank_match.group(1))
                info['validation_rank'] = rank
                # Track timestamp for this validation rank
                if timestamp_match:
                    ts = timestamp_match.group(1)
                    # Only add if this rank hasn't been tracked yet
                    if not any(r == rank for r, _ in info['validation_timestamps']):
                        info['validation_timestamps'].append((rank, ts))
            if target_match:
                info['validation_target'] = target_match.group(1)
            # Track when validation phase started
            if not info['validation_start_time'] and timestamp_match:
                info['validation_start_time'] = timestamp_match.group(1)
        elif 'Drug generation' in line and campaign_started:
            info['current_step'] = 'Generating drug molecules'
            info['status'] = 'Running'
        elif 'Docking' in line and campaign_started:
            info['current_step'] = 'Molecular docking'
            info['status'] = 'Running'
        elif 'COMPLETE' in line and campaign_started:
            info['current_step'] = 'Campaign finished'
            info['status'] = 'Complete'
        elif 'ERROR' in line and campaign_started:
            info['errors'].append(line.strip())
            info['status'] = 'Error'
        elif 'WARNING' in line and campaign_started:
            info['warnings'].append(line.strip())

        # Extract metrics
        # Track completed scenarios for accurate totals
        if 'Scenario' in line and 'complete:' in line:
            info['scenarios_completed'] += 1
            cand_match = re.search(r'(\d+)\s*candidates', line)
            if cand_match:
                # Accumulate candidates from completed scenarios
                info['total_candidates'] += int(cand_match.group(1))
            score_match = re.search(r'best_score=([-+]?\d*\.?\d+)', line)
            if score_match:
                try:
                    score = float(score_match.group(1))
                    if info['best_score'] is None or score > info['best_score']:
                        info['best_score'] = score
                except:
                    pass

        # Progress messages: "Progress: 5/15 trials, 5 candidates, best=0.723"
        # These show CURRENT scenario progress, don't use for total
        if 'Progress:' in line and 'trials' in line:
            trial_match = re.search(r'(\d+)/(\d+)\s*trials', line)
            if trial_match:
                info['variants_tested'] = int(trial_match.group(1))
            score_match = re.search(r'best=([-+]?\d*\.?\d+)', line)
            if score_match:
                try:
                    score = float(score_match.group(1))
                    if info['best_score'] is None or score > info['best_score']:
                        info['best_score'] = score
                except:
                    pass

        # Legacy/fallback metrics
        if 'variants tested' in line.lower() and 'Progress:' not in line:
            match = re.search(r'(\d+)\s*variants', line)
            if match:
                info['variants_tested'] = int(match.group(1))

        if 'receptor' in line.lower():
            match = re.search(r'receptor[:\s]+(\w+)', line, re.IGNORECASE)
            if match:
                info['current_receptor'] = match.group(1)

        # Only parse scores from docking/affinity lines, not from generic "score" mentions
        if ('affinity' in line.lower() or 'binding' in line.lower()) and 'kcal' in line.lower():
            match = re.search(r'[-+]?\d*\.\d+', line)
            if match:
                try:
                    score = float(match.group(0))
                    # Binding scores are typically negative, filter out years/dates
                    if -50 < score < 50:
                        if info['best_score'] is None or score < info['best_score']:
                            info['best_score'] = score
                except:
                    pass

    # Calculate elapsed time
    if start_time and campaign_started:
        elapsed = datetime.now() - start_time
        info['elapsed_time'] = str(elapsed).split('.')[0]
        info['elapsed_seconds'] = elapsed.total_seconds()

    # If campaign hasn't started, return None so we show waiting message
    if not campaign_started:
        return None

    return info

def print_header():
    """Print monitor header."""
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.ENDC}")
    print(f"{Colors.BOLD}{Colors.CYAN}  p53CAD Campaign Monitor{Colors.ENDC}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.ENDC}\n")

def print_status_table(info):
    """Print status information in a beautiful table."""
    if not info:
        print(f"{Colors.YELLOW}⏳ Waiting for campaign to start...{Colors.ENDC}\n")
        return

    # Status color
    status_color = Colors.GREEN if info['status'] == 'Running' else Colors.RED if info['status'] == 'Error' else Colors.YELLOW

    # Print main status box
    print(f"{Colors.BOLD}┌{'─'*78}┐{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Status:{Colors.ENDC} {status_color}{info['status']:<69}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Current Step:{Colors.ENDC} {Colors.CYAN}{info['current_step']:<62}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*78}┤{Colors.ENDC}")

    # Metrics table
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Metric{Colors.ENDC}                    │ {Colors.BOLD}Value{Colors.ENDC}                                     {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*27}┼{'─'*50}┤{Colors.ENDC}")

    # Variants tested
    variants_str = f"{info['variants_tested']:,}"
    print(f"{Colors.BOLD}│{Colors.ENDC} Variants Tested        │ {Colors.GREEN}{variants_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # Total candidates
    candidates_str = f"{info['total_candidates']:,}"
    print(f"{Colors.BOLD}│{Colors.ENDC} Total Candidates       │ {Colors.GREEN}{candidates_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # Best score (oracle score for protein rescue campaigns)
    if info['best_score'] is not None:
        # Check if it's a binding affinity (negative, kcal/mol) or oracle score (0-1)
        if info['best_score'] < 0:
            score_str = f"{info['best_score']:.2f} kcal/mol"
            score_color = Colors.GREEN if info['best_score'] < -7 else Colors.YELLOW
            label = "Best Binding Score"
        else:
            score_str = f"{info['best_score']:.3f}"
            score_color = Colors.GREEN if info['best_score'] > 0.5 else Colors.YELLOW
            label = "Best Oracle Score"
        print(f"{Colors.BOLD}│{Colors.ENDC} {label:<22} │ {score_color}{score_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    else:
        print(f"{Colors.BOLD}│{Colors.ENDC} Best Oracle Score      │ {Colors.YELLOW}{'N/A':<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # Current receptor
    if info['current_receptor']:
        print(f"{Colors.BOLD}│{Colors.ENDC} Current Receptor       │ {Colors.CYAN}{info['current_receptor']:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # Elapsed time
    if info['elapsed_time']:
        print(f"{Colors.BOLD}│{Colors.ENDC} Elapsed Time           │ {Colors.BLUE}{info['elapsed_time']:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # Last update
    if info['last_update']:
        print(f"{Colors.BOLD}│{Colors.ENDC} Last Update            │ {Colors.BLUE}{info['last_update']:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # Validation phase tracking
    if info['in_validation']:
        # Determine which validation phase we're in
        if info['in_physics_validation'] and info['physics_validation_count'] > 0:
            # Physics validation - show energy minimization progress
            progress_pct = (info['physics_validation_count'] / 30) * 100
            validation_str = f"Structure {info['physics_validation_count']} of 30 ({progress_pct:.1f}%)"
            print(f"{Colors.BOLD}│{Colors.ENDC} Physics Validation     │ {Colors.CYAN}{validation_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
        elif info['validation_rank'] > 0:
            # ESMFold validation
            if info['validation_target']:
                target_display = info['validation_target'][:48]
                print(f"{Colors.BOLD}│{Colors.ENDC} Validating Target      │ {Colors.CYAN}{target_display:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

            progress_pct = (info['validation_rank'] / 30) * 100
            validation_str = f"Rank {info['validation_rank']} of 30 ({progress_pct:.1f}%)"
            print(f"{Colors.BOLD}│{Colors.ENDC} Validation Progress    │ {Colors.CYAN}{validation_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

        # Show elapsed time for current structure
        timestamps_to_use = info['physics_timestamps'] if info['in_physics_validation'] else info['validation_timestamps']
        current_count = info['physics_validation_count'] if info['in_physics_validation'] else info['validation_rank']

        if len(timestamps_to_use) >= 1:
            current_ts = timestamps_to_use[-1][1]  # Most recent timestamp
            try:
                current_dt = datetime.strptime(current_ts, '%Y-%m-%d %H:%M:%S')
                elapsed_current = datetime.now() - current_dt
                elapsed_mins, elapsed_secs = divmod(int(elapsed_current.total_seconds()), 60)
                label = "current structure" if info['in_physics_validation'] else "current rank"
                elapsed_str = f"{elapsed_mins}m {elapsed_secs}s ago ({label})"
                print(f"{Colors.BOLD}│{Colors.ENDC} Current Structure Time │ {Colors.BLUE}{elapsed_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
            except:
                pass

        # Calculate time per validation from tracked timestamps
        if len(timestamps_to_use) >= 2:
            # Sort timestamps by count/rank
            sorted_timestamps = sorted(timestamps_to_use, key=lambda x: x[0])
            # Calculate time differences between consecutive validations
            time_diffs = []
            for i in range(1, len(sorted_timestamps)):
                prev_count, prev_ts = sorted_timestamps[i-1]
                curr_count, curr_ts = sorted_timestamps[i]
                try:
                    prev_dt = datetime.strptime(prev_ts, '%Y-%m-%d %H:%M:%S')
                    curr_dt = datetime.strptime(curr_ts, '%Y-%m-%d %H:%M:%S')
                    diff_seconds = (curr_dt - prev_dt).total_seconds()
                    time_diffs.append(diff_seconds)
                except:
                    pass

            if time_diffs:
                # Show last validation time
                last_time_seconds = time_diffs[-1]
                last_mins, last_secs = divmod(int(last_time_seconds), 60)
                last_time_str = f"{last_mins}m {last_secs}s (last validation)"
                print(f"{Colors.BOLD}│{Colors.ENDC} Time per Validation    │ {Colors.BLUE}{last_time_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

                # Calculate average time per validation
                avg_time_seconds = sum(time_diffs) / len(time_diffs)
                avg_mins, avg_secs = divmod(int(avg_time_seconds), 60)
                avg_time_str = f"{avg_mins}m {avg_secs}s (average)"
                print(f"{Colors.BOLD}│{Colors.ENDC} Avg Time per Validation│ {Colors.BLUE}{avg_time_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

                # Calculate ETA based on actual average time
                remaining_count = 30 - current_count
                eta_seconds = remaining_count * avg_time_seconds
                eta_hours, remainder = divmod(int(eta_seconds), 3600)
                eta_mins, _ = divmod(remainder, 60)
                if eta_hours > 0:
                    eta_str = f"~{eta_hours}h {eta_mins}m remaining (based on avg)"
                else:
                    eta_str = f"~{eta_mins}m remaining (based on avg)"
                print(f"{Colors.BOLD}│{Colors.ENDC} Estimated Time Left    │ {Colors.GREEN}{eta_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
            else:
                # Fallback to rough estimate
                remaining_count = 30 - current_count
                default_time = 12 if info['in_physics_validation'] else 20  # Physics is faster
                eta_seconds = remaining_count * default_time * 60
                eta_hours, remainder = divmod(int(eta_seconds), 3600)
                eta_mins, _ = divmod(remainder, 60)
                if eta_hours > 0:
                    eta_str = f"~{eta_hours}h {eta_mins}m remaining (estimated)"
                else:
                    eta_str = f"~{eta_mins}m remaining (estimated)"
                print(f"{Colors.BOLD}│{Colors.ENDC} Estimated Time Left    │ {Colors.YELLOW}{eta_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
        else:
            # Not enough data yet, use fallback estimate
            remaining_count = 30 - current_count
            default_time = 12 if info['in_physics_validation'] else 20  # Physics is faster
            eta_seconds = remaining_count * default_time * 60
            eta_hours, remainder = divmod(int(eta_seconds), 3600)
            eta_mins, _ = divmod(remainder, 60)
            if eta_hours > 0:
                eta_str = f"~{eta_hours}h {eta_mins}m remaining (estimated)"
            else:
                eta_str = f"~{eta_mins}m remaining (estimated)"
            print(f"{Colors.BOLD}│{Colors.ENDC} Estimated Time Left    │ {Colors.YELLOW}{eta_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # Standard scenario tracking (for Pass A/B phases)
    elif info['scenarios_completed'] > 0 and info['elapsed_seconds'] > 0:
        # Scenarios progress
        scenarios_str = f"{info['scenarios_completed']}/{info['total_scenarios']}" if info['total_scenarios'] > 0 else f"{info['scenarios_completed']}"
        print(f"{Colors.BOLD}│{Colors.ENDC} Scenarios Progress     │ {Colors.CYAN}{scenarios_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

        # Time per scenario
        time_per_scenario = info['elapsed_seconds'] / info['scenarios_completed']
        mins, secs = divmod(int(time_per_scenario), 60)
        time_per_scenario_str = f"{mins}m {secs}s per scenario"
        print(f"{Colors.BOLD}│{Colors.ENDC} Avg Time per Scenario  │ {Colors.BLUE}{time_per_scenario_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

        # Time per candidate
        if info['total_candidates'] > 0:
            time_per_candidate = info['elapsed_seconds'] / info['total_candidates']
            mins, secs = divmod(int(time_per_candidate), 60)
            time_per_candidate_str = f"{mins}m {secs}s per candidate"
            print(f"{Colors.BOLD}│{Colors.ENDC} Avg Time per Candidate │ {Colors.BLUE}{time_per_candidate_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

        # ETA (only if not in validation)
        if info['total_scenarios'] > 0 and not info['in_validation']:
            remaining_scenarios = info['total_scenarios'] - info['scenarios_completed']
            eta_seconds = remaining_scenarios * time_per_scenario
            eta_hours, remainder = divmod(int(eta_seconds), 3600)
            eta_mins, eta_secs = divmod(remainder, 60)
            if eta_hours > 0:
                eta_str = f"{eta_hours}h {eta_mins}m remaining"
            else:
                eta_str = f"{eta_mins}m {eta_secs}s remaining"
            print(f"{Colors.BOLD}│{Colors.ENDC} Estimated Time Left    │ {Colors.GREEN}{eta_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # System resources section
    resources = get_system_resources()
    if resources['available'] or MONITORING_AVAILABLE:
        print(f"{Colors.BOLD}├{'─'*78}┤{Colors.ENDC}")
        print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}System Resources{Colors.ENDC}                                                         {Colors.BOLD}│{Colors.ENDC}")
        print(f"{Colors.BOLD}├{'─'*78}┤{Colors.ENDC}")

        # CPU usage
        cpu_str = f"{resources['cpu_percent']:.1f}%"
        cpu_color = Colors.GREEN if resources['cpu_percent'] < 70 else Colors.YELLOW if resources['cpu_percent'] < 90 else Colors.RED
        print(f"{Colors.BOLD}│{Colors.ENDC} CPU Usage              │ {cpu_color}{cpu_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

        # System RAM
        ram_str = f"{resources['ram_used_gb']:.1f} / {resources['ram_total_gb']:.1f} GB ({resources['ram_used_gb']/resources['ram_total_gb']*100:.0f}%)"
        ram_pct = resources['ram_used_gb']/resources['ram_total_gb']*100 if resources['ram_total_gb'] > 0 else 0
        ram_color = Colors.GREEN if ram_pct < 70 else Colors.YELLOW if ram_pct < 90 else Colors.RED
        print(f"{Colors.BOLD}│{Colors.ENDC} System RAM             │ {ram_color}{ram_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

        # GPU usage
        if resources['gpu_util'] > 0 or resources['gpu_mem_total_gb'] > 0:
            gpu_str = f"{resources['gpu_util']:.1f}%"
            gpu_color = Colors.GREEN if resources['gpu_util'] > 40 else Colors.YELLOW if resources['gpu_util'] > 10 else Colors.CYAN
            print(f"{Colors.BOLD}│{Colors.ENDC} GPU Usage              │ {gpu_color}{gpu_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

            # GPU VRAM
            vram_str = f"{resources['gpu_mem_used_gb']:.1f} / {resources['gpu_mem_total_gb']:.1f} GB ({resources['gpu_mem_used_gb']/resources['gpu_mem_total_gb']*100:.0f}%)"
            vram_pct = resources['gpu_mem_used_gb']/resources['gpu_mem_total_gb']*100 if resources['gpu_mem_total_gb'] > 0 else 0
            vram_color = Colors.GREEN if vram_pct < 80 else Colors.YELLOW if vram_pct < 95 else Colors.RED
            print(f"{Colors.BOLD}│{Colors.ENDC} GPU VRAM               │ {vram_color}{vram_str:<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    print(f"{Colors.BOLD}└{'─'*78}┘{Colors.ENDC}\n")

    # Warnings and errors
    if info['warnings']:
        print(f"{Colors.YELLOW}⚠️  Warnings: {len(info['warnings'])}{Colors.ENDC}")
        for warning in info['warnings'][-3:]:  # Show last 3
            print(f"   {Colors.YELLOW}•{Colors.ENDC} {warning[:75]}")
        print()

    if info['errors']:
        print(f"{Colors.RED}❌ Errors: {len(info['errors'])}{Colors.ENDC}")
        for error in info['errors'][-3:]:  # Show last 3
            print(f"   {Colors.RED}•{Colors.ENDC} {error[:75]}")
        print()

def monitor_campaign(log_path, refresh_interval=2):
    """Monitor campaign progress in real-time."""
    print(f"{Colors.CYAN}Monitoring log file: {log_path}{Colors.ENDC}\n")

    try:
        while True:
            clear_screen()
            print_header()

            info = parse_log_file(log_path)
            print_status_table(info)

            print(f"{Colors.BOLD}Press Ctrl+C to exit{Colors.ENDC}")
            print(f"{Colors.BLUE}Refreshing every {refresh_interval} seconds...{Colors.ENDC}")

            time.sleep(refresh_interval)

    except KeyboardInterrupt:
        print(f"\n\n{Colors.GREEN}✓ Monitoring stopped{Colors.ENDC}\n")

def main():
    """Main entry point."""
    project_root = Path(__file__).parent.parent
    log_path = project_root / "logs" / "p53cad_workflow.log"

    if not log_path.parent.exists():
        log_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"{Colors.BOLD}{Colors.GREEN}p53CAD Campaign Monitor{Colors.ENDC}\n")
    print(f"Log file: {log_path}")

    if not log_path.exists():
        print(f"{Colors.YELLOW}Log file not found yet. Waiting for campaign to start...{Colors.ENDC}\n")

    monitor_campaign(log_path, refresh_interval=2)

if __name__ == "__main__":
    main()
