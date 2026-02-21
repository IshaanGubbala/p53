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
        'available': False
    }

    if not MONITORING_AVAILABLE:
        return resources

    try:
        resources['cpu_percent'] = psutil.cpu_percent(interval=0.1)
        mem = psutil.virtual_memory()
        resources['ram_used_gb'] = mem.used / (1024**3)
        resources['ram_total_gb'] = mem.total / (1024**3)

        try:
            pynvml.nvmlInit()
            handle = pynvml.nvmlDeviceGetHandleByIndex(0)
            util = pynvml.nvmlDeviceGetUtilizationRates(handle)
            resources['gpu_util'] = util.gpu
            mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
            resources['gpu_mem_used_gb'] = mem_info.used / (1024**3)
            resources['gpu_mem_total_gb'] = mem_info.total / (1024**3)
            pynvml.nvmlShutdown()
            resources['available'] = True
        except:
            pass

    except:
        pass

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
    last_start_idx = 0
    for i, line in enumerate(lines):
        if 'Pass A [1/' in line:
            last_start_idx = i
    lines = lines[last_start_idx:]

    info = {
        'status': 'Initializing',
        # Phase detection: 'pass_a', 'pass_b', 'validation', 'physics', 'complete'
        'current_phase': 'pass_a',
        'current_step': 'Starting',

        # Pass A tracking
        'pass_a_total': 0,
        'pass_a_completed': 0,

        # Pass B tracking
        'pass_b_total': 0,
        'pass_b_completed': 0,
        'pass_b_current_idx': 0,   # current [X/Y] from log

        # Scenario start timestamps for real-time ETA interpolation
        'last_a_scenario_start_dt': None,
        'last_b_scenario_start_dt': None,

        # Current scenario trial progress
        'current_trial': 0,
        'current_trial_total': 0,

        # Candidates / scores
        'total_candidates': 0,
        'best_score': None,

        # Validation (ESMFold + Physics)
        'validation_rank': 0,
        'validation_target': None,
        'validation_timestamps': [],
        'validation_start_time': None,
        'in_validation': False,
        'physics_validation_count': 0,
        'physics_timestamps': [],
        'in_physics_validation': False,

        # Timing
        'elapsed_time': None,
        'elapsed_seconds': 0,
        'last_update': None,

        'errors': [],
        'warnings': [],
    }

    start_time = None
    campaign_started = False
    # Track timing per pass for ETA
    pass_a_times = []   # seconds between consecutive Pass A scenario completions
    pass_b_times = []
    last_scenario_ts = None

    for line in lines:
        if 'doctor' in line.lower() or 'diagnostic' in line.lower():
            continue

        # Timestamp extraction
        ts_match = re.match(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', line)
        ts_str = ts_match.group(1) if ts_match else None
        ts_dt = None
        if ts_str:
            try:
                ts_dt = datetime.strptime(ts_str, '%Y-%m-%d %H:%M:%S')
            except:
                pass

        if ts_str and campaign_started:
            if not start_time and ts_dt:
                start_time = ts_dt
            info['last_update'] = ts_str

        # ── Phase detection ───────────────────────────────────────────────────

        # Pass A line: "Pass A [idx/total] scenario_id"
        pass_a_match = re.search(r'Pass A \[(\d+)/(\d+)\]', line)
        if pass_a_match:
            info['current_phase'] = 'pass_a'
            info['current_step'] = 'Pass A: Screening all scenarios'
            info['status'] = 'Running'
            campaign_started = True
            if not start_time and ts_dt:
                start_time = ts_dt
            a_idx = int(pass_a_match.group(1))
            a_total = int(pass_a_match.group(2))
            info['pass_a_total'] = a_total
            # current_idx is ahead of completed by 1 (scenario just started)
            info['pass_a_completed'] = max(info['pass_a_completed'], a_idx - 1)
            if ts_dt:
                info['last_a_scenario_start_dt'] = ts_dt

        # Pass B line: "Pass B [idx/total] scenario_id"
        pass_b_match = re.search(r'Pass B \[(\d+)/(\d+)\]', line)
        if pass_b_match:
            info['current_phase'] = 'pass_b'
            info['current_step'] = 'Pass B: Deep optimizing top scenarios'
            info['status'] = 'Running'
            campaign_started = True
            b_idx = int(pass_b_match.group(1))
            b_total = int(pass_b_match.group(2))
            info['pass_b_total'] = b_total
            info['pass_b_current_idx'] = b_idx
            info['pass_b_completed'] = max(info['pass_b_completed'], b_idx - 1)
            if ts_dt:
                info['last_b_scenario_start_dt'] = ts_dt

        # Trial progress: "Progress: 5/6 trials, 5 candidates, best=0.723"
        if 'Progress:' in line and 'trials' in line:
            trial_match = re.search(r'(\d+)/(\d+)\s*trials', line)
            if trial_match:
                info['current_trial'] = int(trial_match.group(1))
                info['current_trial_total'] = int(trial_match.group(2))
            score_match = re.search(r'best=([-+]?\d*\.?\d+)', line)
            if score_match:
                try:
                    score = float(score_match.group(1))
                    if info['best_score'] is None or score > info['best_score']:
                        info['best_score'] = score
                except:
                    pass

        # Scenario complete: "Scenario X complete: N candidates, best_score=Y, baseline=Z"
        if 'Scenario' in line and 'complete:' in line:
            if ts_dt:
                if last_scenario_ts:
                    diff = (ts_dt - last_scenario_ts).total_seconds()
                    if info['current_phase'] == 'pass_a':
                        pass_a_times.append(diff)
                    else:
                        pass_b_times.append(diff)
                last_scenario_ts = ts_dt

            # Increment the correct pass counter
            if info['current_phase'] == 'pass_a':
                info['pass_a_completed'] += 1
            else:
                info['pass_b_completed'] += 1

            cand_match = re.search(r'(\d+)\s*candidates', line)
            if cand_match:
                info['total_candidates'] += int(cand_match.group(1))

            score_match = re.search(r'best_score=([-+]?\d*\.?\d+)', line)
            if score_match:
                try:
                    score = float(score_match.group(1))
                    if info['best_score'] is None or score > info['best_score']:
                        info['best_score'] = score
                except:
                    pass

        # Physics validation
        if 'Energy minimization:' in line:
            info['current_phase'] = 'physics'
            info['current_step'] = 'Physics Validation (Energy Minimization)'
            info['status'] = 'Running'
            info['in_physics_validation'] = True
            info['in_validation'] = True
            campaign_started = True
            info['physics_validation_count'] += 1
            if ts_str:
                info['physics_timestamps'].append((info['physics_validation_count'], ts_str))
                if not info['validation_start_time']:
                    info['validation_start_time'] = ts_str

        # ESMFold validation
        elif 'Tier 1 validation' in line and not info['in_physics_validation']:
            info['current_phase'] = 'validation'
            info['current_step'] = 'ESMFold Structural Validation'
            info['status'] = 'Running'
            info['in_validation'] = True
            campaign_started = True
            rank_match = re.search(r'rank=(\d+)', line)
            target_match = re.search(r'target=(\S+)', line)
            if rank_match:
                rank = int(rank_match.group(1))
                info['validation_rank'] = rank
                if ts_str and not any(r == rank for r, _ in info['validation_timestamps']):
                    info['validation_timestamps'].append((rank, ts_str))
            if target_match:
                info['validation_target'] = target_match.group(1)
            if not info['validation_start_time'] and ts_str:
                info['validation_start_time'] = ts_str

        # Campaign complete
        if 'CAMPAIGN COMPLETE' in line and campaign_started:
            info['current_phase'] = 'complete'
            info['current_step'] = 'Campaign finished'
            info['status'] = 'Complete'

        # Errors / warnings
        if 'ERROR' in line and campaign_started:
            info['errors'].append(line.strip())
            info['status'] = 'Error'
        elif 'WARNING' in line and campaign_started:
            info['warnings'].append(line.strip())

        # Binding affinity fallback score
        if ('affinity' in line.lower() or 'binding' in line.lower()) and 'kcal' in line.lower():
            match = re.search(r'[-+]?\d*\.\d+', line)
            if match:
                try:
                    score = float(match.group(0))
                    if -50 < score < 50:
                        if info['best_score'] is None or score < info['best_score']:
                            info['best_score'] = score
                except:
                    pass

    # Store timing lists for ETA
    info['pass_a_times'] = pass_a_times
    info['pass_b_times'] = pass_b_times

    # Elapsed time
    if start_time and campaign_started:
        elapsed = datetime.now() - start_time
        info['elapsed_time'] = str(elapsed).split('.')[0]
        info['elapsed_seconds'] = elapsed.total_seconds()

    if not campaign_started:
        return None

    return info


def _avg_secs(times: list, default_secs: int) -> float:
    """Return average duration or default."""
    return sum(times) / len(times) if times else float(default_secs)


def _fmt_duration(seconds: int) -> str:
    """Format seconds as Xh Ym Zs."""
    h, rem = divmod(max(seconds, 0), 3600)
    m, s = divmod(rem, 60)
    if h > 0:
        return f"{h}h {m}m"
    if m > 0:
        return f"{m}m {s}s"
    return f"{s}s"


def _eta_str(remaining: int, times: list, default_secs: int,
             elapsed_in_current: float = 0.0) -> str:
    """Format ETA string.

    remaining           -- number of scenarios left (including in-progress one)
    times               -- list of observed per-scenario durations (secs)
    default_secs        -- fallback avg when times is empty
    elapsed_in_current  -- seconds already spent on current in-progress scenario
    """
    avg = _avg_secs(times, default_secs)
    # Total raw ETA: remaining full scenarios * avg, minus time already spent
    eta = max(0, int(remaining * avg - elapsed_in_current))
    return f"{_fmt_duration(eta)} remaining ({remaining} scenarios left)"


def _progress_bar(done: int, total: int, width: int = 30) -> str:
    """Simple ASCII progress bar."""
    if total <= 0:
        return f"{'?':>{width}}"
    frac = min(done / total, 1.0)
    filled = int(frac * width)
    bar = '█' * filled + '░' * (width - filled)
    pct = frac * 100
    return f"{bar} {done}/{total} ({pct:.0f}%)"


def print_header():
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.ENDC}")
    print(f"{Colors.BOLD}{Colors.CYAN}  p53CAD Campaign Monitor{Colors.ENDC}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.ENDC}\n")


def print_status_table(info):
    if not info:
        print(f"{Colors.YELLOW}⏳ Waiting for campaign to start...{Colors.ENDC}\n")
        return

    phase = info['current_phase']
    status_color = (Colors.GREEN if info['status'] == 'Running'
                    else Colors.RED if info['status'] == 'Error'
                    else Colors.YELLOW)

    # ── Header box ───────────────────────────────────────────────────────────
    print(f"{Colors.BOLD}┌{'─'*78}┐{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Status:{Colors.ENDC} {status_color}{info['status']:<69}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Phase: {Colors.ENDC} {Colors.CYAN}{info['current_step']:<70}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*27}┼{'─'*50}┤{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Metric{Colors.ENDC}                    │ {Colors.BOLD}Value{Colors.ENDC}                                     {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*27}┼{'─'*50}┤{Colors.ENDC}")

    def row(label, value, color=Colors.GREEN):
        print(f"{Colors.BOLD}│{Colors.ENDC} {label:<22}    │ {color}{str(value):<48}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # ── Phase-specific metrics ───────────────────────────────────────────────

    if phase in ('pass_a', 'pass_b'):
        a_done = info['pass_a_completed']
        a_total = info['pass_a_total']
        b_done = info['pass_b_completed']
        b_total = info['pass_b_total']
        b_idx = info['pass_b_current_idx']

        if phase == 'pass_a':
            # Pass A progress bar
            bar = _progress_bar(a_done, a_total)
            row("Pass A Progress", bar, Colors.CYAN)

            # Current scenario trial progress
            if info['current_trial_total'] > 0:
                trial_str = f"{info['current_trial']}/{info['current_trial_total']} trials in current scenario"
                row("Current Scenario", trial_str, Colors.BLUE)

        else:  # pass_b
            # Pass A done (summary)
            a_summary = f"{a_done}/{a_total} done" if a_total > 0 else f"{a_done} done"
            row("Pass A (Complete)", a_summary, Colors.GREEN)

            # Pass B progress bar (use current_idx for in-progress display)
            b_display = b_idx if b_idx > 0 else b_done
            bar = _progress_bar(b_display, b_total) if b_total > 0 else f"{b_done} done"
            row("Pass B Progress", bar, Colors.CYAN)

            # Current scenario trial progress
            if info['current_trial_total'] > 0:
                trial_str = f"{info['current_trial']}/{info['current_trial_total']} trials"
                row("Current Scenario", trial_str, Colors.BLUE)

        # Total candidates accumulated
        if info['total_candidates'] > 0:
            row("Candidates Found", f"{info['total_candidates']:,}", Colors.GREEN)

        # Best score
        if info['best_score'] is not None:
            s = info['best_score']
            if s < 0:
                label, sstr = "Best Binding Score", f"{s:.2f} kcal/mol"
                color = Colors.GREEN if s < -7 else Colors.YELLOW
            else:
                label, sstr = "Best Oracle Score", f"{s:.4f}"
                color = Colors.GREEN if s > 0.5 else Colors.YELLOW
            row(label, sstr, color)

        # Elapsed
        if info['elapsed_time']:
            row("Elapsed Time", info['elapsed_time'], Colors.BLUE)

        # Last update
        if info['last_update']:
            row("Last Update", info['last_update'], Colors.BLUE)

        # ETA — interpolated real-time within current in-progress scenario
        now = datetime.now()
        if phase == 'pass_a' and a_total > 0 and a_done > 0:
            elapsed_in_cur = 0.0
            if info['last_a_scenario_start_dt']:
                elapsed_in_cur = (now - info['last_a_scenario_start_dt']).total_seconds()
            remaining = a_total - a_done
            eta = _eta_str(remaining, info['pass_a_times'],
                           default_secs=90, elapsed_in_current=elapsed_in_cur)
            row("ETA (Pass A total)", eta, Colors.GREEN)

            # Show time in current scenario
            if info['last_a_scenario_start_dt']:
                cur_elapsed = int((now - info['last_a_scenario_start_dt']).total_seconds())
                avg = _avg_secs(info['pass_a_times'], 90)
                cur_str = f"{_fmt_duration(cur_elapsed)} / ~{_fmt_duration(int(avg))} avg"
                row("Current Scenario", cur_str, Colors.BLUE)

        elif phase == 'pass_b' and b_total > 0:
            elapsed_in_cur = 0.0
            if info['last_b_scenario_start_dt']:
                elapsed_in_cur = (now - info['last_b_scenario_start_dt']).total_seconds()
            remaining = b_total - b_done
            eta = _eta_str(remaining, info['pass_b_times'],
                           default_secs=300, elapsed_in_current=elapsed_in_cur)
            row("ETA (Pass B total)", eta, Colors.GREEN)

            # Show time in current scenario vs avg
            if info['last_b_scenario_start_dt']:
                cur_elapsed = int((now - info['last_b_scenario_start_dt']).total_seconds())
                avg = _avg_secs(info['pass_b_times'], 300)
                cur_str = f"{_fmt_duration(cur_elapsed)} / ~{_fmt_duration(int(avg))} avg"
                row("Current Scenario", cur_str, Colors.BLUE)

    elif phase in ('validation', 'physics'):
        # Show summarized Pass A + B counts first
        if info['pass_a_total'] > 0:
            row("Pass A", f"{info['pass_a_completed']}/{info['pass_a_total']} done", Colors.GREEN)
        if info['pass_b_total'] > 0:
            row("Pass B", f"{info['pass_b_completed']}/{info['pass_b_total']} done", Colors.GREEN)

        if info['total_candidates'] > 0:
            row("Candidates Found", f"{info['total_candidates']:,}", Colors.GREEN)

        if info['best_score'] is not None:
            s = info['best_score']
            if s < 0:
                row("Best Binding Score", f"{s:.2f} kcal/mol", Colors.GREEN)
            else:
                row("Best Oracle Score", f"{s:.4f}", Colors.GREEN if s > 0.5 else Colors.YELLOW)

        if info['elapsed_time']:
            row("Elapsed Time", info['elapsed_time'], Colors.BLUE)
        if info['last_update']:
            row("Last Update", info['last_update'], Colors.BLUE)

        if phase == 'physics' and info['physics_validation_count'] > 0:
            pct = (info['physics_validation_count'] / 30) * 100
            bar = _progress_bar(info['physics_validation_count'], 30)
            row("Physics Validation", bar, Colors.CYAN)
        elif phase == 'validation' and info['validation_rank'] > 0:
            bar = _progress_bar(info['validation_rank'], 30)
            row("ESMFold Validation", bar, Colors.CYAN)
            if info['validation_target']:
                row("Validating Target", info['validation_target'][:48], Colors.CYAN)

        # ETA for validation
        ts_list = info['physics_timestamps'] if phase == 'physics' else info['validation_timestamps']
        current = info['physics_validation_count'] if phase == 'physics' else info['validation_rank']
        if len(ts_list) >= 2:
            sorted_ts = sorted(ts_list, key=lambda x: x[0])
            diffs = []
            for i in range(1, len(sorted_ts)):
                try:
                    t1 = datetime.strptime(sorted_ts[i-1][1], '%Y-%m-%d %H:%M:%S')
                    t2 = datetime.strptime(sorted_ts[i][1], '%Y-%m-%d %H:%M:%S')
                    diffs.append((t2 - t1).total_seconds())
                except:
                    pass
            if diffs:
                remaining = 30 - current
                eta = _eta_str(remaining, diffs, default_secs=15 if phase == 'physics' else 120)
                row("Estimated ETA", eta, Colors.GREEN)

    elif phase == 'complete':
        row("Status", "Campaign finished!", Colors.GREEN)
        if info['total_candidates'] > 0:
            row("Total Candidates", f"{info['total_candidates']:,}", Colors.GREEN)
        if info['best_score'] is not None:
            s = info['best_score']
            row("Best Score", f"{s:.4f}", Colors.GREEN)
        if info['elapsed_time']:
            row("Total Time", info['elapsed_time'], Colors.BLUE)

    # ── System resources ─────────────────────────────────────────────────────
    resources = get_system_resources()
    if resources['available'] or MONITORING_AVAILABLE:
        print(f"{Colors.BOLD}├{'─'*78}┤{Colors.ENDC}")
        print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}System Resources{Colors.ENDC}                                                         {Colors.BOLD}│{Colors.ENDC}")
        print(f"{Colors.BOLD}├{'─'*78}┤{Colors.ENDC}")

        cpu_color = (Colors.GREEN if resources['cpu_percent'] < 70
                     else Colors.YELLOW if resources['cpu_percent'] < 90
                     else Colors.RED)
        row("CPU Usage", f"{resources['cpu_percent']:.1f}%", cpu_color)

        if resources['ram_total_gb'] > 0:
            ram_pct = resources['ram_used_gb'] / resources['ram_total_gb'] * 100
            ram_color = Colors.GREEN if ram_pct < 70 else Colors.YELLOW if ram_pct < 90 else Colors.RED
            row("System RAM", f"{resources['ram_used_gb']:.1f} / {resources['ram_total_gb']:.1f} GB ({ram_pct:.0f}%)", ram_color)

        if resources['gpu_util'] > 0 or resources['gpu_mem_total_gb'] > 0:
            gpu_color = (Colors.GREEN if resources['gpu_util'] > 40
                         else Colors.YELLOW if resources['gpu_util'] > 10
                         else Colors.CYAN)
            row("GPU Usage", f"{resources['gpu_util']:.1f}%", gpu_color)
            if resources['gpu_mem_total_gb'] > 0:
                vram_pct = resources['gpu_mem_used_gb'] / resources['gpu_mem_total_gb'] * 100
                vram_color = Colors.GREEN if vram_pct < 80 else Colors.YELLOW if vram_pct < 95 else Colors.RED
                row("GPU VRAM", f"{resources['gpu_mem_used_gb']:.1f} / {resources['gpu_mem_total_gb']:.1f} GB ({vram_pct:.0f}%)", vram_color)

    print(f"{Colors.BOLD}└{'─'*78}┘{Colors.ENDC}\n")

    # Warnings / errors (last 3)
    if info['warnings']:
        print(f"{Colors.YELLOW}⚠️  Warnings: {len(info['warnings'])}{Colors.ENDC}")
        for w in info['warnings'][-3:]:
            print(f"   {Colors.YELLOW}•{Colors.ENDC} {w[:75]}")
        print()

    if info['errors']:
        print(f"{Colors.RED}❌ Errors: {len(info['errors'])}{Colors.ENDC}")
        for e in info['errors'][-3:]:
            print(f"   {Colors.RED}•{Colors.ENDC} {e[:75]}")
        print()


def monitor_campaign(log_path, refresh_interval=2):
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
