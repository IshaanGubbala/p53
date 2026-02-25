#!/usr/bin/env python3
"""
Real-time campaign progress monitor with comprehensive detail.

Captures:
- All optimization features (TF32, SDPA, shallow gradients, etc.)
- Loss term breakdown during trials
- Structural awareness metrics
- ESMFold gate status
- Population sharing stats
- GPU/CPU utilization with model memory breakdown

Usage:
    python scripts/monitor_campaign.py

Or run in a separate terminal while the campaign is running.
"""

import os
import sys
import time
import json
from pathlib import Path
from datetime import datetime, timedelta
import re
from collections import defaultdict

# System monitoring
try:
    import psutil
    import pynvml
    import torch
    MONITORING_AVAILABLE = True
except ImportError:
    MONITORING_AVAILABLE = False
    psutil = None
    pynvml = None
    torch = None

# ANSI color codes
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    MAGENTA = '\033[95m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def clear_screen():
    """Clear terminal screen."""
    os.system('cls' if os.name == 'nt' else 'clear')

def get_system_resources():
    """Get detailed system resource usage."""
    resources = {
        'cpu_percent': 0,
        'ram_used_gb': 0,
        'ram_total_gb': 0,
        'gpu_util': 0,
        'gpu_mem_used_gb': 0,
        'gpu_mem_total_gb': 0,
        'gpu_memory_reserved_gb': 0,
        'gpu_compute_capability': '',
        'available': False,
        'cuda_available': False,
    }

    if not MONITORING_AVAILABLE:
        return resources

    try:
        resources['cpu_percent'] = psutil.cpu_percent(interval=0.1)
        mem = psutil.virtual_memory()
        resources['ram_used_gb'] = mem.used / (1024**3)
        resources['ram_total_gb'] = mem.total / (1024**3)

        try:
            import torch
            resources['cuda_available'] = torch.cuda.is_available()
            if resources['cuda_available']:
                resources['available'] = True
                pynvml.nvmlInit()
                handle = pynvml.nvmlDeviceGetHandleByIndex(0)
                
                # GPU utilization
                util = pynvml.nvmlDeviceGetUtilizationRates(handle)
                resources['gpu_util'] = util.gpu
                
                # Memory info
                mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
                resources['gpu_mem_used_gb'] = mem_info.used / (1024**3)
                resources['gpu_mem_total_gb'] = mem_info.total / (1024**3)
                
                # Compute capability - use cuda device properties
                try:
                    import torch as _torch
                    if _torch.cuda.is_available():
                        cc = _torch.cuda.get_device_capability(0)
                        resources['gpu_compute_capability'] = f"{cc[0]}.{cc[1]}"
                        resources['gpu_memory_reserved_gb'] = _torch.cuda.memory_reserved(0) / (1024**3)
                except Exception:
                    resources['gpu_compute_capability'] = 'N/A'
                    resources['gpu_memory_reserved_gb'] = 0
                
                pynvml.nvmlShutdown()
        except Exception as e:
            pass

    except:
        pass

    return resources


def get_runtime_capabilities():
    """Get runtime capabilities from log."""
    caps = {
        'python': '',
        'torch': '',
        'cuda': False,
        'mps': False,
        'transformers': '',
        'tf32': False,
        'cudnn_benchmark': False,
        'sdpa': False,
        'flash_attention': False,
    }
    
    log_path = Path("logs/p53cad_workflow.log")
    if not log_path.exists():
        return caps
    
    try:
        with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        for line in lines[-50:]:  # Check recent lines
            if 'Runtime capabilities' in line:
                # Parse: python=3.11.14 torch=2.3.1+cu121 mps=False cuda=True
                match = re.search(r'python=([\d.]+)', line)
                if match:
                    caps['python'] = match.group(1)
                match = re.search(r'torch=([\d.]+)', line)
                if match:
                    caps['torch'] = match.group(1)
                match = re.search(r'cuda=(\w+)', line)
                if match:
                    caps['cuda'] = match.group(1) == 'True'
                match = re.search(r'mps=(\w+)', line)
                if match:
                    caps['mps'] = match.group(1) == 'True'
                match = re.search(r'transformers=([\d.]+)', line)
                if match:
                    caps['transformers'] = match.group(1)
            
            if 'Optimization config:' in line:
                # Parse: {'shallow_gradients': True, 'batched_trials': False, ...}
                match = re.search(r"'tf32_enabled': (\w+)", line)
                if match:
                    caps['tf32'] = match.group(1) == 'True'
                match = re.search(r"'cudnn_benchmark': (\w+)", line)
                if match:
                    caps['cudnn_benchmark'] = match.group(1) == 'True'
                match = re.search(r"'sdpa_fast_path': (\w+)", line)
                if match:
                    caps['sdpa'] = match.group(1) == 'True'
    except:
        pass
    
    return caps


def parse_log_file(log_path):
    """Parse the log file for comprehensive progress information."""
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
        'current_phase': 'pass_a',
        'current_step': 'Starting',

        # Pass A tracking
        'pass_a_total': 0,
        'pass_a_completed': 0,

        # Pass B tracking
        'pass_b_total': 0,
        'pass_b_completed': 0,
        'pass_b_current_idx': 0,
        'pass_b_from_scenario': None,

        # Pipeline stages tracking
        'pipeline_stage': 'initializing',  # initializing, pass_a, pass_b, ranking, selecting, physics, complete
        'stage_details': {},

        # ESMFold gate
        'esmfold_gate_enabled': False,
        'esmfold_gate_passed': 0,
        'esmfold_gate_filtered': 0,

        # Population sharing
        'population_share_size': 0,

        # Scenario start timestamps
        'last_a_scenario_start_dt': None,
        'last_b_scenario_start_dt': None,

        # Current scenario info
        'current_scenario_id': '',
        'current_trial': 0,
        'current_trial_total': 0,

        # Trial stages within scenario
        'trial_stage': 'initializing',  # loading, scoring, optimizing, evaluating
        'trial_stage_detail': '',
        
        # Deep optimization details
        'deep_optimization_iter': 0,
        'deep_optimization_total': 0,
        'deep_best_score': None,

        # Current loss breakdown
        'loss_total': 0.0,
        'loss_terms': {},

        # Candidates / scores
        'total_candidates': 0,
        'pass_a_candidates': 0,
        'pass_b_candidates': 0,
        'best_score': None,
        'best_identity': 0.0,
        'best_stability': 0.0,
        'best_binding': 0.0,
        'best_mutations': [],
        
        # Top candidates from Pass A (for Pass B selection)
        'top_pass_a_candidates': [],
        
        # Ranking info
        'ranking_in_progress': False,
        'ranking_total': 0,
        'ranking_completed': 0,

        # Selection info
        'selection_in_progress': False,
        'selected_for_pass_b': 0,
        'selection_criteria': '',

        # Identity floor
        'identity_floor': 94.0,

        # Validation
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
    pass_a_times = []
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
            info['pipeline_stage'] = 'pass_a'
            info['status'] = 'Running'
            campaign_started = True
            if not start_time and ts_dt:
                start_time = ts_dt
            a_idx = int(pass_a_match.group(1))
            a_total = int(pass_a_match.group(2))
            info['pass_a_total'] = a_total
            info['pass_a_completed'] = max(info['pass_a_completed'], a_idx - 1)
            if ts_dt:
                info['last_a_scenario_start_dt'] = ts_dt
            
            # Extract scenario ID
            scenario_match = re.search(r'Pass A \[\d+/\d+\]\s+(\S+)', line)
            if scenario_match:
                info['current_scenario_id'] = scenario_match.group(1)

        # Pass B line: "Pass B [idx/total] scenario_id"
        pass_b_match = re.search(r'Pass B \[(\d+)/(\d+)\]', line)
        if pass_b_match:
            info['current_phase'] = 'pass_b'
            info['current_step'] = 'Pass B: Deep optimizing top scenarios'
            info['pipeline_stage'] = 'pass_b'
            info['status'] = 'Running'
            campaign_started = True
            b_idx = int(pass_b_match.group(1))
            b_total = int(pass_b_match.group(2))
            info['pass_b_total'] = b_total
            info['pass_b_current_idx'] = b_idx
            info['pass_b_completed'] = max(info['pass_b_completed'], b_idx - 1)
            if ts_dt:
                info['last_b_scenario_start_dt'] = ts_dt
            
            scenario_match = re.search(r'Pass B \[\d+/\d+\]\s+(\S+)', line)
            if scenario_match:
                info['current_scenario_id'] = scenario_match.group(1)
            
            # Check if Pass B indicates which Pass A scenario this came from
            from_match = re.search(r'from.*scenario.*(\d+)', line, re.IGNORECASE)
            if from_match:
                info['pass_b_from_scenario'] = from_match.group(1)

        # ESMFold gate
        if 'ESMFold gate:' in line:
            info['esmfold_gate_enabled'] = True
            passed_match = re.search(r'(\d+)/(\d+)\s+passed', line)
            if passed_match:
                info['esmfold_gate_passed'] = int(passed_match.group(1))
                info['esmfold_gate_filtered'] = int(passed_match.group(2)) - int(passed_match.group(1))

        # Population sharing
        if 'Population share pool size:' in line:
            match = re.search(r'(\d+)', line)
            if match:
                info['population_share_size'] = int(match.group(1))

        # Identity floor detection
        if '_delivery_identity_floor' in line or 'identity_floor' in line.lower():
            match = re.search(r'(\d+\.?\d*)%?', line)
            if match:
                info['identity_floor'] = float(match.group(1))

        # Trial progress with detailed loss terms
        if 'Progress:' in line and 'trials' in line:
            trial_match = re.search(r'(\d+)/(\d+)\s*trials', line)
            if trial_match:
                info['current_trial'] = int(trial_match.group(1))
                info['current_trial_total'] = int(trial_match.group(2))
            
            # Extract best score
            score_match = re.search(r'best=([-+]?\d*\.?\d+)', line)
            if score_match:
                try:
                    score = float(score_match.group(1))
                    if info['best_score'] is None or score > info['best_score']:
                        info['best_score'] = score
                except:
                    pass
        
        # Trial stage detection (within a scenario)
        if 'Loading variant' in line or 'loading' in line.lower():
            info['trial_stage'] = 'loading'
            info['trial_stage_detail'] = 'Loading protein variant'
        elif 'Scoring variant' in line or 'scoring' in line.lower():
            info['trial_stage'] = 'scoring'
            info['trial_stage_detail'] = 'Scoring with oracle'
        elif 'Optimizing' in line or 'optimizing' in line.lower() or 'gradient' in line.lower():
            info['trial_stage'] = 'optimizing'
            info['trial_stage_detail'] = 'Gradient optimization'
        elif 'Evaluating' in line or 'evaluating' in line.lower():
            info['trial_stage'] = 'evaluating'
            info['trial_stage_detail'] = 'Evaluating candidate'
        
        # Deep optimization iterations
        deep_iter_match = re.search(r'Deep optimization.*?(\d+)/(\d+)', line)
        if deep_iter_match:
            info['deep_optimization_iter'] = int(deep_iter_match.group(1))
            info['deep_optimization_total'] = int(deep_iter_match.group(2))
            info['pipeline_stage'] = 'deep_opt'
        
        # Ranking stage
        if 'Ranking' in line or 'ranking' in line.lower():
            info['ranking_in_progress'] = True
            info['pipeline_stage'] = 'ranking'
            rank_match = re.search(r'(\d+)/(\d+)', line)
            if rank_match:
                info['ranking_completed'] = int(rank_match.group(1))
                info['ranking_total'] = int(rank_match.group(2))
        
        # Selection stage
        if 'Selecting' in line or 'selecting top' in line.lower() or 'Selecting best' in line:
            info['selection_in_progress'] = True
            info['pipeline_stage'] = 'selecting'
            sel_match = re.search(r'(\d+)\s*selected', line)
            if sel_match:
                info['selected_for_pass_b'] = int(sel_match.group(1))
        
        # ESMFold structure prediction
        if 'ESMFold' in line or 'ESM' in line:
            if 'predicting' in line.lower() or 'running' in line.lower():
                info['stage_details']['esmfold'] = 'Predicting structure'
        
        # Pass A candidates collected
        if 'Pass A' in line and 'candidates' in line:
            cand_match = re.search(r'(\d+)\s*candidates', line)
            if cand_match:
                info['pass_a_candidates'] = int(cand_match.group(1))
        
        # Pass B candidates collected  
        if 'Pass B' in line and 'candidates' in line:
            cand_match = re.search(r'(\d+)\s*candidates', line)
            if cand_match:
                info['pass_b_candidates'] = int(cand_match.group(1))

        # Trajectory loss terms - parse all the new loss components
        if '"step":' in line and '"loss_total":' in line:
            try:
                # Try to extract JSON from line
                json_match = re.search(r'\{[^{}]*\}', line)
                if json_match:
                    data = json.loads(json_match.group(0))
                    
                    # Store current loss breakdown
                    info['loss_total'] = data.get('loss_total', 0)
                    
                    loss_keys = [
                        'loss_score_term', 'loss_stability_term', 'loss_binding_term',
                        'loss_hydrophobic_term', 'loss_ood_penalty', 'loss_mutation_penalty',
                        'loss_identity_penalty', 'loss_stability_penalty', 'loss_binding_penalty',
                        'loss_lock_penalty', 'loss_l1_penalty', 'loss_pll_term',
                        'loss_contact_penalty', 'loss_epistasis_penalty', 'loss_dms_penalty',
                        'loss_cancer_pll_term', 'loss_local_score_term', 'loss_cond_rescue_term',
                        'loss_structure_term', 'loss_mut_conf_penalty', 'loss_interface_floor_penalty',
                        'loss_crater_penalty', 'loss_contact_reg_penalty'
                    ]
                    
                    for key in loss_keys:
                        if key in data:
                            info['loss_terms'][key] = data[key]
                    
                    # Extract candidate metrics
                    if 'identity' in data:
                        info['best_identity'] = max(info['best_identity'], data['identity'])
                    if 'stability' in data:
                        info['best_stability'] = max(info['best_stability'], data['stability'])
                    if 'binding' in data:
                        info['best_binding'] = max(info['best_binding'], data['binding'])
                    
                    # Track best mutations
                    if 'n_mutations' in data:
                        if data.get('score', -999) > 0.5:
                            info['best_mutations'].append({
                                'n': data['n_mutations'],
                                'identity': data.get('identity', 0),
                            })
            except:
                pass

        # Scenario complete
        if 'Scenario' in line and 'complete:' in line:
            if ts_dt:
                if last_scenario_ts:
                    diff = (ts_dt - last_scenario_ts).total_seconds()
                    if info['current_phase'] == 'pass_a':
                        pass_a_times.append(diff)
                    else:
                        pass_b_times.append(diff)
                last_scenario_ts = ts_dt

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
            info['current_step'] = 'Physics Validation'
            info['pipeline_stage'] = 'physics'
            info['status'] = 'Running'
            info['in_physics_validation'] = True
            info['in_validation'] = True
            campaign_started = True
            info['physics_validation_count'] += 1
        
        # Physics validation target
        if 'Validating candidate' in line or 'Validating:' in line:
            match = re.search(r'(?:candidate|from)[:\s]+(\S+)', line)
            if match:
                info['validation_target'] = match.group(1)
        
        # Validation rank
        if 'rank' in line.lower() and ('#' in line or 'rank' in line):
            rank_match = re.search(r'rank.*?(\d+)', line, re.IGNORECASE)
            if rank_match:
                info['validation_rank'] = int(rank_match.group(1))
        
        # Energy terms from physics
        if 'delta_delta_g' in line.lower() or 'ddg' in line.lower():
            ddg_match = re.search(r'(?:delta_delta_g|ddg)[:\s=]+([-]?\d+\.?\d*)', line, re.IGNORECASE)
            if ddg_match:
                info['stage_details']['ddg'] = float(ddg_match.group(1))

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

    info['pass_a_times'] = pass_a_times
    info['pass_b_times'] = pass_b_times

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


def _eta_str(remaining: int, times: list, default_secs: int, elapsed_in_current: float = 0.0) -> str:
    """Format ETA string."""
    avg = _avg_secs(times, default_secs)
    eta = max(0, int(remaining * avg - elapsed_in_current))
    return f"{_fmt_duration(eta)} ({remaining} scenarios)"


def _progress_bar(done: int, total: int, width: int = 20) -> str:
    """Simple ASCII progress bar."""
    if total <= 0:
        return f"{'?':>{width}}"
    frac = min(done / total, 1.0)
    filled = int(frac * width)
    bar = '█' * filled + '░' * (width - filled)
    pct = frac * 100
    return f"{bar} {pct:.0f}%"


def print_header():
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*90}{Colors.ENDC}")
    print(f"{Colors.BOLD}{Colors.CYAN}  p53CAD Campaign Monitor - Structural Awareness Pipeline{Colors.ENDC}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*90}{Colors.ENDC}\n")


def print_runtime_info():
    """Print runtime capabilities and optimization settings."""
    caps = get_runtime_capabilities()
    
    print(f"{Colors.BOLD}┌{'─'*88}┐{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Runtime Configuration{Colors.ENDC}" + " " * 65 + f"{Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
    
    def config_row(label, value, color=Colors.GREEN):
        print(f"{Colors.BOLD}│{Colors.ENDC}   {label:<30} │ {color}{str(value):<54}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    
    config_row("Python", caps['python'] or 'N/A')
    config_row("PyTorch", caps['torch'] or 'N/A')
    config_row("Transformers", caps['transformers'] or 'N/A')
    config_row("CUDA Available", "Yes" if caps['cuda'] else "No", Colors.GREEN if caps['cuda'] else Colors.RED)
    config_row("TF32 Enabled", "Yes" if caps['tf32'] else "No", Colors.GREEN if caps['tf32'] else Colors.YELLOW)
    config_row("cuDNN Benchmark", "Yes" if caps['cudnn_benchmark'] else "No", Colors.GREEN if caps['cudnn_benchmark'] else Colors.YELLOW)
    config_row("SDPA Fast Path", "Yes" if caps['sdpa'] else "No", Colors.GREEN if caps['sdpa'] else Colors.YELLOW)
    
    print(f"{Colors.BOLD}└{'─'*88}┘{Colors.ENDC}\n")


def print_status_table(info):
    if not info:
        print(f"{Colors.YELLOW}⏳ Waiting for campaign to start...{Colors.ENDC}\n")
        return

    phase = info['current_phase']
    status_color = (Colors.GREEN if info['status'] == 'Running'
                    else Colors.RED if info['status'] == 'Error'
                    else Colors.YELLOW)

    # ── Header ─────────────────────────────────────────────────────────────────
    print(f"{Colors.BOLD}┌{'─'*88}┐{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Status:{Colors.ENDC} {status_color}{info['status']:<77}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Phase:{Colors.ENDC}  {Colors.CYAN}{info['current_step']:<77}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")

    def row(label, value, color=Colors.GREEN):
        print(f"{Colors.BOLD}│{Colors.ENDC} {label:<28} │ {color}{str(value):<56}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")

    # ── Pipeline Overview ───────────────────────────────────────────────────────
    stage = info.get('pipeline_stage', phase)
    if stage == 'initializing':
        stage = phase
    
    stage_names = {
        'initializing': 'Initializing Campaign',
        'pass_a': 'Pass A - Screening Scenarios',
        'pass_b': 'Pass B - Deep Optimization',
        'ranking': 'Ranking Candidates',
        'selecting': 'Selecting Top Candidates',
        'deep_opt': 'Deep Optimization',
        'physics': 'Physics Validation',
        'complete': 'Campaign Complete'
    }
    current_stage_name = stage_names.get(stage, stage)
    
    # Determine overall pipeline progress
    a_done = info['pass_a_completed']
    a_total = info['pass_a_total']
    b_done = info['pass_b_completed']
    b_total = info['pass_b_total']
    
    # Pipeline stage indicators based on current_phase
    stage_indicators = []
    if phase in ('pass_a', 'pass_b', 'physics', 'complete'):
        stage_indicators.append(('INIT', '✓'))
    else:
        stage_indicators.append(('INIT', '○'))
    
    if a_done == a_total and a_total > 0:
        stage_indicators.append(('PASS_A', '✓'))
    elif phase == 'pass_a':
        stage_indicators.append(('PASS_A', '▓'))
    else:
        stage_indicators.append(('PASS_A', '○'))
    
    # Ranking happens after Pass A completes
    if b_done > 0 or phase in ('pass_b', 'physics', 'complete'):
        stage_indicators.append(('RANK', '✓'))
    elif phase == 'ranking':
        stage_indicators.append(('RANK', '▓'))
    else:
        stage_indicators.append(('RANK', '○'))
    
    if b_done == b_total and b_total > 0:
        stage_indicators.append(('PASS_B', '✓'))
    elif phase == 'pass_b':
        stage_indicators.append(('PASS_B', '▓'))
    else:
        stage_indicators.append(('PASS_B', '○'))
    
    if phase in ('physics', 'complete'):
        stage_indicators.append(('PHYS', '✓'))
    elif phase == 'physics':
        stage_indicators.append(('PHYS', '▓'))
    else:
        stage_indicators.append(('PHYS', '○'))
    
    if phase == 'complete':
        stage_indicators.append(('DONE', '✓'))
    else:
        stage_indicators.append(('DONE', '○'))
    
    stage_str = ' '.join([f"{Colors.CYAN}{name}{Colors.ENDC}:{Colors.GREEN if ind == '✓' else Colors.YELLOW if ind == '▓' else Colors.BLUE}{ind}{Colors.ENDC}" for name, ind in stage_indicators])
    
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Pipeline Stage:{Colors.ENDC} {Colors.YELLOW}{current_stage_name:<60}{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    row("  Pipeline", stage_str, Colors.BLUE)
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")

    # ── Detailed Pass A Progress ───────────────────────────────────────────────
    a_pct = (a_done / a_total * 100) if a_total > 0 else 0
    
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}{Colors.MAGENTA}═══════════════════════════════════════════════════════════════════════{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}  PASS A: SCREENING ALL 108 SCENARIOS{Colors.ENDC}" + " " * 51 + f"{Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}{Colors.MAGENTA}═══════════════════════════════════════════════════════════════════════{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    
    # Detailed scenario progress
    if phase == 'pass_a':
        bar = _progress_bar(a_done, a_total, 40)
        row("  Overall Progress", bar, Colors.CYAN)
        row("  Scenario", f"{a_done + 1:>3} / {a_total:>3} ({a_pct:>5.1f}%)", Colors.GREEN if a_done > 50 else Colors.YELLOW)
        if info['current_scenario_id']:
            row("  Current Scenario ID", f"{info['current_scenario_id'][:55]}", Colors.BLUE)
        # Show ETA
        now = datetime.now()
        if a_total > 0 and a_done > 0:
            elapsed_in_cur = 0.0
            if info['last_a_scenario_start_dt']:
                elapsed_in_cur = (now - info['last_a_scenario_start_dt']).total_seconds()
            remaining = a_total - a_done
            avg_time = _avg_secs(info['pass_a_times'], 90)
            eta_secs = int(remaining * avg_time - elapsed_in_cur)
            row("  ETA", f"{_fmt_duration(max(0, eta_secs))} (scen {remaining} remaining)", Colors.GREEN)
    else:
        row("  Status", f"COMPLETED - All {a_done}/{a_total} scenarios screened", Colors.GREEN)
    
    # Pass A timing stats
    if info['pass_a_times']:
        avg_a = _avg_secs(info['pass_a_times'], 90)
        min_a = min(info['pass_a_times']) if info['pass_a_times'] else 0
        max_a = max(info['pass_a_times']) if info['pass_a_times'] else 0
        row("  Avg Scenario Time", f"{avg_a:.1f}s", Colors.BLUE)

    # ── Detailed Pass B Progress ────────────────────────────────────────────────
    b_done = info['pass_b_completed']
    b_total = info['pass_b_total']
    b_idx = info['pass_b_current_idx']
    b_display = b_idx if b_idx > 0 else b_done
    b_pct = (b_done / b_total * 100) if b_total > 0 else 0
    
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}{Colors.MAGENTA}═══════════════════════════════════════════════════════════════════════{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}  PASS B: DEEP OPTIMIZATION OF TOP SCENARIOS{Colors.ENDC}" + " " * 47 + f"{Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}{Colors.MAGENTA}═══════════════════════════════════════════════════════════════════════{Colors.ENDC} {Colors.BOLD}│{Colors.ENDC}")
    
    if phase == 'pass_b' or (phase in ('pass_a', 'complete') and b_total > 0):
        if phase == 'pass_b':
            bar = _progress_bar(b_display, b_total, 40)
            row("  Overall Progress", bar, Colors.CYAN)
            row("  Scenario", f"{b_done + 1:>3} / {b_total:>3} ({b_pct:>5.1f}%)", Colors.GREEN if b_done > 5 else Colors.YELLOW)
            if info['current_scenario_id']:
                row("  Current Scenario ID", f"{info['current_scenario_id'][:55]}", Colors.BLUE)
            # Show which pass A scenario this came from
            row("  Parent (Pass A)", f"Scenario {info.get('pass_b_from_scenario', 'N/A')}", Colors.MAGENTA)
            # ETA
            now = datetime.now()
            if b_total > 0 and b_done > 0:
                elapsed_in_cur = 0.0
                if info['last_b_scenario_start_dt']:
                    elapsed_in_cur = (now - info['last_b_scenario_start_dt']).total_seconds()
                remaining = b_total - b_done
                avg_time = _avg_secs(info['pass_b_times'], 300)
                eta_secs = int(remaining * avg_time - elapsed_in_cur)
                row("  ETA", f"{_fmt_duration(max(0, eta_secs))} (scen {remaining} remaining)", Colors.GREEN)
        elif phase == 'complete':
            row("  Status", f"COMPLETED - All {b_done}/{b_total} scenarios deep optimized", Colors.GREEN)
        else:
            row("  Status", f"Pending ({b_total} scenarios queued)", Colors.YELLOW)

    # Pass B timing stats
    if info['pass_b_times']:
        avg_b = _avg_secs(info['pass_b_times'], 300)
        row("  Avg Scenario Time", f"{avg_b:.1f}s", Colors.BLUE)

    # ── Trial Progress (within current scenario) ────────────────────────────────
    if info['current_trial_total'] > 0:
        print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
        print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}  Trial Progress (Current Scenario){Colors.ENDC}" + " " * 56 + f"{Colors.BOLD}│{Colors.ENDC}")
        print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
        trial_bar = _progress_bar(info['current_trial'], info['current_trial_total'], 40)
        row("  Trial Iteration", trial_bar, Colors.MAGENTA)
        row("  Progress", f"{info['current_trial']} / {info['current_trial_total']} ({info['current_trial']/info['current_trial_total']*100:.1f}%)", Colors.CYAN)
        
        # Trial stage
        trial_stage = info.get('trial_stage', 'unknown')
        trial_stage_detail = info.get('trial_stage_detail', '')
        stage_color = Colors.GREEN if trial_stage in ('evaluating', 'scoring') else Colors.YELLOW if trial_stage == 'optimizing' else Colors.BLUE
        row("  Current Stage", f"{trial_stage.upper()} - {trial_stage_detail}", stage_color)
        
        # Deep optimization details
        if info['deep_optimization_total'] > 0:
            deep_bar = _progress_bar(info['deep_optimization_iter'], info['deep_optimization_total'], 30)
            row("  Deep Optimization", deep_bar, Colors.CYAN)

    # ── Ranking & Selection ───────────────────────────────────────────────────
    if info['ranking_in_progress'] or info['ranking_total'] > 0:
        print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
        print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}  Ranking & Selection{Colors.ENDC}" + " " * 67 + f"{Colors.BOLD}│{Colors.ENDC}")
        print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
        if info['ranking_in_progress']:
            rank_bar = _progress_bar(info['ranking_completed'], info['ranking_total'], 40)
            row("  Ranking Progress", rank_bar, Colors.CYAN)
        if info['selected_for_pass_b'] > 0:
            row("  Selected for Pass B", f"{info['selected_for_pass_b']} candidates", Colors.GREEN)
        if info['pass_a_candidates'] > 0:
            row("  Pass A Candidates", f"{info['pass_a_candidates']:,} total", Colors.BLUE)

    # ── Physics Validation ─────────────────────────────────────────────────────
    if info['in_physics_validation'] or info['physics_validation_count'] > 0:
        print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
        print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}  Physics Validation (Energy Minimization){Colors.ENDC}" + " " * 51 + f"{Colors.BOLD}│{Colors.ENDC}")
        print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
        row("  Status", "Running" if info['in_physics_validation'] else "Completed", Colors.GREEN if not info['in_physics_validation'] else Colors.YELLOW)
        row("  Validations Run", f"{info['physics_validation_count']}", Colors.CYAN)
        if info['validation_target']:
            row("  Current Target", f"{info['validation_target'][:50]}", Colors.BLUE)
        if info['validation_rank'] > 0:
            row("  Validation Rank", f"#{info['validation_rank']}", Colors.MAGENTA)
        
        # Show delta G if available
        if 'ddg' in info.get('stage_details', {}):
            ddg = info['stage_details']['ddg']
            ddg_color = Colors.GREEN if ddg < 0 else Colors.RED
            row("  Delta Delta G", f"{ddg:.2f} kcal/mol", ddg_color)

    # ── Structural Awareness Features ───────────────────────────────────────
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Structural Awareness Features{Colors.ENDC}" + " " * 63 + f"{Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
    
    row("Identity Floor", f"{info['identity_floor']:.0f}% (surgical rescue)", Colors.YELLOW)
    
    if info['esmfold_gate_enabled']:
        gate_str = f"{info['esmfold_gate_passed']} passed, {info['esmfold_gate_filtered']} filtered"
        row("ESMFold Gate", gate_str, Colors.GREEN if info['esmfold_gate_filtered'] == 0 else Colors.YELLOW)
    
    if info['population_share_size'] > 0:
        row("Population Share", f"{info['population_share_size']} candidates", Colors.CYAN)

    # ── Loss Term Breakdown ─────────────────────────────────────────────────
    if info['loss_terms']:
        print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
        print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Current Loss Breakdown{Colors.ENDC}" + " " * 69 + f"{Colors.BOLD}│{Colors.ENDC}")
        print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
        
        loss_keys = [
            ('loss_score_term', 'Oracle Score'),
            ('loss_stability_term', 'Stability'),
            ('loss_binding_term', 'Binding'),
            ('loss_identity_penalty', 'Identity'),
            ('loss_mutation_penalty', 'Mutation Count'),
            ('loss_contact_penalty', 'Contact'),
            ('loss_epistasis_penalty', 'Epistasis'),
            ('loss_ood_penalty', 'OOD'),
            ('loss_struct_term', 'Structure'),
            ('loss_mut_conf_penalty', 'Mut Confidence'),
            ('loss_interface_floor_penalty', 'Interface Floor'),
            ('loss_crater_penalty', 'Crater'),
            ('loss_contact_reg_penalty', 'Contact Reg'),
            ('loss_cond_rescue_term', 'Conditional Rescue'),
        ]
        
        for key, label in loss_keys:
            if key in info['loss_terms']:
                val = info['loss_terms'][key]
                color = Colors.RED if abs(val) > 1.0 else Colors.GREEN
                row(f"  {label}", f"{val:.4f}", color)

    # ── Candidate Metrics ───────────────────────────────────────────────────
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Candidate Metrics{Colors.ENDC}" + " " * 71 + f"{Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
    
    if info['total_candidates'] > 0:
        row("Total Candidates", f"{info['total_candidates']:,}", Colors.GREEN)
    
    if info['best_score'] is not None:
        s = info['best_score']
        if s < 0:
            row("Best Binding Score", f"{s:.2f} kcal/mol", Colors.GREEN if s < -7 else Colors.YELLOW)
        else:
            row("Best Oracle Score", f"{s:.4f}", Colors.GREEN if s > 0.5 else Colors.YELLOW)
    
    if info['best_identity'] > 0:
        row("Best Identity", f"{info['best_identity']:.1f}%", Colors.GREEN if info['best_identity'] >= 95 else Colors.YELLOW)
    
    if info['best_stability'] > 0:
        row("Best Stability (PLL)", f"{info['best_stability']:.2f}", Colors.GREEN if info['best_stability'] > 0 else Colors.YELLOW)

    # ── Timing ────────────────────────────────────────────────────────────────
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}Timing Summary{Colors.ENDC}" + " " * 71 + f"{Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
    
    if info['elapsed_time']:
        row("Total Elapsed", info['elapsed_time'], Colors.BLUE)
    
    if info['last_update']:
        row("Last Log Update", info['last_update'], Colors.BLUE)
    
    if info['pass_a_times']:
        avg_a = _avg_secs(info['pass_a_times'], 90)
        row("Pass A Avg/Scenario", f"{avg_a:.1f}s", Colors.CYAN)
    
    if info['pass_b_times']:
        avg_b = _avg_secs(info['pass_b_times'], 300)
        row("Pass B Avg/Scenario", f"{avg_b:.1f}s", Colors.CYAN)

    print(f"{Colors.BOLD}└{'─'*88}┘{Colors.ENDC}\n")

    # ── System resources ───────────────────────────────────────────────────
    resources = get_system_resources()
    
    print(f"{Colors.BOLD}┌{'─'*88}┐{Colors.ENDC}")
    print(f"{Colors.BOLD}│{Colors.ENDC} {Colors.BOLD}System Resources{Colors.ENDC}" + " " * 72 + f"{Colors.BOLD}│{Colors.ENDC}")
    print(f"{Colors.BOLD}├{'─'*88}┤{Colors.ENDC}")
    
    cpu_color = (Colors.GREEN if resources['cpu_percent'] < 70
                 else Colors.YELLOW if resources['cpu_percent'] < 90
                 else Colors.RED)
    row("CPU Usage", f"{resources['cpu_percent']:.1f}%", cpu_color)
    
    if resources['ram_total_gb'] > 0:
        ram_pct = resources['ram_used_gb'] / resources['ram_total_gb'] * 100
        ram_color = Colors.GREEN if ram_pct < 70 else Colors.YELLOW if ram_pct < 90 else Colors.RED
        row("System RAM", f"{resources['ram_used_gb']:.1f} / {resources['ram_total_gb']:.1f} GB ({ram_pct:.0f}%)", ram_color)
    
    if resources['cuda_available']:
        gpu_color = (Colors.GREEN if resources['gpu_util'] > 40
                     else Colors.YELLOW if resources['gpu_util'] > 10
                     else Colors.CYAN)
        row("GPU Compute", f"{resources['gpu_util']:.1f}% (CC {resources['gpu_compute_capability']})", gpu_color)
        
        if resources['gpu_mem_total_gb'] > 0:
            vram_pct = resources['gpu_mem_used_gb'] / resources['gpu_mem_total_gb'] * 100
            vram_color = Colors.GREEN if vram_pct < 80 else Colors.YELLOW if vram_pct < 95 else Colors.RED
            row("GPU VRAM Used", f"{resources['gpu_mem_used_gb']:.1f} / {resources['gpu_mem_total_gb']:.1f} GB ({vram_pct:.0f}%)", vram_color)
            row("PyTorch Reserved", f"{resources['gpu_memory_reserved_gb']:.2f} GB", Colors.BLUE)
    else:
        row("GPU", "Not available", Colors.RED)
    
    print(f"{Colors.BOLD}└{'─'*88}┘{Colors.ENDC}\n")

    # Warnings / errors
    if info['warnings']:
        print(f"{Colors.YELLOW}⚠️  Warnings: {len(info['warnings'])}{Colors.ENDC}")
        for w in info['warnings'][-3:]:
            print(f"   {Colors.YELLOW}•{Colors.ENDC} {w[:80]}")
        print()

    if info['errors']:
        print(f"{Colors.RED}❌ Errors: {len(info['errors'])}{Colors.ENDC}")
        for e in info['errors'][-3:]:
            print(f"   {Colors.RED}•{Colors.ENDC} {e[:80]}")
        print()


def monitor_campaign(log_path, refresh_interval=3):
    print(f"{Colors.CYAN}Monitoring: {log_path}{Colors.ENDC}\n")
    try:
        while True:
            clear_screen()
            print_header()
            print_runtime_info()
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

    print(f"{Colors.BOLD}{Colors.GREEN}p53CAD Campaign Monitor{Colors.ENDC}")
    print(f"Log: {log_path}")

    if not log_path.exists():
        print(f"{Colors.YELLOW}Log not found. Waiting for campaign...{Colors.ENDC}\n")

    monitor_campaign(log_path, refresh_interval=3)


if __name__ == "__main__":
    main()
