"""
Analyze functional scoring results across all 4 targets.

Compares:
- Old Pareto fronts (stability-only)
- New Pareto fronts (functional scoring with DNA + interface)

Reports:
- Filtering statistics (how many rescues fail DNA/interface)
- Category distribution (excellent/good/acceptable/poor)
- Ranking changes (which rescues improved/declined)
- Top candidates per target
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json


def load_target_results(target: str) -> dict:
    """Load both old and new results for a target."""
    base_path = Path("Data/processed/rescues") / target
    pareto_path = base_path / "pareto.parquet"

    if not pareto_path.exists():
        return None

    pareto_df = pd.read_parquet(pareto_path)

    # Check if functional scoring was applied
    has_functional = 'functional_score' in pareto_df.columns

    return {
        'target': target,
        'pareto_df': pareto_df,
        'has_functional': has_functional,
        'n_rescues': len(pareto_df)
    }


def analyze_target(data: dict) -> dict:
    """Analyze functional scoring results for one target."""
    if not data or not data['has_functional']:
        return None

    df = data['pareto_df']
    target = data['target']

    # Category distribution
    categories = df['overall_category'].value_counts().to_dict()

    # DNA binding analysis
    binding_good = len(df[df['binding_category'] == 'good'])
    binding_acceptable = len(df[df['binding_category'] == 'acceptable'])
    binding_bad = len(df[df['binding_category'] == 'bad'])

    # Interface analysis
    interface_good = len(df[df['interface_category'] == 'good'])
    interface_acceptable = len(df[df['interface_category'] == 'acceptable'])
    interface_bad = len(df[df['interface_category'] == 'bad'])

    # Top rescues
    top_10 = df.nlargest(10, 'functional_score')[
        ['rescue_mutations', 'functional_score', 'overall_category',
         'ddg_folding', 'ddg_binding', 'ddg_interface', 'risk']
    ].to_dict('records')

    # Failures (rescues that were filtered out)
    failures = df[df['overall_category'] == 'poor']

    return {
        'target': target,
        'total_rescues': len(df),
        'categories': categories,
        'binding': {
            'good': binding_good,
            'acceptable': binding_acceptable,
            'bad': binding_bad,
            'pass_rate': (binding_good + binding_acceptable) / len(df) * 100
        },
        'interface': {
            'good': interface_good,
            'acceptable': interface_acceptable,
            'bad': interface_bad,
            'pass_rate': (interface_good + interface_acceptable) / len(df) * 100
        },
        'top_10': top_10,
        'failures': failures[['rescue_mutations', 'ddg_folding', 'ddg_binding',
                              'ddg_interface', 'binding_category', 'interface_category']
                            ].to_dict('records') if len(failures) > 0 else []
    }


def print_target_report(analysis: dict):
    """Print detailed report for one target."""
    if not analysis:
        return

    target = analysis['target']

    print("=" * 80)
    print(f"{target} - Functional Scoring Results")
    print("=" * 80)
    print()

    # Overview
    print(f"Total Pareto rescues: {analysis['total_rescues']}")
    print()

    # Category distribution
    print("Category Distribution:")
    for cat in ['excellent', 'good', 'acceptable', 'poor']:
        count = analysis['categories'].get(cat, 0)
        pct = count / analysis['total_rescues'] * 100
        print(f"  {cat.capitalize():<12}: {count:3d} ({pct:5.1f}%)")
    print()

    # DNA binding filter
    print("DNA Binding Filter:")
    binding = analysis['binding']
    print(f"  Good:       {binding['good']:3d} ({binding['good']/analysis['total_rescues']*100:5.1f}%)")
    print(f"  Acceptable: {binding['acceptable']:3d} ({binding['acceptable']/analysis['total_rescues']*100:5.1f}%)")
    print(f"  Bad:        {binding['bad']:3d} ({binding['bad']/analysis['total_rescues']*100:5.1f}%)")
    print(f"  Pass rate:  {binding['pass_rate']:.1f}%")
    print()

    # Interface filter
    print("Tetramer Interface Filter:")
    interface = analysis['interface']
    print(f"  Good:       {interface['good']:3d} ({interface['good']/analysis['total_rescues']*100:5.1f}%)")
    print(f"  Acceptable: {interface['acceptable']:3d} ({interface['acceptable']/analysis['total_rescues']*100:5.1f}%)")
    print(f"  Bad:        {interface['bad']:3d} ({interface['bad']/analysis['total_rescues']*100:5.1f}%)")
    print(f"  Pass rate:  {interface['pass_rate']:.1f}%")
    print()

    # Top 10
    print("Top 10 Functional Rescues:")
    print(f"{'Rank':<6}{'Rescue':<15}{'Func':<8}{'Cat':<12}{'ΔΔG_fold':<11}{'ΔΔG_bind':<11}{'ΔΔG_int':<11}")
    print("-" * 80)
    for i, rescue in enumerate(analysis['top_10'], 1):
        print(f"{i:<6}"
              f"{rescue['rescue_mutations']:<15}"
              f"{rescue['functional_score']:<8.3f}"
              f"{rescue['overall_category']:<12}"
              f"{rescue['ddg_folding']:<11.2f}"
              f"{rescue['ddg_binding']:<11.2f}"
              f"{rescue['ddg_interface']:<11.2f}")
    print()

    # Failures
    if analysis['failures']:
        print(f"⚠️  Rescues Filtered Out (category: poor) - {len(analysis['failures'])} total:")
        print()
        for failure in analysis['failures'][:5]:  # Show first 5
            print(f"  - {failure['rescue_mutations']}")
            print(f"      ΔΔG_binding = {failure['ddg_binding']:.2f} ({failure['binding_category']})")
            print(f"      ΔΔG_interface = {failure['ddg_interface']:.2f} ({failure['interface_category']})")
        if len(analysis['failures']) > 5:
            print(f"  ... and {len(analysis['failures']) - 5} more")
        print()
    else:
        print("✅ All rescues pass functional filters!")
        print()


def main():
    print()
    print("=" * 80)
    print("Functional Scoring Analysis - All Targets")
    print("=" * 80)
    print()

    targets = ['R175H', 'R248Q', 'R273H', 'Y220C']

    # Load all results
    all_data = {}
    for target in targets:
        data = load_target_results(target)
        if data and data['has_functional']:
            all_data[target] = data
        elif data:
            print(f"⚠️  {target}: No functional scoring found (old results)")
        else:
            print(f"❌ {target}: No results found")

    if not all_data:
        print()
        print("No functional scoring results found.")
        print("Run: bash experiments/run_all_functional.sh")
        return

    print(f"Found functional scoring for {len(all_data)} targets")
    print()

    # Analyze each target
    analyses = {}
    for target, data in all_data.items():
        analyses[target] = analyze_target(data)

    # Print individual reports
    for target in targets:
        if target in analyses:
            print_target_report(analyses[target])

    # Summary across all targets
    print("=" * 80)
    print("Summary Across All Targets")
    print("=" * 80)
    print()

    total_rescues = sum(a['total_rescues'] for a in analyses.values())
    total_excellent = sum(a['categories'].get('excellent', 0) for a in analyses.values())
    total_good = sum(a['categories'].get('good', 0) for a in analyses.values())
    total_acceptable = sum(a['categories'].get('acceptable', 0) for a in analyses.values())
    total_poor = sum(a['categories'].get('poor', 0) for a in analyses.values())

    print(f"Total Pareto rescues across all targets: {total_rescues}")
    print()
    print("Overall Category Distribution:")
    print(f"  Excellent:  {total_excellent:3d} ({total_excellent/total_rescues*100:5.1f}%)")
    print(f"  Good:       {total_good:3d} ({total_good/total_rescues*100:5.1f}%)")
    print(f"  Acceptable: {total_acceptable:3d} ({total_acceptable/total_rescues*100:5.1f}%)")
    print(f"  Poor:       {total_poor:3d} ({total_poor/total_rescues*100:5.1f}%)")
    print()

    # DNA binding summary
    total_binding_pass = sum(
        a['binding']['good'] + a['binding']['acceptable']
        for a in analyses.values()
    )
    print(f"DNA Binding Pass Rate: {total_binding_pass/total_rescues*100:.1f}%")

    # Interface summary
    total_interface_pass = sum(
        a['interface']['good'] + a['interface']['acceptable']
        for a in analyses.values()
    )
    print(f"Tetramer Interface Pass Rate: {total_interface_pass/total_rescues*100:.1f}%")
    print()

    # Save summary JSON
    summary_path = Path("Data/processed/functional_scores/summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        'targets': list(analyses.keys()),
        'total_rescues': total_rescues,
        'category_distribution': {
            'excellent': total_excellent,
            'good': total_good,
            'acceptable': total_acceptable,
            'poor': total_poor
        },
        'pass_rates': {
            'dna_binding': total_binding_pass / total_rescues * 100,
            'interface': total_interface_pass / total_rescues * 100
        },
        'per_target': {
            target: {
                'total': a['total_rescues'],
                'excellent': a['categories'].get('excellent', 0),
                'dna_pass_rate': a['binding']['pass_rate'],
                'interface_pass_rate': a['interface']['pass_rate']
            }
            for target, a in analyses.items()
        }
    }

    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"Summary saved to: {summary_path}")
    print()


if __name__ == "__main__":
    main()
