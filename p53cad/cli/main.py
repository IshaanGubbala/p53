import click
from pathlib import Path
import json
from p53cad.core.runtime import bootstrap_runtime, get_runtime_capabilities, log_runtime_capabilities, select_device

bootstrap_runtime()

from p53cad.core.logging import setup_logging, get_logger

@click.group()
def cli():
    """p53CAD: Generative Protein Design Platform"""
    log_path = setup_logging()
    logger = get_logger("p53cad.cli")
    logger.info("Workflow log file: %s", log_path)
    log_runtime_capabilities("p53cad.cli.runtime")


@cli.command()
def doctor():
    """Print runtime capability matrix and dependency readiness."""
    logger = get_logger("p53cad.cli.doctor")
    caps = get_runtime_capabilities()
    logger.info("Runtime diagnostic report:")
    for key in sorted(caps.keys()):
        logger.info("  - %s: %s", key, caps[key])

    # Optional docking diagnostics if module is available.
    try:
        from p53cad.engine.drug_generator import DrugGeneratorEngine, P53_BINDING_POCKETS

        engine = DrugGeneratorEngine()
        manifest = engine.get_receptor_manifest(allow_wt_receptor_fallback=False)
        logger.info("Drug runtime diagnostics:")
        logger.info("  - receptor_collisions: %s", manifest.get("duplicate_receptor_paths", {}))
        for pocket_key in sorted(P53_BINDING_POCKETS.keys()):
            mode = engine.get_mode_capabilities(
                pocket_key=pocket_key,
                method="docking",
                allow_wt_receptor_fallback=False,
            )
            logger.info(
                "  - %s: ready=%s receptor=%s reason=%s",
                pocket_key,
                mode.get("ready"),
                mode.get("receptor_path"),
                mode.get("reason", ""),
            )
    except Exception as exc:
        logger.warning("Drug diagnostics unavailable: %s", exc)


@cli.command("campaign-run")
@click.option("--run-id", default=None, help="Optional run identifier. Auto-generated if omitted.")
@click.option(
    "--budget",
    default="high",
    show_default=True,
    type=click.Choice(["quick", "fast", "medium", "high"]),
    help="Campaign compute budget profile (quick ~4-6h, fast ~8-12h, medium ~20h, high ~44h).",
)
@click.option("--seed", default=42, show_default=True, type=int, help="Base random seed.")
@click.option("--resume/--no-resume", default=True, show_default=True, help="Resume an existing run id if present.")
@click.option("--include-pairs/--no-include-pairs", default=True, show_default=True, help="Include hotspot pair scenarios.")
@click.option("--max-scenarios", default=None, type=int, help="Cap scenario count for dry-runs.")
@click.option("--shortlist-n", default=30, show_default=True, type=int, help="Top-N shortlist for presentation.")
@click.option(
    "--output-dir",
    default="data/campaigns",
    show_default=True,
    type=click.Path(),
    help="Campaign artifact base directory.",
)
@click.option("--no-clinical", is_flag=True, help="Skip clinical scoring pass.")
@click.option(
    "--device",
    default="auto",
    show_default=True,
    type=click.Choice(["auto", "cpu", "cuda", "mps", "xpu"]),
    help="Compute device for models (auto = CUDA > MPS > XPU > CPU).",
)
def campaign_run(
    run_id,
    budget,
    seed,
    resume,
    include_pairs,
    max_scenarios,
    shortlist_n,
    output_dir,
    no_clinical,
    device,
):
    """Run results-first campaign across hotspot scenarios and delivery methods."""
    logger = get_logger("p53cad.cli.campaign_run")
    from p53cad.engine.campaign import CampaignRunner
    from p53cad.results.store import CampaignStore

    resolved_device = str(select_device(device))
    store = CampaignStore(base_dir=Path(output_dir))
    runner = CampaignRunner(store=store, device=resolved_device)
    logger.info(
        "Starting campaign run: run_id=%s budget=%s seed=%d include_pairs=%s max_scenarios=%s device=%s",
        run_id or "auto",
        budget,
        int(seed),
        bool(include_pairs),
        max_scenarios if max_scenarios is not None else "all",
        resolved_device,
    )
    result = runner.run(
        run_id=run_id,
        seed=seed,
        resume=resume,
        include_pairs=include_pairs,
        max_scenarios=max_scenarios,
        shortlist_n=shortlist_n,
        budget=budget,
        with_clinical=not no_clinical,
    )
    logger.info(
        "Campaign complete: run_id=%s scenarios=%d candidates=%d shortlist=%d run_dir=%s",
        result.get("run_id"),
        int(result.get("n_scenarios", 0)),
        int(result.get("n_candidates", 0)),
        int(result.get("n_shortlist", 0)),
        result.get("run_dir"),
    )


@cli.command("campaign-report")
@click.option("--run-id", default=None, help="Run identifier to report. Defaults to latest run.")
@click.option("--shortlist-n", default=30, show_default=True, type=int, help="Top-N shortlist size.")
@click.option(
    "--output-dir",
    default="data/campaigns",
    show_default=True,
    type=click.Path(),
    help="Campaign artifact base directory.",
)
@click.option(
    "--oracle",
    "oracle_checkpoint",
    default=None,
    type=click.Path(exists=True),
    help=(
        "Path to an alternative oracle checkpoint (.pt) for re-scoring candidates "
        "before building the shortlist. Use this to apply a fine-tuned oracle "
        "(e.g. trained with multi-mutation augmentation) to an existing campaign "
        "without running a new optimization campaign."
    ),
)
def campaign_report(run_id, shortlist_n, output_dir, oracle_checkpoint):
    """Regenerate Top-N report artifacts from an existing campaign run.

    Optionally re-scores all candidates with a new oracle checkpoint before
    shortlist selection, allowing you to apply an improved model without
    re-running the full optimization campaign (~1 hr vs ~22 hrs).
    """
    logger = get_logger("p53cad.cli.campaign_report")
    from p53cad.engine.campaign import CampaignRunner
    from p53cad.results.store import CampaignStore

    store = CampaignStore(base_dir=Path(output_dir))
    resolved_run_id = run_id or store.latest_run_id()
    if not resolved_run_id:
        logger.error("No campaign runs found in %s", output_dir)
        return

    runner = CampaignRunner(store=store)
    report = runner.report_run(
        resolved_run_id,
        shortlist_n=shortlist_n,
        oracle_checkpoint=oracle_checkpoint,
    )
    logger.info(
        "Campaign report written: run_id=%s shortlist=%d run_dir=%s",
        report.get("run_id"),
        int(report.get("n_shortlist", 0)),
        report.get("run_dir"),
    )


@cli.command("campaign-list")
@click.option(
    "--output-dir",
    default="data/campaigns",
    show_default=True,
    type=click.Path(),
    help="Campaign artifact base directory.",
)
def campaign_list(output_dir):
    """List campaign runs and statuses from artifact index."""
    logger = get_logger("p53cad.cli.campaign_list")
    from p53cad.results.store import CampaignStore

    store = CampaignStore(base_dir=Path(output_dir))
    runs = store.list_runs()
    if runs.empty:
        logger.info("No campaign runs found in %s", output_dir)
        return
    logger.info("Campaign runs (%d):", len(runs))
    for _, row in runs.iterrows():
        logger.info(
            "  - %s | status=%s | updated=%s | path=%s",
            row.get("run_id", "n/a"),
            row.get("status", "n/a"),
            row.get("updated_at_utc", "n/a"),
            row.get("path", "n/a"),
        )

@cli.command()
@click.option('--dms', type=click.Path(exists=True), default="data/raw/p53_DMS_Giacomelli_2018.csv", help="Path to DMS data")
@click.option('--epochs', default=5, help="Training epochs")
@click.option('--output', type=click.Path(), default="data/models", help="Output directory")
@click.option('--val-split', default=0.1, show_default=True, type=float, help="Validation split fraction (0 disables).")
@click.option('--patience', default=8, show_default=True, type=int, help="Early stopping patience (validation epochs).")
@click.option('--min-delta', default=1e-4, show_default=True, type=float, help="Minimum validation improvement to reset patience.")
@click.option('--batch-size', default=32, show_default=True, type=int, help="Training batch size.")
@click.option('--seed', default=42, show_default=True, type=int, help="Random seed used for train/validation split.")
@click.option('--hidden-dim', default=128, show_default=True, type=int, help="Oracle first hidden layer width.")
@click.option('--num-layers', default=2, show_default=True, type=int, help="Number of oracle hidden layers.")
@click.option('--dropout', default=0.2, show_default=True, type=float, help="Dropout probability between hidden layers.")
@click.option('--arch', default='legacy_mlp', show_default=True, type=click.Choice(['legacy_mlp', 'attention_pooling']),
              help="Oracle architecture: legacy_mlp (mean-pooled MLP) or attention_pooling (learnable position weighting).")
def train(dms, epochs, output, val_split, patience, min_delta, batch_size, seed, hidden_dim, num_layers, dropout, arch):
    """Train the Functional Oracle on DMS data."""
    logger = get_logger("p53cad.cli.train")
    logger.info("Starting training pipeline...")
    
    from p53cad.data.dms import load_dms_data
    from p53cad.engine.latent import ManifoldEmbedder
    from p53cad.engine.oracle import FunctionalOracle
    
    # Logic ported from run_functional_training.py
    df = load_dms_data(dms)
    embedder = ManifoldEmbedder()
    oracle = FunctionalOracle(
        input_dim=embedder.hidden_size,
        hidden_dim=hidden_dim,
        num_layers=num_layers,
        dropout=dropout,
        arch=arch,
    )
    
    # Hydrate sequences from mutation names
    if "sequence" not in df.columns:
        logger.info("Hydrating full protein sequences from mutation names...")
        from p53cad.data.dms import hydrate_sequences
        df = hydrate_sequences(df)
        logger.info(f"Hydrated {len(df)} sequences (filtered complex variants).")
        
    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)
    save_path = out_dir / "functional_oracle.pt"
    
    oracle.train(
        df,
        embedder,
        epochs=epochs,
        save_path=save_path,
        val_split=val_split,
        early_stopping_patience=patience,
        min_delta=min_delta,
        batch_size=batch_size,
        seed=seed,
    )

@cli.command("train-multimut")
@click.option('--dms', type=click.Path(exists=True), default="data/raw/p53_DMS_Giacomelli_2018.csv",
              show_default=True, help="Path to DMS data (same as used for original training).")
@click.option('--checkpoint', type=click.Path(exists=True), default="data/models/functional_oracle.pt",
              show_default=True, help="Existing attention_pooling oracle checkpoint to fine-tune.")
@click.option('--output', type=click.Path(), default="data/models/functional_oracle_multimut.pt",
              show_default=True, help="Output path for fine-tuned checkpoint.")
@click.option('--n-synthetic', default=50_000, show_default=True, type=int,
              help="Number of synthetic k-mutation combinations to generate.")
@click.option('--k-min', default=2, show_default=True, type=int, help="Minimum mutations per synthetic example.")
@click.option('--k-max', default=15, show_default=True, type=int, help="Maximum mutations per synthetic example.")
@click.option('--epochs', default=10, show_default=True, type=int, help="Fine-tuning epochs.")
@click.option('--lr', default=2e-4, show_default=True, type=float,
              help="Learning rate (lower than original to prevent catastrophic forgetting).")
@click.option('--batch-size', default=32, show_default=True, type=int)
@click.option('--val-split', default=0.1, show_default=True, type=float)
@click.option('--patience', default=5, show_default=True, type=int, help="Early stopping patience.")
@click.option('--seed', default=42, show_default=True, type=int)
def train_multimut(dms, checkpoint, output, n_synthetic, k_min, k_max, epochs, lr, batch_size, val_split, patience, seed):
    """Fine-tune oracle on synthetic multi-mutation data to fix single-mutation extrapolation.

    The base oracle is trained exclusively on single-residue DMS variants, so during
    inference on 18-30-mutation rescue candidates the attention mechanism must
    extrapolate far outside its training distribution.

    This command generates N synthetic k-mutation combinations using the thermodynamic
    additivity pseudo-label model (Z_pseudo = additive sum with soft saturation), mixes
    them with the original single-mutation data, and fine-tunes the oracle at a reduced
    learning rate to prevent catastrophic forgetting.

    Literature-validated multi-mutation rescues (N239Y+R249S, H168R+R249S, T284R+R175H,
    etc.) are injected as high-confidence anchors (5x upweighted).

    Example:
        p53cad train-multimut --n-synthetic 50000 --k-max 20 --epochs 15
    """
    logger = get_logger("p53cad.cli.train_multimut")
    from p53cad.data.dms import load_dms_data, hydrate_sequences
    from p53cad.engine.latent import ManifoldEmbedder
    from p53cad.engine.oracle import FunctionalOracle

    logger.info("Loading DMS data from %s", dms)
    df = load_dms_data(dms)
    if "sequence" not in df.columns:
        logger.info("Hydrating sequences from mutation names...")
        df = hydrate_sequences(df)
        logger.info("Hydrated %d sequences.", len(df))

    logger.info("Loading embedder (ESM-2 650M)...")
    embedder = ManifoldEmbedder()

    logger.info("Loading oracle checkpoint: %s", checkpoint)
    oracle = FunctionalOracle(
        model_path=checkpoint,
        input_dim=embedder.hidden_size,
        arch="attention_pooling",
    )

    if oracle.arch_name != "attention_pooling":
        raise click.ClickException(
            "train-multimut requires an attention_pooling oracle. "
            f"Loaded checkpoint has arch='{oracle.arch_name}'. "
            "Re-train with: p53cad train --arch attention_pooling"
        )

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    oracle.train_multimut(
        dms_data=df,
        embedder=embedder,
        n_synthetic=n_synthetic,
        k_range=(k_min, k_max),
        epochs=epochs,
        save_path=output_path,
        val_split=val_split,
        early_stopping_patience=patience,
        batch_size=batch_size,
        lr=lr,
        seed=seed,
    )
    logger.info("Multi-mutation fine-tuned oracle saved to %s", output_path)
    logger.info(
        "To use: p53cad campaign-run --oracle %s  (or copy to data/models/functional_oracle.pt)",
        output_path,
    )


@cli.command()
@click.argument('targets', nargs=-1)
@click.option('--samples', default=20, help="Number of candidates per target")
@click.option('--lock', default="", help="Residue positions to lock (e.g. 248,273)")
def design(targets, samples, lock):
    """Run Latent Manifold Rescue on TARGETS (e.g. R175H)."""
    logger = get_logger("p53cad.cli.design")
    
    # Parse locked residues
    locked_indices = []
    if lock:
        try:
            locked_indices = [int(x.strip()) - 1 for x in lock.split(",") if x.strip()]
            logger.info(f"Locking critical residues: {[i+1 for i in locked_indices]}")
        except ValueError:
            logger.error("Invalid lock format. Use comma-separated integers.")
    if not targets:
        logger.warning("No targets specified. Usage: p53cad design R175H Y220C")
        return

    from p53cad.engine.latent import ManifoldEmbedder, ManifoldWalker
    from p53cad.engine.oracle import FunctionalOracle
    import torch
    import pandas as pd
    from tqdm import tqdm
    
    logger.info(f"Designing rescue candidates for: {targets}")
    
    embedder = ManifoldEmbedder()
    walker = ManifoldWalker(embedder)
    
    # Load Oracle
    oracle_path = Path("data/models/functional_oracle.pt")
    oracle = None
    if oracle_path.exists():
        oracle = FunctionalOracle(model_path=oracle_path)
    else:
        logger.warning("Functional Oracle not found. Running in unsupervised mode.")
    
    results = []
    from p53cad.data.dms import P53_WT
    wt_seq = P53_WT 
    
    for target in targets:
        logger.info(f"Target: {target}")
        from p53cad.data.dms import apply_mutation
        cancer_seq = apply_mutation(wt_seq, target)
        if not cancer_seq:
            logger.warning(f"Could not parse target mutation: {target}. Using WT baseline.")
            cancer_seq = wt_seq
            
        z_cancer = embedder.encode(cancer_seq)
        
        # Guided Steering with Live Feedback
        if oracle:
            logger.info("Running Gradient Steering...")
            # We implement the loop here manually to show progress
            # z is (1, L, D)
            z = z_cancer.clone().detach().requires_grad_(True)
            optimizer = torch.optim.Adam([z], lr=0.05)
            
            pbar = tqdm(range(samples), desc="Steering Latent Vector")
            
            path_guided = []
            for i in pbar:
                optimizer.zero_grad()
                pooled = z.mean(dim=1)
                score = oracle.model(pooled)
                loss = -score
                
                # Constraint: Discourage mutations in locked residues
                # We can do this by adding a penalty if the latent vector drifts 
                # too far from WT in those specific positions
                if locked_indices:
                    z_wt = embedder.encode(wt_seq)
                    lock_loss = torch.norm(z[:, locked_indices, :] - z_wt[:, locked_indices, :])
                    loss += 10.0 * lock_loss # Strong penalty
                
                loss.backward()
                optimizer.step()
                
                # Identify current mutation (fast check)
                with torch.no_grad():
                     seq = embedder.decode(z)
                     path_guided.append(seq)
                     
                     # Simple logic to find mutation relative to WT
                     # We assume a reference WT sequence is available
                     from p53cad.engine.latent import ManifoldWalker
                     # For more detailed feedback, we could use extract_mutations if we had it
                     # But for 'live' display, just showing the current score is good.
                     # Let's add a custom detection:
                     if len(seq) == len(wt_seq):
                         muts = [f"{wt_seq[j]}{j+1}{seq[j]}" for j in range(len(wt_seq)) if wt_seq[j] != seq[j]]
                         mut_str = ",".join(muts[:2]) + ("..." if len(muts) > 2 else "")
                     else:
                         mut_str = "Aligning..."
                
                # Update progress bar
                current_score = score.item()
                pbar.set_postfix({"FuncScore": f"{current_score:.4f}", "Mut": mut_str})
                     
            for p in path_guided:
                muts = [f"{wt_seq[j]}{j+1}{p[j]}" for j in range(len(wt_seq)) if wt_seq[j] != p[j]]
                results.append({"target": target, "sequence": p, "mutations": ",".join(muts), "strategy": "guided"})
        
        # Linear Interpolation (Fast)
        logger.info("Generating interpolation baseline...")
        path = walker.interpolate(cancer_seq, wt_seq, steps=samples)
        for p in path:
            muts = [f"{wt_seq[j]}{j+1}{p[j]}" for j in range(len(wt_seq)) if wt_seq[j] != p[j]]
            results.append({"target": target, "sequence": p, "mutations": ",".join(muts), "strategy": "linear"})
            
    # Save
    out_path = Path("data/processed/candidates.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(out_path, index=False)
    logger.info(f"Saved {len(results)} candidates to {out_path}")

@cli.command()
def analyze():
    """Run Grassmannian Analysis."""
    logger = get_logger("p53cad.cli.analyze")
    logger.info("Running Grassmannian Analysis...")
    from p53cad.analysis.grassmann import GrassmannMetric
    from p53cad.engine.latent import ManifoldEmbedder
    import pandas as pd
    
    input_path = Path("data/processed/candidates.csv")
    if not input_path.exists():
        logger.error("No candidates found. Run 'p53cad design' first.")
        return
        
    df = pd.read_csv(input_path)
    embedder = ManifoldEmbedder()
    metric = GrassmannMetric(embedder)
    
    # Mock analysis loop
    dists = []
    from p53cad.data.dms import P53_WT
    wt_seq = P53_WT
    for seq in df["sequence"]:
        try:
            d = metric.grassmann_distance(wt_seq, seq)
            dists.append(d)
        except Exception as e:
            logger.warning(f"Error: {e}")
            dists.append(0.0)
            
    df["grassmann_dist"] = dists
    
    # Pareto Ranking: Multi-objective score
    # Normalize Grassmann Distance (0-1)
    if not df.empty and df["grassmann_dist"].max() > df["grassmann_dist"].min():
        norm_dist = (df["grassmann_dist"] - df["grassmann_dist"].min()) / (df["grassmann_dist"].max() - df["grassmann_dist"].min())
        # Higher score is better, Lower dist is better -> score + (1-dist)
        # We need Predicted Score here too if available
        # But analyze is generic. Let's stick to Grassmann dist for now or add a combined score.
        logger.info("Ranking candidates based on physics-informed novelty...")
        
    df.to_csv(input_path.with_name("candidates_analyzed.csv"), index=False)
    logger.info("Analysis complete.")


@cli.command()
@click.option("--pocket", default="Y220C_cavity", show_default=True, help="Target pocket key.")
@click.option("--n-candidates", default=10, show_default=True, type=int, help="Number of drug candidates.")
@click.option("--method", default="template", show_default=True, type=click.Choice(["template", "docking", "docking_md"]), help="Drug generation mode.")
@click.option("--allow-wt-receptor-fallback", is_flag=True, help="Allow fallback to shared WT receptor if pocket-specific receptor is missing.")
@click.option("--output", type=click.Path(), default="data/processed/drug_candidates.json", show_default=True, help="Output JSON path.")
def drug(pocket, n_candidates, method, allow_wt_receptor_fallback, output):
    """Generate drug candidates with strict capability checks."""
    logger = get_logger("p53cad.cli.drug")
    from p53cad.engine.drug_generator import DrugGeneratorEngine

    engine = DrugGeneratorEngine(allow_wt_receptor_fallback=allow_wt_receptor_fallback)
    caps = engine.get_mode_capabilities(
        pocket_key=pocket,
        method=method,
        allow_wt_receptor_fallback=allow_wt_receptor_fallback,
    )
    if not caps.get("ready", False):
        logger.error("Drug mode unavailable: %s", caps.get("reason", "unknown reason"))
        return

    candidates = engine.generate_for_pocket(
        pocket_name=pocket,
        n_candidates=n_candidates,
        method=method,
        allow_wt_receptor_fallback=allow_wt_receptor_fallback,
    )
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "pocket": pocket,
        "method": method,
        "allow_wt_receptor_fallback": bool(allow_wt_receptor_fallback),
        "count": len(candidates),
        "candidates": [c.to_dict() for c in candidates],
    }
    out_path.write_text(json.dumps(payload, indent=2))
    logger.info("Wrote %d candidates to %s", len(candidates), out_path)

@cli.command()
@click.option('--pdb', type=click.Path(exists=True), default="data/raw/p53_wt.pdb", help="Path to WT PDB")
@click.option('--csv', type=click.Path(exists=True), default="data/processed/candidates.csv", help="Candidates CSV")
@click.option('--output', type=click.Path(), default="p53_rescue_session.pml", help="Output PyMol script")
def visualize(pdb, csv, output):
    """Generate PyMol session to visualize candidates."""
    logger = get_logger("p53cad.cli.visualize")
    from p53cad.viz.pymol import PyMolGenerator
    
    logger.info(f"Generating visualization script for {csv}...")
    
    viz = PyMolGenerator(pdb_path=pdb)
    viz.generate_session(Path(csv), Path(output))
    
    logger.info(f"Done! Launching PyMol...")
    import subprocess
    try:
        subprocess.run(["pymol", str(output)], check=True)
    except FileNotFoundError:
        logger.error("PyMol not found in PATH. Please run 'p53cad visualize' when PyMol is installed or available.")
    except Exception as e:
        logger.error(f"Failed to launch PyMol: {e}")

@cli.command()
def explain():
    """Run Explainable AI (Saliency) on top candidates."""
    logger = get_logger("p53cad.cli.explain")
    from p53cad.engine.explain import SaliencyMap
    from p53cad.engine.oracle import FunctionalOracle
    from p53cad.engine.latent import ManifoldEmbedder
    from p53cad.data.dms import P53_WT
    
    oracle_path = Path("data/models/functional_oracle.pt")
    if not oracle_path.exists():
        logger.error("Oracle not found. Run 'p53cad train' first.")
        return
        
    embedder = ManifoldEmbedder()
    oracle = FunctionalOracle(model_path=oracle_path)
    explainer = SaliencyMap(oracle, embedder)
    
    logger.info("Analyzing p53 Wild-Type saliency...")
    hotspots = explainer.get_top_hotspots(P53_WT)
    for h in hotspots:
        logger.info(f"Position {h['pos']}{h['aa']}: Importance {h['importance']:.4f}")

@cli.command()
@click.option("--run-id", default=None, help="Campaign run identifier. Defaults to latest run.")
@click.option("--top-n", default=None, type=int, help="Only validate top N candidates (by oracle rank). Default: all.")
@click.option("--tier2-top-n", default=3, show_default=True, type=int, help="Number of top candidates for Tier 2 MD simulation.")
@click.option("--simulation-ns", default=0.2, show_default=True, type=float, help="MD production length in nanoseconds.")
@click.option("--skip-md", is_flag=True, help="Skip Tier 2 MD stability simulations (fastest mode).")
@click.option("--skip-esmfold", is_flag=True, help="Skip local ESMFold structure prediction.")
@click.option("--skip-energy", is_flag=True, help="Skip OpenMM energy minimization.")
@click.option("--skip-dna", is_flag=True, help="Skip DNA-binding interface analysis.")
@click.option("--skip-binding-sim", is_flag=True, help="Skip MM-GBSA DNA-binding simulation (expensive).")
@click.option("--device", default="auto", show_default=True, type=click.Choice(["auto", "cpu", "cuda", "mps"]), help="Device for ESMFold inference (auto = CUDA > MPS > CPU).")
@click.option(
    "--output-dir",
    default="data/campaigns",
    show_default=True,
    type=click.Path(),
    help="Campaign artifact base directory.",
)
def validate(run_id, top_n, tier2_top_n, simulation_ns, skip_md, skip_esmfold, skip_energy, skip_dna, skip_binding_sim, device, output_dir):
    """Run physics-based validation on a campaign's top candidates."""
    logger = get_logger("p53cad.cli.validate")
    from p53cad.results.store import CampaignStore

    store = CampaignStore(base_dir=Path(output_dir))
    resolved_run_id = run_id or store.latest_run_id()
    if not resolved_run_id:
        logger.error("No campaign runs found in %s", output_dir)
        return

    bundle = store.load_run_bundle(resolved_run_id)
    top_df = bundle["top30"]
    if top_df.empty:
        logger.error("No top-30 candidates found for run %s", resolved_run_id)
        return

    device = str(select_device(device))
    logger.info("Physics validation: run_id=%s candidates=%d device=%s", resolved_run_id, len(top_df), device)

    from p53cad.engine.physics_validation import PhysicsValidationPipeline
    from p53cad.data.dms import P53_WT, apply_mutation

    # Build cancer sequence map from targets in top_df
    cancer_sequences = {}
    if "target_label" in top_df.columns:
        for target_label in top_df["target_label"].unique():
            # Extract first target mutation from the label (e.g. "R175H+gene_therapy" -> "R175H")
            target_mut = str(target_label).split("+")[0].strip()
            cancer_seq = apply_mutation(P53_WT, target_mut)
            if cancer_seq:
                cancer_sequences[str(target_label)] = cancer_seq

    run_dir = bundle["run_dir"]
    pipeline = PhysicsValidationPipeline(
        device=device,
        cache_dir=Path(run_dir) / "esmfold_cache",
        wt_pdb_path=Path("data/raw/p53_wt.pdb"),
    )
    report = pipeline.validate_campaign(
        run_id=resolved_run_id,
        top_df=top_df,
        output_dir=Path(run_dir),
        wt_sequence=P53_WT,
        cancer_sequences=cancer_sequences,
        tier1_top_n=top_n,
        tier2_top_n=tier2_top_n,
        simulation_ns=simulation_ns,
        skip_esmfold=skip_esmfold,
        skip_energy=skip_energy,
        skip_md=skip_md,
        skip_dna=skip_dna,
        skip_binding_sim=skip_binding_sim,
    )

    # Save report
    val_path = store.write_validation(resolved_run_id, "physics_validation.json", report.to_dict())
    logger.info("Physics validation saved to %s", val_path)

    # Print summary table
    has_binding_sim = any(c.get("dna_binding_sim") for c in report.candidates)
    logger.info("Physics Validation Summary (%d candidates, %.1f min):",
                report.n_candidates, report.elapsed_total_sec / 60)
    if has_binding_sim:
        logger.info("  %-4s %-20s %-8s %-8s %-10s %-10s %-8s %-8s %s",
                    "Rank", "Target", "pLDDT", "DDG(WT)", "DNA Geom", "Bind dG", "H-bonds", "Score", "Verdict")
    else:
        logger.info("  %-4s %-20s %-8s %-8s %-10s %-8s %s",
                    "Rank", "Target", "pLDDT", "DDG(WT)", "DNA Score", "Score", "Verdict")
    for c in report.candidates:
        ef = c.get("esmfold", {})
        ddg_info = c.get("ddg", {})
        dna = c.get("dna_binding", {})
        bsim = c.get("dna_binding_sim", {})
        if has_binding_sim:
            logger.info("  %-4d %-20s %-8s %-8s %-10s %-10s %-8s %-8s %s",
                         c.get("candidate_rank", 0),
                         c.get("target_label", "?")[:20],
                         f"{ef.get('dbd_plddt', 0):.1f}" if ef else "-",
                         f"{ddg_info.get('ddg_vs_wt_kcal', 0):+.1f}" if ddg_info else "-",
                         f"{dna.get('interface_preservation_score', 0):.2f}" if dna else "-",
                         f"{bsim.get('delta_binding_kcal', 0):+.1f}" if bsim else "-",
                         f"{bsim.get('n_hbonds_protein_dna', 0)}" if bsim else "-",
                         f"{c.get('overall_physics_score', 0):.0f}",
                         c.get("verdict", "?"))
        else:
            logger.info("  %-4d %-20s %-8s %-8s %-10s %-8s %s",
                         c.get("candidate_rank", 0),
                         c.get("target_label", "?")[:20],
                         f"{ef.get('dbd_plddt', 0):.1f}" if ef else "-",
                         f"{ddg_info.get('ddg_vs_wt_kcal', 0):+.1f}" if ddg_info else "-",
                         f"{dna.get('interface_preservation_score', 0):.2f}" if dna else "-",
                         f"{c.get('overall_physics_score', 0):.0f}",
                         c.get("verdict", "?"))


@cli.command()
def lab():
    """Launch the p53CAD Generative Lab Dashboard."""
    logger = get_logger("p53cad.cli.lab")
    logger.info("Starting Bio-CAD Laboratory...")
    import subprocess
    app_path = Path(__file__).parent.parent / "app" / "main.py"
    subprocess.run(["streamlit", "run", str(app_path)])

if __name__ == '__main__':
    cli()
