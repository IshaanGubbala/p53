import click
from pathlib import Path
import os

# Fix for OpenMP error on Mac (common with PyTorch + MKL)
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from p53cad.core.logging import setup_logging, get_logger

@click.group()
def cli():
    """p53CAD: Generative Protein Design Platform"""
    setup_logging()

@cli.command()
@click.option('--dms', type=click.Path(exists=True), default="data/raw/p53_DMS_Giacomelli_2018.csv", help="Path to DMS data")
@click.option('--epochs', default=5, help="Training epochs")
@click.option('--output', type=click.Path(), default="data/models", help="Output directory")
def train(dms, epochs, output):
    """Train the Functional Oracle on DMS data."""
    logger = get_logger("p53cad.cli.train")
    logger.info("Starting training pipeline...")
    
    from p53cad.data.dms import load_dms_data
    from p53cad.engine.latent import LatentEmbedder
    from p53cad.engine.oracle import FunctionalOracle
    
    # Logic ported from run_functional_training.py
    df = load_dms_data(dms)
    embedder = LatentEmbedder()
    oracle = FunctionalOracle(input_dim=320)
    
    # Hydrate sequences from mutation names
    if "sequence" not in df.columns:
        logger.info("Hydrating full protein sequences from mutation names...")
        from p53cad.data.dms import hydrate_sequences
        df = hydrate_sequences(df)
        logger.info(f"Hydrated {len(df)} sequences (filtered complex variants).")
        
    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)
    save_path = out_dir / "functional_oracle.pt"
    
    oracle.train(df, embedder, epochs=epochs, save_path=save_path)

@cli.command()
@click.argument('targets', nargs=-1)
@click.option('--samples', default=20, help="Number of candidates per target")
def design(targets, samples):
    """Run Latent Manifold Rescue on TARGETS (e.g. R175H)."""
    logger = get_logger("p53cad.cli.design")
    if not targets:
        logger.warning("No targets specified. Usage: p53cad design R175H Y220C")
        return

    from p53cad.engine.latent import LatentEmbedder, ManifoldWalker
    from p53cad.engine.oracle import FunctionalOracle
    import torch
    import pandas as pd
    from tqdm import tqdm
    
    logger.info(f"Designing rescue candidates for: {targets}")
    
    embedder = LatentEmbedder()
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
    from p53cad.engine.latent import LatentEmbedder
    import pandas as pd
    
    input_path = Path("data/processed/candidates.csv")
    if not input_path.exists():
        logger.error("No candidates found. Run 'p53cad design' first.")
        return
        
    df = pd.read_csv(input_path)
    embedder = LatentEmbedder()
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
    df.to_csv(input_path.with_name("candidates_analyzed.csv"), index=False)
    logger.info("Analysis complete.")

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

if __name__ == '__main__':
    cli()
