from __future__ import annotations

import argparse
import yaml
from pathlib import Path
import sys
import pandas as pd

# Add project root to path
sys.path.append(str(Path.cwd()))

from src.core.logging import get_logger, setup_logging
from src.design.generative import get_generator, extract_mutations
# Reuse utilities
from experiments.run_design_rescues import _paths_from_config, _load_sequence, _parse_mutation
from experiments.run_generative_design import apply_mutation

def main():
    parser = argparse.ArgumentParser(description="Run Latent Manifold Rescue")
    parser.add_argument("--config", type=Path, default="configs/scoring.yaml")
    parser.add_argument("--targets", type=str, nargs="+", default=["R175H", "Y220C"], help="Cancer mutations to rescue")
    parser.add_argument("--n_samples", type=int, default=20, help="Number of interpolation steps")
    parser.add_argument("--model", type=str, default="facebook/esm2_t6_8M_UR50D", help="ESM model")
    
    args = parser.parse_args()
    setup_logging()
    logger = get_logger(__name__)
    
    # Load config
    if not args.config.exists():
        logger.error(f"Config not found: {args.config}")
        return 1
    with args.config.open() as f:
        config = yaml.safe_load(f)
        
    paths = _paths_from_config(config.get("paths", {}))
    
    # Load WT
    p53_cfg = config.get("p53", {})
    uniprot_id = p53_cfg.get("uniprot_id", "P04637")
    fasta_path = paths["interim"] / "uniprot" / f"{uniprot_id}.fasta"
    wt_sequence = _load_sequence(fasta_path)
    
    # Init generator
    gen_config = {
        "model_name": args.model
    }
    generator = get_generator("latent", gen_config)
    
    # Init Oracle
    from src.scoring.functional_predictor import FunctionalOracle
    oracle_path = Path("data/models/functional_oracle.pt")
    if oracle_path.exists():
        logger.info(f"Loading Functional Oracle from {oracle_path}")
        oracle = FunctionalOracle(model_path=oracle_path, input_dim=320)
    else:
        logger.warning(f"Functional Oracle not found at {oracle_path}. Skipping guided rescue.")
        oracle = None
    
    results = []
    
    for target in args.targets:
        logger.info(f"Processing Target: {target}")
        
        # 1. Create Cancer Sequence
        try:
            cancer_seq = apply_mutation(wt_sequence, target)
        except Exception as e:
            logger.error(f"Error applying mutation {target}: {e}")
            continue
            
        # 2. Run Rescue (Interpolate back to WT)
        # We pass wt_sequence as target_sequence for interpolation
        logger.info("Strategy 1: Linear Interpolation (Unsupervised)")
        candidates_linear = generator.generate(cancer_seq, n_samples=args.n_samples, target_sequence=wt_sequence)
        
        # 3. Strategy 2: Functional Steering (Guided)
        candidates_guided = []
        if oracle:
            logger.info("Strategy 2: Gradient Steering (Functional Guidance)")
            # Get latent of cancer seq
            z_cancer = generator.embedder.encode(cancer_seq)
            # Steer
            guided_path = generator.walker.steer_with_gradient(z_cancer, oracle.model, steps=args.n_samples, step_size=0.05)
            candidates_guided = guided_path
            
        # Combine
        combined_candidates = []
        for c in candidates_linear:
            combined_candidates.append((c, "Linear_Interpolation"))
        for c in candidates_guided:
            combined_candidates.append((c, "Gradient_Steering"))
            
        logger.info(f"Generated {len(combined_candidates)} total candidates.")
        
        for cand, strategy in combined_candidates:
            muts = extract_mutations(wt_sequence, cand)
            
            # Predict Score if Oracle exists
            pred_score = oracle.predict(generator.embedder.encode(cand)) if oracle else 0.0
            
            if not muts:
                cand_type = "Reversion"
            elif target in [m for m in muts]:
                if len(muts) == 1:
                    cand_type = "Self"
                else:
                    cand_type = "Rescue_Attempt"
            else:
                 cand_type = "Alternative_WT"
            
            mut_str = ",".join(sorted(muts))
            results.append({
                "target": target,
                "sequence": cand,
                "mutations": mut_str,
                "n_mutations": len(muts),
                "type": cand_type,
                "strategy": strategy,
                "predicted_function": pred_score
            })

    if not results:
        logger.warning("No results generated.")
        return 0
        
    df = pd.DataFrame(results)
    
    # Save
    out_csv = paths["processed"] / "latent_rescue_candidates.csv"
    df.to_csv(out_csv, index=False)
    logger.info(f"Saved results to {out_csv}")
    
    # Print distinct mutations
    logger.info("Distinct Candidates Types:")
    print(df["type"].value_counts())
    
    return 0

if __name__ == "__main__":
    main()
