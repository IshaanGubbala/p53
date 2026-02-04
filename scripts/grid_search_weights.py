"""
Grid Search for Optimal Loss Function Weights
Systematically tests combinations of identity preservation and stability weights
to find the best balance between functional rescue and sequence conservation.
"""

import torch
import torch.nn.functional as F
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import product
from tqdm import tqdm

from p53cad.engine.latent import ManifoldEmbedder
from p53cad.engine.oracle import FunctionalOracle
from p53cad.data.dms import P53_WT, apply_mutation

# Grid search parameters
IDENTITY_WEIGHTS = [
    (10, 30),   # Original
    (25, 75),   # Moderate
    (50, 150),  # Strong
    (75, 225),  # Very Strong
    (100, 300), # Maximum
]

STABILITY_COEFFS = [
    3.0,   # Low
    6.0,   # Original
    9.0,   # Moderate
    12.0,  # High
    15.0,  # Very High
]

# Test configuration
TEST_MUTATION = "R175H"
TEST_STEPS = 100  # Shorter for grid search
LEARNING_RATE = 0.04

def run_optimization(embedder, oracle, identity_l1, identity_mse, stability_coeff):
    """Run a single optimization with given weights."""
    
    # Setup
    cancer_seq = apply_mutation(P53_WT, TEST_MUTATION)
    AA_IDS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    
    emb = embedder.get_embeddings(cancer_seq).detach().requires_grad_(True)
    emb_wt = embedder.get_embeddings(P53_WT).detach()
    
    optimizer = torch.optim.Adam([emb], lr=LEARNING_RATE)
    
    # Run optimization
    for step in range(TEST_STEPS):
        optimizer.zero_grad()
        
        z, logits, probs = embedder.latent_forward_ascent(emb)
        pooled = z.mean(dim=1)
        score = oracle.model(pooled)
        
        logits_aa = logits[:, :, AA_IDS]
        log_probs = F.log_softmax(logits_aa, dim=-1)
        stability = log_probs.max(dim=-1).values.mean()
        
        # Loss function with current weights
        loss = -score * 2.0
        loss -= stability_coeff * stability
        loss += identity_l1 * torch.norm(emb - emb_wt, p=1) / emb.numel()
        loss += identity_mse * F.mse_loss(emb, emb_wt)
        
        loss.backward()
        optimizer.step()
    
    # Evaluate final result
    with torch.no_grad():
        z_final, logits_final, probs_final = embedder.latent_forward_ascent(emb)
        final_score = oracle.model(z_final.mean(dim=1)).item()
        final_stability = F.log_softmax(logits_final[:, :, AA_IDS], dim=-1).max(dim=-1).values.mean().item()
        
        # Decode sequence
        top_ids_aa = torch.argmax(logits_final[:, :, AA_IDS], dim=-1)[0]
        top_ids = torch.tensor([AA_IDS[i] for i in top_ids_aa]).to(emb.device)
        tokens = embedder.tokenizer.convert_ids_to_tokens(top_ids)
        final_seq = "".join(tokens)[:len(P53_WT)]
        
        # Calculate identity
        mutations = sum(1 for i in range(len(P53_WT)) if P53_WT[i] != final_seq[i])
        identity = 100.0 * (1.0 - mutations / len(P53_WT))
        
        # Count amino acid usage
        aa_counts = {}
        for aa in final_seq:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        most_common_aa = max(aa_counts.items(), key=lambda x: x[1])
        aa_diversity = len(aa_counts)
    
    return {
        'final_score': final_score,
        'final_stability': final_stability,
        'identity_pct': identity,
        'mutation_count': mutations,
        'aa_diversity': aa_diversity,
        'most_common_aa': most_common_aa[0],
        'most_common_count': most_common_aa[1],
    }

def main():
    print("🔬 Starting Grid Search for Optimal Loss Weights")
    print(f"Testing {len(IDENTITY_WEIGHTS)} identity settings × {len(STABILITY_COEFFS)} stability settings")
    print(f"Total combinations: {len(IDENTITY_WEIGHTS) * len(STABILITY_COEFFS)}")
    print(f"Test mutation: {TEST_MUTATION}, Steps: {TEST_STEPS}\n")
    
    # Load models
    print("Loading models...")
    embedder = ManifoldEmbedder()
    oracle_path = Path("data/models/functional_oracle.pt")
    if not oracle_path.exists():
        print("❌ Error: Oracle model not found. Please train it first.")
        return
    oracle = FunctionalOracle(model_path=oracle_path)
    print("✅ Models loaded\n")
    
    # Run grid search
    results = []
    
    for (id_l1, id_mse), stab_coeff in tqdm(
        list(product(IDENTITY_WEIGHTS, STABILITY_COEFFS)),
        desc="Grid Search Progress"
    ):
        try:
            result = run_optimization(embedder, oracle, id_l1, id_mse, stab_coeff)
            results.append({
                'identity_l1_weight': id_l1,
                'identity_mse_weight': id_mse,
                'stability_coeff': stab_coeff,
                **result
            })
        except Exception as e:
            print(f"\n⚠️  Error with weights (L1={id_l1}, MSE={id_mse}, Stab={stab_coeff}): {e}")
            continue
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Save results
    output_path = Path("results/grid_search_weights.csv")
    output_path.parent.mkdir(exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"\n✅ Results saved to {output_path}")
    
    # Analysis
    print("\n" + "="*80)
    print("📊 GRID SEARCH RESULTS ANALYSIS")
    print("="*80)
    
    # Find best configurations for different objectives
    print("\n🎯 Best for Functional Rescue Score:")
    best_score = df.loc[df['final_score'].idxmax()]
    print(f"   Identity L1: {best_score['identity_l1_weight']}, MSE: {best_score['identity_mse_weight']}")
    print(f"   Stability Coeff: {best_score['stability_coeff']}")
    print(f"   → Score: {best_score['final_score']:.4f}, Identity: {best_score['identity_pct']:.1f}%")
    
    print("\n🧬 Best for Sequence Identity (>90%):")
    high_identity = df[df['identity_pct'] >= 90.0]
    if len(high_identity) > 0:
        best_identity = high_identity.loc[high_identity['final_score'].idxmax()]
        print(f"   Identity L1: {best_identity['identity_l1_weight']}, MSE: {best_identity['identity_mse_weight']}")
        print(f"   Stability Coeff: {best_identity['stability_coeff']}")
        print(f"   → Score: {best_identity['final_score']:.4f}, Identity: {best_identity['identity_pct']:.1f}%")
    else:
        print("   ⚠️  No configurations achieved >90% identity")
    
    print("\n⚖️  Best Balanced (Score × Identity):")
    df['balance_metric'] = df['final_score'] * (df['identity_pct'] / 100.0)
    best_balanced = df.loc[df['balance_metric'].idxmax()]
    print(f"   Identity L1: {best_balanced['identity_l1_weight']}, MSE: {best_balanced['identity_mse_weight']}")
    print(f"   Stability Coeff: {best_balanced['stability_coeff']}")
    print(f"   → Score: {best_balanced['final_score']:.4f}, Identity: {best_balanced['identity_pct']:.1f}%")
    print(f"   → Balance Metric: {best_balanced['balance_metric']:.4f}")
    
    print("\n🌈 Best for Amino Acid Diversity:")
    best_diversity = df.loc[df['aa_diversity'].idxmax()]
    print(f"   Identity L1: {best_diversity['identity_l1_weight']}, MSE: {best_diversity['identity_mse_weight']}")
    print(f"   Stability Coeff: {best_diversity['stability_coeff']}")
    print(f"   → AA Diversity: {best_diversity['aa_diversity']} unique amino acids")
    print(f"   → Most common: {best_diversity['most_common_aa']} ({best_diversity['most_common_count']} occurrences)")
    
    print("\n" + "="*80)
    print(f"📈 Full results available in: {output_path}")
    print("="*80)

if __name__ == "__main__":
    main()
