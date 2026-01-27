#!/bin/bash
# Run all 4 p53 cancer hotspot targets with functional scoring
# This will re-run the design pipeline with EvoEF2-based functional scoring

set -e  # Exit on error

echo "================================================================================"
echo "Running All Targets with Functional Scoring"
echo "================================================================================"
echo ""
echo "Targets: R175H, R248Q, R273H, Y220C"
echo "Features:"
echo "  - Monomer stability (ΔΔG_folding)"
echo "  - DNA binding affinity (ΔΔG_binding via EvoEF2 ComputeBinding)"
echo "  - Tetramer interface stability (ΔΔG_interface via EvoEF2 ComputeBinding)"
echo "  - Risk score (MSA conservation)"
echo ""
echo "Expected runtime: 3-5 hours total (using cached stability calculations)"
echo "  - R175H: ~1.5 hours (162 Pareto rescues)"
echo "  - R248Q: ~1 hour"
echo "  - R273H: ~1 hour"
echo "  - Y220C: ~1 hour"
echo ""
echo "================================================================================"
echo ""

# Create log directory
mkdir -p logs/functional_scoring

# Run each target
for TARGET in R175H R248Q R273H Y220C; do
    echo "================================================================================"
    echo "Starting: $TARGET"
    echo "================================================================================"

    LOG_FILE="logs/functional_scoring/${TARGET}_$(date +%Y%m%d_%H%M%S).log"

    echo "Running: python -m experiments.run_design_rescues --targets $TARGET --functional-scoring"
    echo "Log file: $LOG_FILE"
    echo ""

    python -m experiments.run_design_rescues \
        --targets $TARGET \
        --functional-scoring \
        2>&1 | tee "$LOG_FILE"

    if [ $? -eq 0 ]; then
        echo ""
        echo "✅ $TARGET completed successfully"
        echo ""
    else
        echo ""
        echo "❌ $TARGET failed (see log: $LOG_FILE)"
        echo ""
        exit 1
    fi
done

echo "================================================================================"
echo "All targets completed successfully!"
echo "================================================================================"
echo ""
echo "Results saved to:"
echo "  - Data/processed/rescues/R175H/pareto.parquet"
echo "  - Data/processed/rescues/R248Q/pareto.parquet"
echo "  - Data/processed/rescues/R273H/pareto.parquet"
echo "  - Data/processed/rescues/Y220C/pareto.parquet"
echo ""
echo "Each pareto.parquet file now contains functional score columns:"
echo "  - functional_score (composite 0-1 score)"
echo "  - ddg_binding (DNA binding affinity, kcal/mol)"
echo "  - ddg_interface (tetramer interface stability, kcal/mol)"
echo "  - overall_category (excellent/good/acceptable/poor)"
echo "  - binding_category, interface_category"
echo ""
echo "Next steps:"
echo "  1. Analyze results: python experiments/analyze_functional_results.py"
echo "  2. Compare old vs new rankings"
echo "  3. Generate validation report"
echo ""
