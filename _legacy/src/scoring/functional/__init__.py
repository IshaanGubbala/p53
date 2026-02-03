"""
Functional scoring module for p53 rescue mutations.

Evaluates rescue quality across multiple dimensions:
- Monomer stability (ΔΔG_folding)
- DNA binding affinity (ΔΔG_binding)
- Tetramer interface stability (ΔΔG_interface)
- Risk score (MSA conservation, etc.)

Two implementations available:
1. Heuristic (fast, qualitative) - binding_affinity.py, interface_stability.py
2. EvoEF2-based (slower, quantitative) - evoef2_binding.py, functional_score_evoef2.py
"""

# Heuristic implementations (fast screening)
from .binding_affinity import (
    calculate_binding_affinity,
    calculate_binding_affinity_cached,
    BindingAffinityResult
)

from .interface_stability import (
    calculate_interface_stability,
    calculate_interface_stability_cached,
    InterfaceStabilityResult
)

from .functional_score import (
    calculate_functional_score,
    score_rescue_candidates,
    FunctionalScoreResult
)

# EvoEF2-based implementations (accurate energetics)
from .evoef2_binding import (
    calculate_binding_evoef2,
    calculate_interface_evoef2,
    EvoEF2BindingResult
)

from .functional_score_evoef2 import (
    calculate_functional_score_evoef2,
    score_rescue_candidates_evoef2,
    FunctionalScoreEvoEF2Result
)

__all__ = [
    # Heuristic
    'calculate_binding_affinity',
    'calculate_binding_affinity_cached',
    'BindingAffinityResult',
    'calculate_interface_stability',
    'calculate_interface_stability_cached',
    'InterfaceStabilityResult',
    'calculate_functional_score',
    'score_rescue_candidates',
    'FunctionalScoreResult',
    # EvoEF2
    'calculate_binding_evoef2',
    'calculate_interface_evoef2',
    'EvoEF2BindingResult',
    'calculate_functional_score_evoef2',
    'score_rescue_candidates_evoef2',
    'FunctionalScoreEvoEF2Result'
]
