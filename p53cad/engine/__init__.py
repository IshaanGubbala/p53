"""
p53-proteoMgCAD Engine Module

Core algorithms and analysis engines for rescue mutation design.
"""

from .latent import ManifoldEmbedder, ManifoldWalker
from .oracle import FunctionalOracle
from .explain import SaliencyMap

# Enhanced modules (may require additional dependencies)
try:
    from .drug_generator import DrugGeneratorEngine, P53_BINDING_POCKETS
except ImportError:
    pass

try:
    from .explainability import ExplainabilityEngine, P53_FUNCTIONAL_SITES
except ImportError:
    pass

try:
    from .experimental import ExperimentalPipeline
except ImportError:
    pass

try:
    from .multi_target import MultiTargetPlatform, TUMOR_SUPPRESSORS
except ImportError:
    pass

try:
    from .md_validation import MDValidationEngine
except ImportError:
    pass

__all__ = [
    'ManifoldEmbedder',
    'ManifoldWalker',
    'FunctionalOracle',
    'SaliencyMap',
    'DrugGeneratorEngine',
    'ExplainabilityEngine',
    'ExperimentalPipeline',
    'MultiTargetPlatform',
    'MDValidationEngine',
]
