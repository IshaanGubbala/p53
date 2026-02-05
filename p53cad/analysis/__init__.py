"""
p53-proteoMgCAD Analysis Module

Analysis tools for rescue mutation evaluation and clinical impact assessment.
"""

from .grassmann import GrassmannMetric

# Enhanced modules
try:
    from .clinical_impact import (
        ClinicalImpactEngine,
        PatientStratifier,
        TherapeuticIndexCalculator,
        DeliveryFeasibilityAssessor,
        ImmunogenicityAssessor
    )
except ImportError:
    pass

__all__ = [
    'GrassmannMetric',
    'ClinicalImpactEngine',
    'PatientStratifier',
]
