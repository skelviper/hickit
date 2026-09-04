"""Auditable diagnostics for SNP-free diploid reconstruction.

Training-facing modules in this package never read phase or CHARM/3DG truth.
Truth-dependent evaluation is kept in the experiment drivers and marked ORACLE.
"""

from .likelihood import ScoreConfig, posterior_from_distances
from .state import STATE_LABELS, STATE_TO_INDEX

__all__ = [
    "STATE_LABELS",
    "STATE_TO_INDEX",
    "ScoreConfig",
    "posterior_from_distances",
]
