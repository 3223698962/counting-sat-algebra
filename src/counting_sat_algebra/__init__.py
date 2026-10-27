"""counting_sat_algebra: small algebraic #SAT testbed (Steps 0, 1, and 2).

Steps 0 and 1 check the counting polynomial, univariate interpolation, and a
bounded-occurrence degree cap. Step 2 checks a simplified (r, τ, D)-design and
a blocked sum against the full sum H(t) on hand-built replacements.

Every check still enumerates an exponential set: the cube {0,1}^n, or a
measurement space that is exponential in the number of blocks for the popcount
design. The matrix M0 and Theorem 1.1 are not implemented. Nothing here runs
in 2^{o(n)} time.
"""

from .blocked_sum import (
    BlockedExample,
    BlockedSumResult,
    GrowthRow,
    StructuredReplacement,
    blocked_sum,
    hand_built_examples,
    hand_built_replacement,
    popcount_group_instance,
    popcount_growth,
)
from .design import BlockDesign, make_design, popcount_design

__all__ = [
    "BlockDesign",
    "BlockedExample",
    "BlockedSumResult",
    "GrowthRow",
    "StructuredReplacement",
    "blocked_sum",
    "hand_built_examples",
    "hand_built_replacement",
    "make_design",
    "popcount_design",
    "popcount_group_instance",
    "popcount_growth",
]

__version__ = "0.3.0"
