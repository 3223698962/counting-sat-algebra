"""counting_sat_algebra: small algebraic #SAT testbed (Steps 0 through 6).

Steps 0 and 1 check the counting polynomial, univariate interpolation, and a
bounded-occurrence degree cap. Step 2 checks a simplified (r, τ, D)-design and
a blocked sum against the full sum H(t) on hand-built replacements. Step 3
builds the replacement matrix M0 = V0 f_I(W0) only at side lengths up to 64.
Step 4 recovers H(0) from Hasse coefficients at fewer distinct sample points
than the D+1 value samples of Lagrange interpolation. The condition count
stays D+1, and each sample point still enumerates the cube. Step 5 copies
variables through an implication cycle so each occurs in at most three
clauses, with #F' = #F on instances small enough to enumerate. Step 6 scans
a machine-limited grid and records which of those pieces match brute force
and which stop on an unfinished construction or a budget.

Every check still enumerates an exponential set: the cube {0,1}^n, a
measurement space that is exponential in the number of blocks for the popcount
design, or a matrix whose side length is r**(2Q) on a tiny parameter. Copying
can make the reduced n larger than the source n, and past that point the
count is recorded as unchecked. The matrix at the parameter list R ≥ 32,
6R ≤ Q is not built. The scan does not time that matrix. Theorem 1.1 is not
implemented. Nothing here runs in 2^{o(n)} time.
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
from .hermite import (
    HermiteRecovery,
    Step4Case,
    analyze_step4_case,
    eval_H_hasse_jet,
    hermite_eval,
    step4_cases,
)
from .replacement_matrix import (
    AuxiliarySlice,
    PlacedFormula,
    Step3Example,
    build_auxiliary_slice,
    place_formula,
    step3_examples,
)
from .parameter_scan import (
    DISTANCE_TO_FULL_ALGORITHM,
    TIMING_DISCLAIMER,
    ScanCell,
    ScanReport,
    run_parameter_scan,
    scan_specs,
)
from .sparsify import (
    Sparsification,
    Step5Case,
    Step5Report,
    analyze_step5_case,
    sparsify,
    step5_cases,
)

__all__ = [
    "AuxiliarySlice",
    "BlockDesign",
    "BlockedExample",
    "BlockedSumResult",
    "GrowthRow",
    "HermiteRecovery",
    "PlacedFormula",
    "Step3Example",
    "Sparsification",
    "Step4Case",
    "DISTANCE_TO_FULL_ALGORITHM",
    "ScanCell",
    "ScanReport",
    "Step5Case",
    "Step5Report",
    "StructuredReplacement",
    "TIMING_DISCLAIMER",
    "analyze_step4_case",
    "analyze_step5_case",
    "blocked_sum",
    "build_auxiliary_slice",
    "eval_H_hasse_jet",
    "hand_built_examples",
    "hand_built_replacement",
    "hermite_eval",
    "make_design",
    "place_formula",
    "popcount_design",
    "run_parameter_scan",
    "popcount_group_instance",
    "popcount_growth",
    "scan_specs",
    "sparsify",
    "step3_examples",
    "step4_cases",
    "step5_cases",
]

__version__ = "0.7.0"
