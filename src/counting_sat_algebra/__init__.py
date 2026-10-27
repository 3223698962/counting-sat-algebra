"""Exact algebraic #SAT components aligned with the local reference.

Steps 0–8 preserve the original finite checks. Steps 9–11 implement the
reference's conditional-expectation schedules, exact rational design checks,
streamed shared-measurement Taylor recovery, and structured auxiliary entry
queries even at R=32, Q=192. Step 12 connects those entries to sparse admitted
placements. Step 13 counts independent components and supports DIMACS input.
Step 14 contracts disjoint supports of the actual replaced factors. Step 15
uses finite Taylor tensor elimination on connected supports. Step 16 chooses
Lemma 4.8's integer schedule and runs a budget-feasible finite member when
Lemma 4.7 accepts it. Step 17 builds finite Section 5.1–5.2 objects: sparse
F2 systems, integer lifts, Lemma 5.1's coefficient division, and the truncated
product. Those checks are exhaustive on a budgeted cube. The dense matrix
remains limited to side 256.

The target is reference Theorem 1.1's exact count with its time/space bounds.
The general-formula reduction and complete asymptotic pipeline remain open.
Finite comparisons and numerical admission do not establish those bounds.
Step 17 does not implement Theorem 1.1.
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
from .common_count import CommonCountPlan, CommonCountResult, count_placed_formula
from .conditional_balance import ConditionalSchedules, choose_schedules
from .reference_parameters import (
    AsymptoticChoice, ChoiceBudgets, Lemma48Proposal, ParameterTuple,
    ReferenceParameters, check_reference_parameters, choose_asymptotic_parameters,
    choose_asymptotic_parameters_for_cnf, lemma_4_8_integers, lemma_4_8_proposal,
)
from .structured_auxiliary import StructuredAuxiliary, build_structured_auxiliary
from .sparse_placement import (
    SparsePlacedFormula, SparsePlacementBudgets, SparseScheduleOutcome,
    check_sparse_design, place_sparse_formula, schedule_sparse_placement,
)
from .component_count import ComponentCountResult, count_cnf_components
from .support_count import SupportCountPlan, SupportCountResult, count_placed_formula_by_support
from .tensor_count import ContractionBudgets, ContractionPlan, ContractionResult, count_placed_formula_by_contraction
from .sparse_systems import (
    STEP17_SCOPE, BudgetRefusal, CommonZeroResult, ParityPolynomials, Section5Budgets,
    SparseSystem, TruncatedProduct, brute_common_zeros, build_truncated_product,
    construct_parity_polynomials, count_common_zeros, integer_lift, monomial_allowance,
    sparse_polynomial, sparse_system,
)
from .dimacs import parse_dimacs, read_dimacs, format_dimacs
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
from .scheduler import (
    BalanceCheck,
    ScheduleOutcome,
    Step7Example,
    check_balance,
    schedule_placement,
    step7_examples,
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
    "SparsePlacedFormula",
    "SparsePlacementBudgets",
    "SparseScheduleOutcome",
    "ComponentCountResult",
    "SupportCountPlan",
    "SupportCountResult",
    "ContractionBudgets",
    "ContractionPlan",
    "ContractionResult",
    "count_placed_formula_by_contraction",
    "STEP17_SCOPE",
    "BudgetRefusal",
    "CommonZeroResult",
    "ParityPolynomials",
    "Section5Budgets",
    "SparseSystem",
    "TruncatedProduct",
    "brute_common_zeros",
    "build_truncated_product",
    "construct_parity_polynomials",
    "count_common_zeros",
    "integer_lift",
    "monomial_allowance",
    "sparse_polynomial",
    "sparse_system",
    "count_placed_formula_by_support",
    "check_sparse_design",
    "place_sparse_formula",
    "schedule_sparse_placement",
    "count_cnf_components",
    "parse_dimacs",
    "read_dimacs",
    "format_dimacs",
    "CommonCountPlan",
    "CommonCountResult",
    "ConditionalSchedules",
    "AsymptoticChoice",
    "ChoiceBudgets",
    "Lemma48Proposal",
    "ParameterTuple",
    "ReferenceParameters",
    "StructuredAuxiliary",
    "build_structured_auxiliary",
    "check_reference_parameters",
    "choose_asymptotic_parameters",
    "choose_asymptotic_parameters_for_cnf",
    "lemma_4_8_integers",
    "lemma_4_8_proposal",
    "choose_schedules",
    "count_placed_formula",
    "AuxiliarySlice",
    "BalanceCheck",
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
    "ScheduleOutcome",
    "Step7Example",
    "Step5Case",
    "Step5Report",
    "StructuredReplacement",
    "TIMING_DISCLAIMER",
    "analyze_step4_case",
    "analyze_step5_case",
    "blocked_sum",
    "build_auxiliary_slice",
    "check_balance",
    "eval_H_hasse_jet",
    "hand_built_examples",
    "hand_built_replacement",
    "hermite_eval",
    "make_design",
    "place_formula",
    "popcount_design",
    "run_parameter_scan",
    "schedule_placement",
    "popcount_group_instance",
    "popcount_growth",
    "scan_specs",
    "sparsify",
    "step3_examples",
    "step4_cases",
    "step5_cases",
    "step7_examples",
]

__version__ = "0.18.0"
