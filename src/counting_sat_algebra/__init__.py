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
product. Those checks are exhaustive on a budgeted cube. Step 18 builds the
Section 5.3–5.5 replacement factors and sums them with the Section 3 routine
in A_ell on a checked finite design. Lemma 4.8 does not form D on these
branches. The dense matrix remains limited to side 256. Step 19
checks finite Section 6 objects: violation polynomials, a budgeted Lemma 6.1
matrix, greedy row selection, and disjoint sparsification that keeps the
declared n. Proposition 6.4's integers are checked. The rho bound is not
claimed.
Step 20 adds Appendix D's finite coefficient rings, CNF tensor networks,
direct modular leaves, Cauchy local coordinates and tapered Hermite recovery.
Step 21 constructs finite collision labels and packet layouts, and checks
the prospective child topologies. Step 22 checks finite weighted modes.
Step 23 supplies Lean-compatible syntax, exact dyadic parameter search and
a streamed disjoint counting frontier on the original declared n. Step 24
constructs ternary Newton digits and exact accepted/marked query counts
through shared-marker truncation and finite squarefree polynomial sums.
Step 25 constructs the scheduled sparse ternary code with regenerated rows
and selects history rows by accumulated marked-violation mass before the
separate final accepted query. Step 26 computes the early high-incidence
threshold and the high-variable frontend, and routes the hash leaf queries
through the streamed conditioned rank loop. Step 27 constructs the counting
extension ring, the fixed nodes and the Q*J > degree Hermite placement.
Step 28 streams the correction pairs and re-reads the per-node Taylor
coefficients by Newton Hermite. The DIMACS Lean backend connects the native
and hash leaves to the same frontier.

The active roadmap follows the supplied SqrtLogSpace counting program.
Packing, native placement and Python time/space refinement remain open.
The supplied Lean library has theorem proof source; local kernel validation
is pending. Finite comparisons do not establish general asymptotic bounds.
Appendix D packet tensors and recursive contraction remain a later route.
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
from .algebra_count import (
    STEP18_SCOPE, AlgebraBudgets, AlgebraDesign, IncidenceCut, ReplacementSumResult,
    count_replacement_sum, design_problems, finite_admitted_design, grouped_unit_design,
    h_degree_bound, incidence_cut, input_degree_bound, sum_supplied_design,
)
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
from .lean_code import (
    CodeBudgets, CodeMeter, CodeRefusal, CodeSearchResult, SparseCodeCertificate,
    SparseCodeParameters, SparseCodeResult, SparseCodeRows, SparseTannerCode,
    canonical_hash_candidates, certify_sparse_code, construct_sparse_code, radix_words,
    search_local_ternary_matrix, search_walsh_shifts, sparse_code_parameters,
)
from .lean_encoding import (
    cnf_numbers, count_word, decode_cnf, decode_count_word, encode_cnf, formula_bit_length,
    literal_code,
)
from .lean_greedy import (
    STEP25_SCOPE, GreedyBudgets, GreedyCountResult, count_cnf_greedy, occurrence_bound,
    select_history_rows,
)
from .lean_ext import (
    STEP27_SCOPE, CountingExtension, ExtensionBudgets, ExtensionRefusal,
    counting_extension, extension_schedule, hermite_reconstruct, is_irreducible_f3,
    layout_cutoff, layout_energy, layout_q, layout_r, node_code_words, node_element,
    node_stream, search_irreducible_polynomial, taylor_jets,
)
from .lean_high import (
    STEP26_SCOPE, HighBudgets, HighFrontend, HighLoopResult, count_cnf_high_loop,
    early_log, early_threshold, high_frontend,
)
from .lean_taylor import (
    STEP28_SCOPE, TaylorBudgets, TaylorRefusal, decode_pair, decode_stream,
    difference_coefficient, direct_sum, encode_pair, iter_pairs,
    node_coefficient_stream, pair_count, reconstruct_from_stream,
    streamed_boolean_sum, streamed_taylor_count,
)
from .lean_parameters import SparseParameter, choose_sparse_parameter, parameter_test
from .lean_ring import (
    DenseElement, LeanAlgebraBudgets, LeanAlgebraRefusal, TernaryIndicator, TruncatedRing,
    construct_ternary_indicator,
)
from .lean_query import (
    STEP24_SCOPE, TernaryHashRow, TernaryQueryResult, count_cnf_ternary, count_ternary_query,
)
from .lean_sparsify import (
    STEP23_SCOPE, PriorityChoice, StreamingBudgets, StreamingCountResult, count_cnf_streaming,
    normalize_lean_cnf, select_priority, split_priority,
)
from .support_count import SupportCountPlan, SupportCountResult, count_placed_formula_by_support
from .tensor_count import ContractionBudgets, ContractionPlan, ContractionResult, count_placed_formula_by_contraction
from .modular_ring import ModularBudgetRefusal, ModularRing, extension_ring_for_size
from .modular_network import (
    STEP20_SCOPE, BooleanTensorNetwork, CnfNetworkEncoding, ModularBudgets,
    ModularCountResult, NetworkContractionResult, TensorVertex,
    contract_network, count_cnf_modular, encode_cnf_network,
)
from .modular_recovery import (
    DeformationResult, ModularRecoveryResult, PacketLayout, TaperedJet, TaperedNode, TaperedPlan,
    cauchy_local_polynomial, embed_network, evaluate_deformation,
    plan_tapered_recovery, recover_network_modular, recover_tapered_jets,
)
from .modular_packing import (
    STEP21_SCOPE, AppendixDIntegers, AppendixDProposal, ChildTopology, CollisionChoice,
    CollisionLabels, PackedNetwork, PackingBudgets, PackingCheck, PackingHypotheses,
    PackingPlan, PackingResult, appendix_d_integers, appendix_d_proposal,
    check_packing, child_topology, choose_collision_labels, pack_network, packing_hypotheses,
)
from .modular_modes import (
    STEP22_SCOPE, FrozenModeDescription, ModeBudgets, ModeRequest, ModeWitness,
    SupportRecord, WeightedMode, admit_exponent_tuple, admit_support, build_mode_table,
    recover_crossing_residue, recover_network_with_modes, retained_degree_ok,
    support_admissible,
)
from .equation_selection import (
    STEP19_SCOPE, DisjointCountResult, GreedySelectionResult, InitialMatrixResult,
    Proposition64Schedule, Section6Budgets, SparsificationTree,
    build_disjoint_sparsification, clause_violation_polynomial, construct_initial_matrix,
    count_via_disjoint_sparsification, expand_binary_row, greedy_select_equations,
    lemma_6_1_parameters, normalize_cnf, proposition_6_4_schedule,
)
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
    "STEP23_SCOPE", "PriorityChoice", "SparseParameter", "StreamingBudgets",
    "StreamingCountResult", "choose_sparse_parameter", "cnf_numbers", "count_cnf_streaming",
    "count_word", "decode_cnf", "decode_count_word", "encode_cnf", "formula_bit_length",
    "literal_code", "normalize_lean_cnf", "parameter_test", "select_priority", "split_priority",
    "DenseElement", "LeanAlgebraBudgets", "LeanAlgebraRefusal", "TernaryIndicator", "TruncatedRing",
    "construct_ternary_indicator", "STEP24_SCOPE", "TernaryHashRow", "TernaryQueryResult",
    "count_cnf_ternary", "count_ternary_query", "STEP25_SCOPE", "GreedyBudgets",
    "GreedyCountResult", "count_cnf_greedy", "occurrence_bound", "select_history_rows",
    "CodeBudgets", "CodeMeter", "CodeRefusal", "CodeSearchResult", "SparseCodeCertificate",
    "SparseCodeParameters", "SparseCodeResult", "SparseCodeRows", "SparseTannerCode",
    "canonical_hash_candidates", "certify_sparse_code", "construct_sparse_code", "radix_words",
    "search_local_ternary_matrix", "search_walsh_shifts", "sparse_code_parameters",
    "STEP26_SCOPE", "HighBudgets", "HighFrontend", "HighLoopResult", "count_cnf_high_loop",
    "early_log", "early_threshold", "high_frontend",
    "STEP27_SCOPE", "CountingExtension", "ExtensionBudgets", "ExtensionRefusal",
    "counting_extension", "extension_schedule", "hermite_reconstruct",
    "is_irreducible_f3", "layout_cutoff", "layout_energy", "layout_q", "layout_r",
    "node_code_words", "node_element", "node_stream", "search_irreducible_polynomial",
    "taylor_jets",
    "STEP28_SCOPE", "TaylorBudgets", "TaylorRefusal", "decode_pair", "decode_stream",
    "difference_coefficient", "direct_sum", "encode_pair", "iter_pairs",
    "node_coefficient_stream", "pair_count", "reconstruct_from_stream",
    "streamed_boolean_sum", "streamed_taylor_count",
    "STEP22_SCOPE",
    "FrozenModeDescription",
    "ModeBudgets",
    "ModeRequest",
    "ModeWitness",
    "SupportRecord",
    "WeightedMode",
    "admit_exponent_tuple",
    "admit_support",
    "build_mode_table",
    "recover_crossing_residue",
    "recover_network_with_modes",
    "retained_degree_ok",
    "support_admissible",
    "STEP21_SCOPE",
    "AppendixDIntegers",
    "AppendixDProposal",
    "ChildTopology",
    "CollisionChoice",
    "CollisionLabels",
    "PackedNetwork",
    "PackingBudgets",
    "PackingCheck",
    "PackingHypotheses",
    "PackingPlan",
    "PackingResult",
    "appendix_d_integers",
    "appendix_d_proposal",
    "check_packing",
    "child_topology",
    "choose_collision_labels",
    "pack_network",
    "packing_hypotheses",
    "STEP20_SCOPE",
    "ModularBudgetRefusal",
    "ModularRing",
    "extension_ring_for_size",
    "BooleanTensorNetwork",
    "CnfNetworkEncoding",
    "ModularBudgets",
    "ModularCountResult",
    "NetworkContractionResult",
    "TensorVertex",
    "contract_network",
    "count_cnf_modular",
    "encode_cnf_network",
    "DeformationResult",
    "ModularRecoveryResult",
    "PacketLayout",
    "TaperedJet",
    "TaperedNode",
    "TaperedPlan",
    "cauchy_local_polynomial",
    "embed_network",
    "evaluate_deformation",
    "plan_tapered_recovery",
    "recover_network_modular",
    "recover_tapered_jets",
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
    "STEP19_SCOPE",
    "DisjointCountResult",
    "GreedySelectionResult",
    "InitialMatrixResult",
    "Proposition64Schedule",
    "Section6Budgets",
    "SparsificationTree",
    "build_disjoint_sparsification",
    "clause_violation_polynomial",
    "construct_initial_matrix",
    "count_via_disjoint_sparsification",
    "expand_binary_row",
    "greedy_select_equations",
    "lemma_6_1_parameters",
    "normalize_cnf",
    "proposition_6_4_schedule",
    "STEP18_SCOPE",
    "AlgebraBudgets",
    "AlgebraDesign",
    "IncidenceCut",
    "ReplacementSumResult",
    "count_replacement_sum",
    "design_problems",
    "finite_admitted_design",
    "grouped_unit_design",
    "h_degree_bound",
    "incidence_cut",
    "input_degree_bound",
    "sum_supplied_design",
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

__version__ = "0.28.0"
