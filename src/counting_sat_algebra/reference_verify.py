"""Executable finite alignment checks for reference Sections 3–6 and Appendix D.

Section 5 checks are the finite monomial, Lemma 5.1, and truncated-product
objects, plus the finite Section 5.3–5.5 replacement sum in A_ell. Section 6
checks are budgeted equation selection and disjoint sparsification on a fixed
n. Step 20 adds Appendix D coefficient rings, networks and local recovery.
Step 21 adds finite collision packing and prospective child topologies.
Recursive modular contraction and Theorem 1.1 remain open; the checks do not
claim Proposition 6.4's running time.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from fractions import Fraction
from pathlib import Path

from .cnf import CNF, brute_count
from .component_count import count_cnf_components
from .common_count import count_placed_formula
from .hermite import eval_H_hasse_jet
from .replacement_matrix import build_auxiliary_slice, place_formula, step3_examples
from .reference_parameters import (
    ChoiceBudgets, check_reference_parameters, choose_asymptotic_parameters_for_cnf,
    lemma_4_8_integers,
)
from .scheduler import schedule_placement
from .structured_auxiliary import build_structured_auxiliary
from .sparse_placement import place_sparse_formula, schedule_sparse_placement
from .support_count import _support_components, count_placed_formula_by_support
from .tensor_count import ContractionBudgets, _contract_component, _plan_components, count_placed_formula_by_contraction
from .sparse_systems import (
    Section5Budgets, SparseSystem, beigel_tarui_complement, brute_common_zeros, build_truncated_product,
    construct_parity_polynomials, count_common_zeros, evaluate_truncated_product, integer_lift,
    monomial_allowance, sparse_polynomial, sparse_system, truncated_indicator_residue,
)
from .algebra_count import (
    AlgebraBudgets, count_replacement_sum, h_degree_bound, incidence_cut, input_degree_bound,
)
from .equation_selection import (
    STEP19_SCOPE, Section6Budgets, construct_initial_matrix, count_via_disjoint_sparsification,
    count_violations_common_zeros, greedy_select_equations, kernel_support_check, largest_s,
    lemma_6_1_parameters, proposition_6_4_schedule, choose_branching_M,
)
from .modular_ring import ModularRing, extension_ring_for_size, evaluate_polynomial, interpolate_polynomial
from .modular_network import (
    STEP20_SCOPE, BooleanTensorNetwork, ModularBudgets, TensorVertex,
    contract_network, count_cnf_modular,
)
from .modular_recovery import (
    PacketLayout, cauchy_local_polynomial, embed_network, evaluate_deformation,
    recover_network_modular,
)
from .modular_packing import (
    STEP21_SCOPE, PackingBudgets, appendix_d_integers, appendix_d_proposal,
    child_topology, choose_collision_labels, pack_network, packing_hypotheses,
)
from .modular_ring import ModularBudgetRefusal


def _row(name, n, m, ok, label, **details):
    return dict(name=name, n=n, m=m, brute=details.pop('brute', -1), ok=ok, label=label if ok else 'FAIL', **details)


def run_reference_design_suite() -> list[dict]:
    print('=== (K) Reference conditional expectation and rational design ===', flush=True)
    cnf = CNF(1024, [[(v, bool(v % 2))] for v in range(1024)])
    outcome = schedule_placement(cnf, r=2, Q=2, R=1)
    ok = (outcome.accepted and outcome.design_method == 'linear' and outcome.balance.admitted
          and outcome.multiplicity_ok and outcome.group_size_ok and outcome.final_potential <= outcome.initial_potential <= 2)
    rows = [_row('admit-n1024-linear', cnf.n, cnf.m, bool(ok), 'CHECKED',
                 multiplicity=outcome.multiplicity, multiplicity_bound=str(outcome.multiplicity_bound),
                 max_group_size=outcome.max_group_size, trials=outcome.trials,
                 initial_potential=str(outcome.initial_potential), final_potential=str(outcome.final_potential),
                 count_checked=False)]
    print(f"  [K] n=1024: placement {'CHECKED' if ok else 'FAIL'}; multiplicity={outcome.multiplicity}, largest group={outcome.max_group_size}", flush=True)
    for n, Q, R, admitted in ((1024, 2, 1, False), (16384, 3, 2, True)):
        parameters = check_reference_parameters(n, tuple((v,) for v in range(n)), r=2, Q=Q, R=R, tau=Fraction(1, 2))
        correct = parameters.admitted == admitted
        rows.append(_row(f'lemma4.7-n{n}', n, n, correct, 'CHECKED',
                         numerical_hypotheses_admitted=parameters.admitted,
                         tail=str(parameters.tail_probability), tail_limit=str(parameters.tail_limit),
                         placement_checked=False, count_checked=False, reasons=list(parameters.reasons)))
        print(f"  [K] n={n}, Q={Q}, R={R}: tail={parameters.tail_probability}, limit={parameters.tail_limit}; hypotheses {'admitted' if parameters.admitted else 'rejected'}", flush=True)
    return rows


def run_common_count_suite() -> list[dict]:
    print('=== (L) Section 3 shared-measurement counting and Taylor recovery ===', flush=True)
    rows = []
    for example in step3_examples():
        placed = example.placed
        result = count_placed_formula(placed)
        expected = brute_count(placed.cnf)
        jets_match = result.counted and all(
            jet == eval_H_hasse_jet(placed.cnf, placed.family, Fraction(t), result.plan.q)
            for t, jet in enumerate(result.jets, 1)
        )
        ok = result.count == expected and jets_match
        rows.append(_row(example.name, placed.cnf.n, placed.cnf.m, bool(ok), 'COUNTED',
                         brute=expected, count=result.count, jets=[list(map(str, jet)) for jet in result.jets],
                         plan=asdict(result.plan), boolean_vectors=result.boolean_vectors,
                         measurement_tuples=result.measurement_tuples, factor_evaluations=result.factor_evaluations))
        print(f"  [L] {example.name}: #F={result.count}, pairs={result.plan.substitution_pairs}, local vectors={result.boolean_vectors}/{result.plan.boolean_vector_bound}; {'COUNTED' if ok else 'FAIL'}", flush=True)
    return rows


def run_structured_auxiliary_suite() -> list[dict]:
    print('=== (M) Structured auxiliary operator, including R=32, Q=192 ===', flush=True)
    rows = []
    for r, Q, R in ((2, 2, 2), (2, 3, 2)):
        dense = build_auxiliary_slice(r, Q, R)
        structured = build_structured_auxiliary(r, Q, R)
        ok = all(structured.matrix_entry(row, col) == dense.M0[u][y]
                 and structured.replacement_entry(row, col) == dense.rows[u][y]
                 for u, row in enumerate(dense.labels) for y, col in enumerate(dense.labels))
        rows.append(_row(f'structured-d{dense.dimension}', 0, 0, ok, 'CHECKED', dimension=dense.dimension, entry_comparisons=dense.dimension ** 2))
        print(f"  [M] D={dense.dimension}: all matrix and replacement entries {'CHECKED' if ok else 'FAIL'}", flush=True)
    full = build_structured_auxiliary(2, 192, 32)
    row = ((0,) * 192, (1,) * 192)
    ok = (full.dimension == 2 ** 384 and full.replacement_entry(row, row) == (Fraction(1),)
          and full.replacement_entry(row, (row[1], row[0])) == (Fraction(0),)
          and all(full.f_values[j] == (-1) ** j for j in full.interval))
    rows.append(_row('full-parameter-coordinate-queries', 0, 0, ok, 'CHECKED',
                     dimension=str(full.dimension), measurement_modes=str(full.measurement_mode_count),
                     dense_matrix_allocated=False, modes_enumerated=False, count_checked=False))
    print(f"  [M] R=32, Q=192, D=2**384: coordinate queries {'CHECKED' if ok else 'FAIL'}", flush=True)
    return rows


def run_sparse_placement_suite() -> list[dict]:
    print('=== (N) Sparse placement satisfying all Lemma 4.7 hypotheses ===', flush=True)
    cnf = CNF(16384, [[(v, bool(v % 2))] for v in range(16384)])
    outcome = schedule_sparse_placement(cnf, r=2, Q=3, R=2)
    ok = outcome.accepted and outcome.placed.statistics.polynomial_entries == 59392
    rows = [_row('admit-n16384-sparse', cnf.n, cnf.m, ok, 'CHECKED',
                 numerical_hypotheses_admitted=outcome.parameters.admitted,
                 tail=str(outcome.parameters.tail_probability),
                 storage=asdict(outcome.placed.statistics) if outcome.placed else None,
                 multiplicity=outcome.placed.copies if outcome.placed else None,
                 group_ranks=outcome.design.group_ranks if outcome.design else (),
                 count_checked=False)]
    print(f"  [N] n=16384: hypotheses and sparse design {'CHECKED' if ok else 'FAIL'}", flush=True)
    if not ok:
        print(outcome.reason, flush=True)
        return rows
    storage = outcome.placed.statistics
    print(f"  [N] polynomial entries={storage.polynomial_entries}; dense slots={storage.dense_polynomial_slots}; measurement entries={storage.measurement_entries}", flush=True)
    counted = count_placed_formula(outcome.placed)
    refused = counted.status == 'refused' and counted.count is None and counted.boolean_vectors == 0
    rows.append(_row('n16384-common-budget', cnf.n, cnf.m, refused, 'REFUSED',
                     count=None, q=counted.plan.q, J=counted.plan.J,
                     substitution_pair_bound_bits=counted.plan.substitution_pairs.bit_length(),
                     measurement_tuple_bound_bits=[value.bit_length() for value in counted.plan.measurement_tuples],
                     local_cube_bound_bits=[[value.bit_length() for value in sizes] for sizes in counted.plan.local_cube_sizes],
                     boolean_vector_bound_bits=counted.plan.boolean_vector_bound.bit_length(),
                     boolean_vectors=counted.boolean_vectors, reason=counted.reason))
    print(f"  [N] global common sum: {'REFUSED before enumeration' if refused else 'FAIL'}", flush=True)
    return rows


def run_support_count_suite() -> list[dict]:
    print('=== (P) Global placed polynomial by disjoint replacement supports ===', flush=True)
    cnf = CNF(16384, [[(v, bool(v % 2))] for v in range(16384)])
    outcome = schedule_sparse_placement(cnf, r=2, Q=3, R=2)
    if not outcome.accepted:
        return [_row('n16384-support-contraction', cnf.n, cnf.m, False, 'COUNTED', reason=outcome.reason)]
    result = count_placed_formula_by_support(outcome.placed, vector_budget=49152)
    ok = (result.counted and result.count == 1 and result.plan.components == 10240
          and result.plan.largest_component == 4 and result.plan.max_component_degree == 4
          and result.plan.boolean_vector_bound == 49152 and result.boolean_vectors == 96
          and result.components_counted == result.plan.components
          and result.local_evaluations + result.memo_hits == result.components_counted)
    rows = [_row('n16384-support-contraction', cnf.n, cnf.m, ok, 'COUNTED',
                 count=result.count, expected_count=1,
                 count_validation='one assignment for independent alternating-polarity unit clauses',
                 global_placement_checked=True, numerical_hypotheses_admitted=outcome.parameters.admitted,
                 same_replacement_polynomials=True, memo_enabled=True,
                 plan=asdict(result.plan) if result.plan else None,
                 components_counted=result.components_counted, local_evaluations=result.local_evaluations,
                 memo_hits=result.memo_hits, boolean_vectors=result.boolean_vectors,
                 factor_evaluations=result.factor_evaluations)]
    print(f"  [P] n=16384: #F={result.count}; support size <= 4; vectors={result.boolean_vectors}/49152; {'COUNTED' if ok else 'FAIL'}", flush=True)

    # Reuse the exact same placement, but evaluate every support component.
    # This checks the advertised bound without relying on repeated local factors.
    uncached = count_placed_formula_by_support(outcome.placed, vector_budget=49152, memo_entries=0)
    uncached_ok = (uncached.counted and uncached.count == result.count == 1
                   and uncached.plan == result.plan and uncached.boolean_vectors == 49152
                   and uncached.components_counted == uncached.local_evaluations == 10240
                   and uncached.memo_hits == 0)
    rows.append(_row('n16384-support-no-memo', cnf.n, cnf.m, uncached_ok, 'COUNTED',
                     count=uncached.count, expected_count=1,
                     count_validation='one assignment for independent alternating-polarity unit clauses',
                     global_placement_checked=True, numerical_hypotheses_admitted=outcome.parameters.admitted,
                     same_replacement_polynomials=True, memo_enabled=False,
                     count_matches_memoized=uncached.counted and result.counted and uncached.count == result.count,
                     plan=asdict(uncached.plan) if uncached.plan else None,
                     components_counted=uncached.components_counted, local_evaluations=uncached.local_evaluations,
                     memo_hits=uncached.memo_hits, boolean_vectors=uncached.boolean_vectors,
                     factor_evaluations=uncached.factor_evaluations))
    print(f"  [P] same placement, no memo: #F={uncached.count}; vectors={uncached.boolean_vectors}/49152; {'COUNTED' if uncached_ok else 'FAIL'}", flush=True)

    cnf = CNF(13, [[(v, False), (v + 1, False)] for v in range(12)])
    placed = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 1),
                                  ((0, 0),) * 13, ((0, 0),) * 12)
    result = count_placed_formula_by_support(placed)
    refused = result.status == 'refused' and result.count is None and result.boolean_vectors == 0
    rows.append(_row('connected-support-budget', cnf.n, cnf.m, refused, 'REFUSED',
                     count=None, plan=asdict(result.plan), boolean_vectors=result.boolean_vectors,
                     reason=result.reason))
    print(f"  [P] connected support n=13: {'REFUSED before enumeration' if refused else 'FAIL'}", flush=True)
    return rows


def run_component_count_suite() -> list[dict]:
    print('=== (O) Exact component counts on supplied CNFs ===', flush=True)
    cnf = CNF(8, [[(0, False), (1, False), (2, True)], [(3, False), (4, True), (5, False)]])
    result = count_cnf_components(cnf)
    ok = result.counted and result.count == brute_count(cnf) == 196
    rows = [_row('component-demo-n8', 8, 2, ok, 'COUNTED', brute=196, count=result.count,
                 components=result.components, free_variables=result.free_variables,
                 boolean_vectors=result.boolean_vectors)]
    print(f"  [O] n=8: #F={result.count}, declared free variables={result.free_variables}; {'COUNTED' if ok else 'FAIL'}", flush=True)
    cnf = CNF(32769, [[(3*i, False), (3*i+1, True), (3*i+2, False)] for i in range(10923)])
    parameters = check_reference_parameters(cnf.n, tuple(tuple(v for v, _neg in clause) for clause in cnf.clauses),
                                             r=2, Q=3, R=2, tau=Fraction(1, 2))
    result = count_cnf_components(cnf)
    ok = result.counted and result.count == 7 ** 10923 and parameters.admitted
    rows.append(_row('disjoint-3cnf-n32769', cnf.n, cnf.m, ok, 'COUNTED',
                     count_hex=hex(result.count) if result.counted else None,
                     expected_count='7**10923', count_validation='closed form for disjoint three-literal clauses',
                     numerical_hypotheses_admitted=parameters.admitted,
                     global_placement_checked=False, components=result.components,
                     local_evaluations=result.local_evaluations, memo_hits=result.memo_hits,
                     boolean_vectors=result.boolean_vectors))
    print(f"  [O] n=32769: #F=7**10923; local evaluations={result.local_evaluations}; memo hits={result.memo_hits}; {'COUNTED' if ok else 'FAIL'}", flush=True)
    return rows


def run_tensor_count_suite() -> list[dict]:
    print('=== (Q) Connected replacement supports by finite Taylor contraction ===', flush=True)

    def details(result):
        return dict(plan=asdict(result.plan) if result.plan else None,
                    components_counted=result.components_counted, boolean_vectors=result.boolean_vectors,
                    factor_evaluations=result.factor_evaluations, jet_products=result.jet_products,
                    peak_table_entries=result.peak_table_entries,
                    peak_table_coefficients=result.peak_table_coefficients)

    def path_placement(n, free=0):
        cnf = CNF(n + free, [[(v, False), (v + 1, False)] for v in range(n - 1)])
        placed = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 1),
                                      ((0, 0),) * cnf.n, ((0, 0),) * cnf.m)
        return cnf, placed

    cnf, placed = path_placement(13)
    legacy = count_placed_formula_by_support(placed)
    result = count_placed_formula_by_contraction(placed)
    expected = brute_count(cnf)
    ok = (result.counted and result.count == expected == 610 and legacy.status == 'refused'
          and result.plan.max_bucket_variables == 2 and result.boolean_vectors == 98)
    rows = [_row('connected-support-n13', cnf.n, cnf.m, ok, 'COUNTED',
                 brute=expected, count=result.count, previous_support_status=legacy.status,
                 same_replacement_polynomials=True, **details(result))]
    print(f"  [Q] connected n=13: #F={result.count}, local vectors={result.boolean_vectors}; {'COUNTED' if ok else 'FAIL'}", flush=True)

    cnf, placed = path_placement(1024, free=2)
    a, b = 0, 1
    for _ in range(1026):
        a, b = b, a + b
    result = count_placed_formula_by_contraction(placed)
    supplied = count_cnf_components(cnf, backend='tensor', max_component_variables=1024)
    ok = (result.counted and result.count == 4 * a and result.plan.largest_component == 1024
          and result.plan.free_variables == 2 and result.plan.max_bucket_variables == 2
          and result.boolean_vectors == 8186 and result.peak_table_entries == 4094
          and supplied.counted and supplied.count == result.count and supplied.boolean_vectors == 8186)
    rows.append(_row('connected-path-n1024-free2', cnf.n, cnf.m, ok, 'COUNTED',
                     count_hex=hex(result.count) if result.counted else None,
                     expected_count='4*Fibonacci(1026)',
                     count_validation='independent recurrence for bit strings with no adjacent zeros',
                     component_backend_count_matches=supplied.counted and supplied.count == result.count,
                     component_backend_boolean_vectors=supplied.boolean_vectors,
                     same_replacement_polynomials=True, **details(result)))
    print(f"  [Q] connected n=1024 plus two free bits: #F=4*Fibonacci(1026), local vectors={result.boolean_vectors}; {'COUNTED' if ok else 'FAIL'}", flush=True)

    cnf = CNF(2, [[(0, False)], [(1, True)]])
    dense = place_formula(cnf, build_auxiliary_slice(2, 2, 2),
                          ((0, 0), (1, 0)), ((1, 1), (0, 1)))
    placed = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 2),
                                  dense.variable_schedules, dense.clause_schedules)
    result = count_placed_formula_by_contraction(placed)
    components, free = _support_components(placed)
    plan, local_plans = _plan_components(components, free, ContractionBudgets())
    stats = dict(boolean_vectors=0, factor_evaluations=0, jet_products=0,
                 peak_table_entries=0, peak_table_coefficients=0)
    jet = _contract_component(components[0], local_plans[0], stats)
    expected_jet = eval_H_hasse_jet(cnf, dense.family, Fraction(1), plan.max_component_degree)
    ok = result.counted and result.count == 1 and len(components) == 1 and any(jet[1:]) and jet == expected_jet
    rows.append(_row('coupled-support-full-taylor', cnf.n, cnf.m, bool(ok), 'COUNTED',
                     brute=1, count=result.count, jet=list(map(str, jet)),
                     independent_jet=list(map(str, expected_jet)),
                     coefficient_validation_vectors=stats['boolean_vectors'], **details(result)))
    print(f"  [Q] coupled replacements: nonconstant Taylor coefficients {'CHECKED' if ok else 'FAIL'}", flush=True)

    cnf = CNF(13, [[(v, False) for v in range(13)]])
    placed = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 1),
                                  ((0, 0),) * cnf.n, ((0, 0),))
    result = count_placed_formula_by_contraction(placed)
    refused = result.status == 'refused' and result.count is None and result.boolean_vectors == 0
    rows.append(_row('wide-tensor-budget', cnf.n, cnf.m, refused, 'REFUSED',
                     count=None, boolean_vectors=result.boolean_vectors, reason=result.reason))
    print(f"  [Q] width-13 factor: {'REFUSED before numerical tables' if refused else 'FAIL'}", flush=True)
    return rows


def _alternating_units(n: int) -> CNF:
    return CNF(n, [[(v, bool(v % 2))] for v in range(n)])


def run_asymptotic_family_suite() -> list[dict]:
    """Choose Lemma 4.8's schedule and run a budget-feasible finite member.

    The closed form is checked exactly. On these n it fails before D is formed.
    The executed tuple is whatever ``check_reference_parameters`` accepts inside
    the sparse-placement budgets. Counting uses the existing backends and returns
    no partial count when a budget fails. This is not an asymptotic runtime proof.
    """
    print('=== (R) Asymptotic parameter family: Lemma 4.8 and a finite admitted member ===', flush=True)
    rows = []
    large_n = 16384
    nu, r, Q, R = lemma_4_8_integers(1 << 19, 1)
    guard_fails = r >= 2 and Q <= R
    rows.append(_row(
        'lemma48-integers-lambda-2**19', 0, 0, guard_fails, 'CHECKED',
        log_n=1 << 19, log_log_n=nu, r=r, Q=Q, R=R,
        preliminary_guard_passed=False, dimension_formed=False, matrix_allocated=False,
        numerical_hypotheses_admitted=False,
        note='Q <= R, so tau and D are not formed; this lambda is not an executed instance',
    ))
    print(f'  [R] lambda=2**19: r={r}, Q={Q}, R={R}; preliminary guard {"fails before D" if guard_fails else "FAIL"}', flush=True)

    cnf = _alternating_units(large_n)
    choice = choose_asymptotic_parameters_for_cnf(cnf)
    proposal = choice.proposal
    proposal_ok = (
        proposal.log_n == 14 and proposal.log_log_n == 4 and proposal.beta == 1
        and proposal.r == 0 and proposal.Q == 0 and proposal.R == 16
        and proposal.dimension_formed is False and proposal.dimension is None
        and proposal.tau is None and proposal.admitted is False and proposal.lemma_4_7_admitted is False
        and any('r=0 < 2' in reason for reason in proposal.reasons)
    )
    rows.append(_row(
        'lemma48-units-n16384', cnf.n, cnf.m, proposal_ok, 'CHECKED',
        log_n=proposal.log_n, log_log_n=proposal.log_log_n, beta=proposal.beta,
        r=proposal.r, Q=proposal.Q, R=proposal.R, tau=None, dimension=None,
        dimension_formed=proposal.dimension_formed, matrix_allocated=False,
        numerical_hypotheses_admitted=proposal.admitted, reasons=list(proposal.reasons),
        placement_started=False,
    ))
    print(f'  [R] n=16384 Lemma 4.8: r={proposal.r}, Q={proposal.Q}, R={proposal.R}; D not formed; hypotheses rejected', flush=True)

    smaller = _alternating_units(8192)
    smaller_choice = choose_asymptotic_parameters_for_cnf(smaller)
    smaller_ok = (
        smaller_choice.source == 'none' and smaller_choice.selected is None
        and smaller_choice.admitted is False
        and any('binomial tail 7/16' in item for item in smaller_choice.obstacles)
        and any('2(b+w)D^2 = 16384' in item for item in smaller_choice.obstacles)
    )
    rows.append(_row(
        'finite-member-units-n8192', smaller.n, smaller.m, smaller_ok, 'CHECKED',
        source=smaller_choice.source, numerical_hypotheses_admitted=False,
        examined=smaller_choice.examined, admitted_candidates=smaller_choice.admitted_candidates,
        obstacles=list(smaller_choice.obstacles), reasons=list(smaller_choice.reasons),
        placement_started=False, matrix_allocated=False,
    ))
    print(f'  [R] n=8192: no Lemma 4.7 tuple inside the budgets; {"CHECKED failure" if smaller_ok else "FAIL"}', flush=True)

    selected = choice.selected
    choice_ok = (
        proposal_ok and choice.source == 'finite-budget-member' and selected is not None
        and selected.placeable and choice.admitted and selected.reference.admitted
        and (selected.r, selected.Q, selected.R, selected.tau, selected.A0, selected.dimension) == (2, 3, 2, Fraction(5, 16), 0, 64)
        and selected.tail_probability == Fraction(5, 32) and selected.tail_limit == Fraction(5, 32)
        and selected.meets_r_at_least_32 is False and selected.meets_six_r_at_most_q is False
        and 2 * selected.tau * cnf.n == 10240
    )
    rows.append(_row(
        'finite-member-units-n16384', cnf.n, cnf.m, bool(choice_ok), 'CHECKED',
        source=choice.source, numerical_hypotheses_admitted=bool(selected and selected.reference.admitted),
        r=None if selected is None else selected.r, Q=None if selected is None else selected.Q,
        R=None if selected is None else selected.R,
        tau=None if selected is None else str(selected.tau),
        hand_tau='1/2', measurement_allowance=None if selected is None else int(2 * selected.tau * cnf.n),
        hand_measurement_allowance=16384, interval=None if selected is None else list(selected.interval),
        dimension=None if selected is None else selected.dimension,
        tail=None if selected is None else str(selected.tail_probability),
        tail_limit=None if selected is None else str(selected.tail_limit),
        alphabet=None if selected is None else selected.alphabet,
        trials=None if selected is None else selected.trials,
        examined=choice.examined, admitted_candidates=choice.admitted_candidates,
        meets_r_at_least_32=None if selected is None else selected.meets_r_at_least_32,
        meets_six_r_at_most_q=None if selected is None else selected.meets_six_r_at_most_q,
        count_checked=False,
        note='Least tau meeting the exact tail. Not the Lemma 4.8 closed form, and not R>=32, 6R<=Q.',
    ))
    print(f'  [R] n=16384 finite member: r=2, Q=3, R=2, tau=5/16; {"CHECKED" if choice_ok else "FAIL"}', flush=True)

    tight = choose_asymptotic_parameters_for_cnf(cnf, budgets=ChoiceBudgets(max_trials=1))
    tight_ok = (
        tight.admitted and tight.selected is not None and tight.placeable is False
        and tight.selected.r == 2 and tight.selected.trials > 1
    )
    rows.append(_row(
        'finite-member-trials-budget', cnf.n, cnf.m, tight_ok, 'REFUSED',
        source=tight.source, numerical_hypotheses_admitted=tight.admitted, placeable=tight.placeable,
        trials=None if tight.selected is None else tight.selected.trials, trial_budget=1,
        placement_started=False, matrix_allocated=False, reasons=list(tight.reasons),
    ))
    print(f'  [R] trial budget 1: {"REFUSED before schedule search" if tight_ok else "FAIL"}', flush=True)

    try:
        build_auxiliary_slice(2, 192, 32)
        dense_refused = False
    except ValueError as exc:
        dense_refused = '2**384' in str(exc) and 'not built' in str(exc)
    rows.append(_row(
        'dense-full-m0-still-refused', 0, 0, dense_refused, 'REFUSED',
        r=2, Q=192, R=32, matrix_allocated=False,
        note='The dense side-length cap is unchanged. Step 16 does not allocate full M0.',
    ))
    print(f'  [R] dense M0 at R=32, Q=192: {"REFUSED before allocation" if dense_refused else "FAIL"}', flush=True)

    if not choice_ok or selected is None:
        rows.append(_row('chosen-placement-n16384', cnf.n, cnf.m, False, 'CHECKED', reason='chooser did not admit a placeable tuple'))
        return rows

    outcome = schedule_sparse_placement(cnf, r=selected.r, Q=selected.Q, R=selected.R, tau=selected.tau, A0=selected.A0)
    storage = outcome.placed.statistics if outcome.placed is not None else None
    rank_sum = sum(outcome.design.group_ranks[0]) if outcome.design and outcome.design.group_ranks else None
    placed_ok = (
        outcome.accepted and outcome.parameters.admitted and outcome.parameters.tau == Fraction(5, 16)
        and storage is not None and storage.polynomial_entries == 59392 and storage.measurement_entries == 14336
        and storage.dense_polynomial_slots == 536870912 and outcome.placed.copies == 256
        and outcome.design.max_group_size == 8192 and outcome.design.group_ranks == ((2048, 2048),) * 3
        and rank_sum is not None and rank_sum <= 2 * outcome.parameters.tau * cnf.n
    )
    rows.append(_row(
        'chosen-placement-n16384', cnf.n, cnf.m, bool(placed_ok), 'CHECKED',
        numerical_hypotheses_admitted=outcome.parameters.admitted, tau=str(outcome.parameters.tau),
        storage=asdict(storage) if storage else None,
        multiplicity=outcome.placed.copies if outcome.placed else None,
        group_ranks=[list(ranks) for ranks in outcome.design.group_ranks] if outcome.design else [],
        measurement_allowance=int(2 * outcome.parameters.tau * cnf.n),
        independent_measurements=rank_sum, count_checked=False,
        reason=outcome.reason,
    ))
    print(f'  [R] n=16384 sparse placement at the chosen tau: {"CHECKED" if placed_ok else "FAIL"}', flush=True)
    if not placed_ok or outcome.placed is None:
        return rows

    counted = count_placed_formula(outcome.placed)
    common_refused = counted.status == 'refused' and counted.count is None and counted.boolean_vectors == 0
    rows.append(_row(
        'chosen-n16384-common-budget', cnf.n, cnf.m, common_refused, 'REFUSED',
        count=None, q=counted.plan.q, J=counted.plan.J,
        boolean_vector_bound_bits=counted.plan.boolean_vector_bound.bit_length(),
        boolean_vectors=counted.boolean_vectors, reason=counted.reason,
    ))
    print(f'  [R] common backend: {"REFUSED before enumeration" if common_refused else "FAIL"}', flush=True)

    support = count_placed_formula_by_support(outcome.placed, vector_budget=49152)
    support_ok = (
        support.counted and support.count == 1 and support.plan is not None
        and support.plan.components == 10240 and support.plan.largest_component == 4
        and support.plan.max_component_degree == 4 and support.plan.boolean_vector_bound == 49152
        and support.boolean_vectors == 96 and support.local_evaluations == 13 and support.memo_hits == 10227
    )
    rows.append(_row(
        'chosen-n16384-support', cnf.n, cnf.m, support_ok, 'COUNTED',
        count=support.count, expected_count=1,
        count_validation='one assignment for independent alternating-polarity unit clauses',
        same_chosen_placement=True, memo_enabled=True,
        plan=asdict(support.plan) if support.plan else None,
        boolean_vectors=support.boolean_vectors, local_evaluations=support.local_evaluations,
        memo_hits=support.memo_hits,
        note='Disjoint replacement supports are not a general connected-instance bound.',
    ))
    print(f'  [R] support backend: #F={support.count}; vectors={support.boolean_vectors}; {"COUNTED" if support_ok else "FAIL"}', flush=True)

    tensor = count_placed_formula_by_contraction(outcome.placed, budgets=ContractionBudgets())
    if tensor.counted:
        tensor_ok = (
            tensor.count == 1 and tensor.plan is not None and tensor.plan.components == 10240
            and tensor.plan.largest_component == 4 and tensor.plan.max_bucket_variables == 4
            and tensor.boolean_vectors == 176128 and tensor.factor_evaluations == 98304
            and tensor.jet_products == 262144 and tensor.peak_table_entries == 72
            and tensor.plan.planning_steps == 268288 and tensor.plan.peak_table_coefficients == 360
        )
        rows.append(_row(
            'chosen-n16384-tensor', cnf.n, cnf.m, bool(tensor_ok), 'COUNTED',
            count=tensor.count, expected_count=1,
            count_validation='one assignment for independent alternating-polarity unit clauses',
            same_chosen_placement=True, memo_enabled=False,
            plan=asdict(tensor.plan) if tensor.plan else None,
            boolean_vectors=tensor.boolean_vectors, factor_evaluations=tensor.factor_evaluations,
            jet_products=tensor.jet_products, peak_table_entries=tensor.peak_table_entries,
            note='Finite Taylor elimination on these supports. Not Appendix D.',
        ))
        print(f'  [R] tensor backend: #F={tensor.count}; vectors={tensor.boolean_vectors}; {"COUNTED" if tensor_ok else "FAIL"}', flush=True)
    else:
        tensor_refused = tensor.status == 'refused' and tensor.count is None and tensor.boolean_vectors == 0
        rows.append(_row(
            'chosen-n16384-tensor', cnf.n, cnf.m, tensor_refused, 'REFUSED',
            count=None, boolean_vectors=tensor.boolean_vectors, reason=tensor.reason,
            plan=asdict(tensor.plan) if tensor.plan else None,
            note='Refused before numerical tables. Not a partial count.',
        ))
        print(f'  [R] tensor backend: {"REFUSED before numerical tables" if tensor_refused else "FAIL"}', flush=True)
    return rows


def _height(y: int) -> int:
    return y * (y - 1)


def _paper_ell2_truncated(y1: int, y2: int) -> int:
    """The three-term product displayed after Lemma 5.2."""
    return (
        (1 - y1) * (1 - y2)
        + (2 * y1 - 1) * _height(y1) * (1 - y2)
        + (1 - y1) * (2 * y2 - 1) * _height(y2)
    )


def _paper_ell2_omitted(y1: int, y2: int) -> int:
    return (2 * y1 - 1) * (2 * y2 - 1) * _height(y1) * _height(y2)


def run_section5_suite() -> list[dict]:
    """Finite Section 5.1–5.2 checks. Not Theorem 2.2's running time."""
    print('=== (S) Section 5 sparse systems, Lemma 5.1, truncated product ===', flush=True)
    rows = []
    max_bits = 0
    identity_ok = True
    samples = range(-6, 13)
    for ell in range(1, 5):
        parity = construct_parity_polynomials(ell)
        again = construct_parity_polynomials(ell)
        raw = beigel_tarui_complement(ell)
        max_bits = max(max_bits, parity.coefficient_bits)
        identity_ok = (
            identity_ok
            and parity.coefficients == again.coefficients
            and raw[0] == 1
            and all(coefficient == 0 for coefficient in raw[1:ell])
            and all(parity.indicator_agrees(y) for y in samples)
        )
    rows.append(_row(
        'lemma-5.1-integers', 0, 0, identity_ok, 'CHECKED',
        ell_values=[1, 2, 3, 4], integer_samples=len(samples),
        max_coefficient_bits=max_bits, count_checked=False,
        note='Even y is congruent to 1 and odd y to 0 modulo 2^ell. Finite check only.',
    ))
    print(f'  [S] Lemma 5.1 on ell=1..4: {"CHECKED" if identity_ok else "FAIL"}', flush=True)

    parity = construct_parity_polynomials(2)
    raw = beigel_tarui_complement(2)
    product = build_truncated_product(2, 2)
    ell2_ok = (
        parity.coefficients == ((1, -1), (-1, 2))
        and raw == (1, 0, -3, 2)
        and [list(term.powers) for term in product.terms] == [[0, 0], [0, 1], [1, 0]]
        and product.shape_ok()
        and all(term.linear_factors <= 2 and term.quadratic_factors <= 1 for term in product.terms)
        and (1, 1) not in {term.powers for term in product.terms}
    )
    for y1 in range(-4, 7):
        for y2 in range(-4, 7):
            got = evaluate_truncated_product(product, parity, (y1, y2))
            full = parity.evaluate(y1) * parity.evaluate(y2)
            if (
                got != _paper_ell2_truncated(y1, y2)
                or got + _paper_ell2_omitted(y1, y2) != full
                or _paper_ell2_omitted(y1, y2) % 4
                or truncated_indicator_residue(parity, (y1, y2)) != got % 4
            ):
                ell2_ok = False
    rows.append(_row(
        'lemma-5.1-ell-2-product', 0, 2, ell2_ok, 'CHECKED',
        coefficients=[[1, -1], [-1, 2]], p_coefficients=[1, 0, -3, 2],
        retained_powers=[[0, 0], [0, 1], [1, 0]], omitted_powers=[1, 1],
        linear_factor_bound=2, quadratic_factor_bound=1, count_checked=False,
        note='The omitted h(y1)h(y2) term is divisible by 4 at integer inputs.',
    ))
    print(f'  [S] ell=2 displayed product: {"CHECKED" if ell2_ok else "FAIL"}', flush=True)

    cancelled = sparse_polynomial(4, [[0], [0], [1]], degree_bound=2, monomial_bound=3)
    constant = sparse_polynomial(2, [[]], degree_bound=1, monomial_bound=1)
    zero = sparse_polynomial(4, [], degree_bound=2, monomial_bound=2)
    collapsed = sparse_polynomial(4, [[0, 0, 1]], degree_bound=2, monomial_bound=1)
    lift_system = sparse_system(4, [[[0, 1], [2]]], k=2, c=1)
    pointwise_lift = True
    for mask in range(16):
        assignment = [(mask >> variable) & 1 for variable in range(4)]
        lift = integer_lift(lift_system.equations[0], assignment)
        x0, x1, x2 = assignment[:3]
        pointwise_lift = pointwise_lift and (lift % 2 == 0) == ((x0 * x1 + x2) % 2 == 0)
    rejected_long = False
    try:
        sparse_system(4, [[[0], [1], [2]]], k=1, c=1)
    except ValueError:
        rejected_long = True
    rejected_degree = False
    try:
        sparse_system(4, [[[0, 1, 2]]], k=2, c=1)
    except ValueError:
        rejected_degree = True
    normalization_ok = (
        cancelled.monomials == (frozenset({1}),)
        and constant.is_one and integer_lift(constant, [0, 1]) == 1
        and zero.is_zero and integer_lift(zero, [1, 1, 1, 1]) == 0
        and collapsed.monomials == (frozenset({0, 1}),) and collapsed.degree == 2
        and monomial_allowance(1, 1) == 0
        and monomial_allowance(3, 1) == 1 and monomial_allowance(3, 2) == 2
        and monomial_allowance(4, 1) == 2 and monomial_allowance(1024, 3) == 1000
        and lift_system.equations[0].degree == 2 and pointwise_lift
        and rejected_long and rejected_degree
        and lift_system.meets_theorem_size is False
    )
    rows.append(_row(
        'sparse-f2-lift', 4, 1, normalization_ok, 'CHECKED',
        monomial_allowance_n4_c1=2, monomial_allowance_n1024_c3=1000,
        duplicates_cancelled=True, empty_monomial_is_one=True, count_checked=False,
        note='(log2 n)^c is an exact integer comparison. n=4 is below the theorem size.',
    ))
    print(f'  [S] sparse normalization and lift parity: {"CHECKED" if normalization_ok else "FAIL"}', flush=True)

    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)
    counted = count_common_zeros(system)
    brute = brute_common_zeros(system)
    counted_parity = construct_parity_polynomials(5)
    counted_product = build_truncated_product(2, 5)
    pointwise = counted_product.shape_ok()
    for mask in range(16):
        assignment = [(mask >> variable) & 1 for variable in range(4)]
        lifts = [integer_lift(equation, assignment) for equation in system.equations]
        residue = evaluate_truncated_product(counted_product, counted_parity, lifts) % 32
        indicator = 0 if any(lift % 2 for lift in lifts) else 1
        pointwise = pointwise and residue == indicator
        pointwise = pointwise and truncated_indicator_residue(counted_parity, lifts) == indicator
    counted_ok = (
        counted.counted and counted.count == brute == 4
        and counted.boolean_vectors == 16 and counted.terms == 15 and counted.ell == 5
        and counted.coefficient_bits == 7 and counted.input_degree_bound == 20
        and counted.linear_factor_bound == 2 and counted.quadratic_factor_bound == 4
        and counted.term_table_allocated and counted.coefficients_constructed and counted.enumerated
        and pointwise and counted.count is not None
    )
    rows.append(_row(
        'tiny-common-zeros-n4', system.n, system.m_test, counted_ok, 'COUNTED',
        brute=brute, count=counted.count, ell=counted.ell, terms=counted.terms,
        boolean_vectors=counted.boolean_vectors, coefficient_bits=counted.coefficient_bits,
        input_degree_bound=counted.input_degree_bound,
        linear_factor_bound=counted.linear_factor_bound,
        quadratic_factor_bound=counted.quadratic_factor_bound,
        pointwise_indicators=pointwise, term_table_allocated=True,
        note='Exhaustive residue on 16 points. Not Section 5.3, Section 6, or Theorem 1.1.',
    ))
    print(f'  [S] n=4 common zeros: N(P)={counted.count}; {"COUNTED" if counted_ok else "FAIL"}', flush=True)

    empty = count_common_zeros(sparse_system(5, [], k=2, c=1))
    empty_ok = (
        empty.counted and empty.count == brute_common_zeros(sparse_system(5, [], k=2, c=1)) == 32
        and empty.boolean_vectors == 0 and empty.term_table_allocated is False
        and empty.coefficients_constructed is False and empty.enumerated is False
    )
    rows.append(_row(
        'empty-system-n5', 5, 0, empty_ok, 'COUNTED',
        brute=32, count=empty.count, boolean_vectors=0, term_table_allocated=False,
        note='An empty equation list is 2^n. No cube was scanned.',
    ))
    print(f'  [S] empty system: N(P)={empty.count}; {"COUNTED" if empty_ok else "FAIL"}', flush=True)

    constant_system = sparse_system(4, [[[]]], k=2, c=1)
    constant = count_common_zeros(constant_system)
    constant_ok = (
        constant_system.equations[0].is_one
        and constant.counted and constant.count == 0 and constant.boolean_vectors == 0
        and constant.term_table_allocated is False and constant.enumerated is False
        and brute_common_zeros(constant_system) == 0
    )
    rows.append(_row(
        'constant-one-equation', 4, 1, constant_ok, 'COUNTED',
        brute=0, count=constant.count, boolean_vectors=0, term_table_allocated=False,
        note='The constant equation 1 = 0 has no solutions. No partial search.',
    ))
    print(f'  [S] constant 1 = 0: N(P)={constant.count}; {"COUNTED" if constant_ok else "FAIL"}', flush=True)

    term_refusal = count_common_zeros(system, Section5Budgets(max_terms=2))
    term_ok = (
        term_refusal.refused and term_refusal.count is None and term_refusal.boolean_vectors == 0
        and term_refusal.term_table_allocated is False and term_refusal.coefficients_constructed is False
        and 'term table not allocated' in term_refusal.reason
    )
    rows.append(_row(
        'term-table-budget', system.n, system.m_test, term_ok, 'REFUSED',
        count=None, boolean_vectors=0, term_table_allocated=False, coefficients_constructed=False,
        term_budget=2, reason=term_refusal.reason,
    ))
    print(f'  [S] term budget: {"REFUSED before allocation" if term_ok else "FAIL"}', flush=True)

    bit_refusal = count_common_zeros(system, Section5Budgets(max_coefficient_bits=1))
    bit_ok = (
        bit_refusal.refused and bit_refusal.count is None and bit_refusal.boolean_vectors == 0
        and bit_refusal.term_table_allocated is False and bit_refusal.coefficients_constructed is False
        and 'not constructed' in bit_refusal.reason
    )
    rows.append(_row(
        'coefficient-bit-budget', system.n, system.m_test, bit_ok, 'REFUSED',
        count=None, boolean_vectors=0, term_table_allocated=False, coefficients_constructed=False,
        coefficient_bit_budget=1, reason=bit_refusal.reason,
    ))
    print(f'  [S] coefficient bits: {"REFUSED before the term table" if bit_ok else "FAIL"}', flush=True)

    cube_refusal = count_common_zeros(system, Section5Budgets(max_boolean_vectors=1))
    cube_ok = (
        cube_refusal.refused and cube_refusal.count is None and cube_refusal.boolean_vectors == 0
        and cube_refusal.term_table_allocated is False and cube_refusal.coefficients_constructed is False
        and cube_refusal.enumerated is False and 'no assignments scanned' in cube_refusal.reason
    )
    rows.append(_row(
        'cube-budget', system.n, system.m_test, cube_ok, 'REFUSED',
        count=None, boolean_vectors=0, term_table_allocated=False, coefficients_constructed=False,
        boolean_vector_budget=1, reason=cube_refusal.reason,
    ))
    print(f'  [S] cube budget: {"REFUSED before the scan" if cube_ok else "FAIL"}', flush=True)
    return rows


def _boolean_coefficient_sum(system, ell: int) -> tuple[int, ...]:
    """Cube sum of the truncated factors in Z[V]/(V^ell). Not the Section 3 routine."""
    parity = construct_parity_polynomials(ell)
    active = [equation for equation in system.equations if not equation.is_zero]
    total = [0] * ell
    for mask in range(1 << system.n):
        bits = [(mask >> variable) & 1 for variable in range(system.n)]
        accumulator = [0] * ell
        accumulator[0] = 1
        for equation in active:
            lift = integer_lift(equation, bits)
            height = lift * (lift - 1)
            column = []
            height_power = 1
            for constant, slope in parity.coefficients:
                column.append((constant + slope * lift) * height_power)
                height_power *= height
            combined = [0] * ell
            for left, left_value in enumerate(accumulator):
                if left_value == 0:
                    continue
                for right, right_value in enumerate(column):
                    if left + right < ell and right_value:
                        combined[left + right] += left_value * right_value
            accumulator = combined
        for slot, value in enumerate(accumulator):
            total[slot] += value
    return tuple(total)


def run_section5_replacement_suite() -> list[dict]:
    """Finite Sections 5.3–5.5. Not Section 6, Appendix D, or Theorem 2.2."""
    print('=== (T) Section 5.3–5.5 replacement sum in A_ell ===', flush=True)
    rows = []
    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)
    cut = incidence_cut(system)
    degree = input_degree_bound(system.k, 2, 5)
    replaced = h_degree_bound(system.k, 2, 5, 1)
    bounds_ok = (
        degree == 20 and replaced == 20
        and cut.incidence_total == 8 and cut.occurrence_cap == 4 and cut.heavy == ()
        and cut.neutral_scopes == ((0,), (1,), (2,), (3,))
        and cut.equation_scopes == ((0, 1, 2), (3,))
        and cut.log_n == 2 and cut.sqrt_log_n == 2
    )
    rows.append(_row(
        'replacement-degree-and-incidence', system.n, system.m_test, bounds_ok, 'CHECKED',
        input_degree_bound=degree, h_degree_bound=replaced, incidence_total=cut.incidence_total,
        occurrence_cap=cut.occurrence_cap, heavy=[], neutral_factors=4,
        section_6_implemented=False, appendix_d_implemented=False, count_checked=False,
        note='Equation (5.5) and the Section 5.5 occurrence cap. No sum is run in this row.',
    ))
    print(f'  [T] degree d={degree}, d_H={replaced}, b0={cut.occurrence_cap}: {"CHECKED" if bounds_ok else "FAIL"}', flush=True)

    counted = count_replacement_sum(system)
    oracle = _boolean_coefficient_sum(system, 5)
    brute = brute_common_zeros(system)
    counted_ok = (
        counted.counted and counted.count == brute == 4
        and counted.h0_coefficients == oracle
        and counted.unreduced_sum == sum(oracle) == 2308
        and counted.unreduced_sum % 32 == 4
        and counted.boolean_vectors == 399 and counted.vector_bound == 714
        and counted.q == 0 and counted.jet_length == 1 and counted.design_Q == 21
        and counted.replacement_degree == 1 and counted.slope == '1/21'
        and counted.input_degree_bound == 20 and counted.h_degree_bound == 20
        and counted.neutral_factors == 4 and counted.equations_used == 2
        and counted.branches == 1 and counted.branch_counts == (4,)
        and counted.sample_at_one != '4' and counted.enumerated and counted.jets_allocated
        and counted.lemma_48_guard_passed is False and counted.lemma_48_dimension_formed is False
        and counted.lemma_48_r == 1 and counted.lemma_48_Q == 0 and counted.lemma_48_R == 10
    )
    rows.append(_row(
        'replacement-sum-n4', system.n, system.m_test, counted_ok, 'COUNTED',
        brute=brute, count=counted.count, ell=counted.ell,
        h0_coefficients=list(counted.h0_coefficients), unreduced_sum=counted.unreduced_sum,
        boolean_vectors=counted.boolean_vectors, vector_bound=counted.vector_bound,
        factor_evaluations=counted.factor_evaluations, q=counted.q, jet_length=counted.jet_length,
        design_Q=counted.design_Q, replacement_degree=counted.replacement_degree, slope=counted.slope,
        input_degree_bound=counted.input_degree_bound, h_degree_bound=counted.h_degree_bound,
        sample_at_one=counted.sample_at_one, neutral_factors=counted.neutral_factors,
        note='S_ell(H(0)) residue 4 matches the cube oracle and the F2 scan. The row is not Lemma 4.8.',
    ))
    print(
        f'  [T] n=4 replacement sum: N(P)={counted.count}; vectors={counted.boolean_vectors}; '
        f'{"COUNTED" if counted_ok else "FAIL"}',
        flush=True,
    )

    branched_system = sparse_system(4, [[[0]], [[0]], [[0]], [[0]], [[0]], [[0], [1]]], k=1, c=1)
    branched = count_replacement_sum(branched_system)
    branched_brute = brute_common_zeros(branched_system)
    branched_ok = (
        branched.counted and branched.count == branched_brute == 4
        and branched.heavy_variables == (0,) and branched.occurrence_cap == 6
        and branched.branches == 2 and branched.branch_counts == (4, 0)
        and branched.boolean_vectors == 110 and branched.branches_enumerated == 1
        and branched.enumerated and branched.jets_allocated
    )
    rows.append(_row(
        'high-incidence-branches-n4', branched_system.n, branched_system.m_test, branched_ok, 'COUNTED',
        brute=branched_brute, count=branched.count, branch_counts=list(branched.branch_counts),
        heavy=list(branched.heavy_variables), occurrence_cap=branched.occurrence_cap,
        boolean_vectors=branched.boolean_vectors, branches_enumerated=branched.branches_enumerated,
        note='Fixing the heavy variable splits the cube. The branches add to N(P).',
    ))
    print(
        f'  [T] high-incidence branches: N(P)={branched.count}; counts={list(branched.branch_counts)}; '
        f'{"COUNTED" if branched_ok else "FAIL"}',
        flush=True,
    )

    guard_ok = (
        branched.lemma_48_guard_passed is False
        and branched.lemma_48_dimension_formed is False
        and branched.lemma_48_r == 1 and branched.lemma_48_Q == 0 and branched.lemma_48_R == 8
        and branched.counted
    )
    rows.append(_row(
        'lemma-4.8-guard-on-branch', 3, 1, guard_ok, 'CHECKED',
        r=branched.lemma_48_r, Q=branched.lemma_48_Q, R=branched.lemma_48_R,
        guard_passed=False, dimension_formed=False, count_checked=False,
        note='The restricted branch has n\'=3. Tau and D are not formed. The finite design still summed.',
    ))
    print(
        f'  [T] Lemma 4.8 on the branch: r={branched.lemma_48_r}, Q={branched.lemma_48_Q}, '
        f'R={branched.lemma_48_R}; {"CHECKED" if guard_ok else "FAIL"}',
        flush=True,
    )

    vector_refusal = count_replacement_sum(system, AlgebraBudgets(max_boolean_vectors=1))
    vector_ok = (
        vector_refusal.refused and vector_refusal.count is None
        and vector_refusal.boolean_vectors == 0 and vector_refusal.jets_allocated is False
        and vector_refusal.coefficients_constructed is False and vector_refusal.enumerated is False
        and vector_refusal.vector_bound == 714 and 'no enumeration started' in vector_refusal.reason
    )
    rows.append(_row(
        'replacement-vector-budget', system.n, system.m_test, vector_ok, 'REFUSED',
        count=None, boolean_vectors=0, vector_bound=vector_refusal.vector_bound,
        jets_allocated=False, coefficients_constructed=False, boolean_vector_budget=1,
        reason=vector_refusal.reason,
    ))
    print(f'  [T] vector budget: {"REFUSED before enumeration" if vector_ok else "FAIL"}', flush=True)

    branch_refusal = count_replacement_sum(branched_system, AlgebraBudgets(max_branches=1))
    branch_ok = (
        branch_refusal.refused and branch_refusal.count is None
        and branch_refusal.boolean_vectors == 0 and branch_refusal.branches == 2
        and branch_refusal.branches_enumerated == 0 and branch_refusal.enumerated is False
        and branch_refusal.jets_allocated is False and branch_refusal.coefficients_constructed is False
        and 'no branch enumerated' in branch_refusal.reason
    )
    rows.append(_row(
        'replacement-branch-budget', branched_system.n, branched_system.m_test, branch_ok, 'REFUSED',
        count=None, boolean_vectors=0, branches=2, branches_enumerated=0,
        jets_allocated=False, branch_budget=1, reason=branch_refusal.reason,
    ))
    print(f'  [T] branch budget: {"REFUSED before enumeration" if branch_ok else "FAIL"}', flush=True)

    bit_refusal = count_replacement_sum(system, AlgebraBudgets(max_coefficient_bits=1))
    bit_ok = (
        bit_refusal.refused and bit_refusal.count is None and bit_refusal.boolean_vectors == 0
        and bit_refusal.coefficients_constructed is False and bit_refusal.jets_allocated is False
        and bit_refusal.enumerated is False and 'not constructed' in bit_refusal.reason
    )
    rows.append(_row(
        'replacement-coefficient-budget', system.n, system.m_test, bit_ok, 'REFUSED',
        count=None, boolean_vectors=0, jets_allocated=False, coefficients_constructed=False,
        coefficient_bit_budget=1, reason=bit_refusal.reason,
    ))
    print(f'  [T] coefficient bits: {"REFUSED before the jets" if bit_ok else "FAIL"}', flush=True)
    return rows


def run_section6_suite() -> list[dict]:
    """Finite Section 6. Not Appendix D, not Proposition 6.4's runtime, not Theorem 1.1."""
    print('=== (U) Section 6 equation selection and disjoint sparsification ===', flush=True)
    rows = []
    parameters = lemma_6_1_parameters(3, 3, 1)
    built = construct_initial_matrix(3, 3, 1)
    matrix = built.matrix
    kernel = None if matrix is None else kernel_support_check(
        matrix, numerator=parameters.m, denominator=parameters.P, max_boolean_vectors=64,
    )
    kernel_ok = (
        built.constructed and matrix is not None and kernel is not None and kernel.ok
        and parameters.mu == 2 and parameters.P == 80_000 and parameters.D0 == 2_560_000
        and parameters.q_G == 1 and parameters.h == 0 and parameters.local_search_bits == 0
        and matrix.nrows == 2 and matrix.nrows <= 3 * 3 and matrix.max_row_weight == 3
        and matrix.max_row_weight <= parameters.D0
        and kernel.min_nonzero_weight == 2 and kernel.nonzero_kernel_vectors == 3
        and kernel.vectors_scanned == 8
        and kernel.min_nonzero_weight * parameters.P >= parameters.m
        and built.rows_allocated and built.search_started is False
    )
    rows.append(_row(
        'lemma-6.1-tiny-kernel', 3, 3, kernel_ok, 'CHECKED',
        mu=parameters.mu, P=parameters.P, D0=parameters.D0, rows=matrix.nrows if matrix else 0,
        max_row_weight=matrix.max_row_weight if matrix else None,
        min_nonzero_weight=None if kernel is None else kernel.min_nonzero_weight,
        nonzero_kernel_vectors=None if kernel is None else kernel.nonzero_kernel_vectors,
        vectors_scanned=0 if kernel is None else kernel.vectors_scanned,
        row_bound=9, count_checked=False, appendix_d_implemented=False,
        rho_bound_claimed=False, theorem_1_1_complete=False, lean_proof_complete=False,
        note='Even-parity A(0) on m=3. Every nonzero kernel vector has weight 2 >= 3/80000.',
    ))
    print(
        f'  [U] Lemma 6.1 m=3: min weight {None if kernel is None else kernel.min_nonzero_weight}; '
        f'{"CHECKED" if kernel_ok else "FAIL"}',
        flush=True,
    )

    search = construct_initial_matrix(2048, 2048, 1)
    search_bits = 0 if search.parameters is None else search.parameters.local_search_bits
    search_ok = (
        search.refused and search.matrix is None and search.rows_allocated is False
        and search.search_started is False and search_bits == 524_288
        and 'no local matrix allocated' in search.reason
    )
    rows.append(_row(
        'lemma-6.1-search-budget', 2048, 2048, search_ok, 'REFUSED',
        count=None, local_search_bits=search_bits, rows_allocated=False, search_started=False,
        vectors_scanned=0, appendix_d_implemented=False, rho_bound_claimed=False,
        theorem_1_1_complete=False, lean_proof_complete=False,
        reason=search.reason,
    ))
    print(f'  [U] Lemma 6.1 search: {"REFUSED before allocation" if search_ok else "FAIL"}', flush=True)

    tiny = CNF(2, [[(0, False)], [(1, False)]])
    greedy = greedy_select_equations(tiny, list_length=1, B=1)
    final = None if greedy.equations is None else count_common_zeros(
        SparseSystem(tiny.n, 2, 1, greedy.equations),
    )
    brute = brute_count(tiny)
    greedy_ok = (
        greedy.counted and greedy.count == brute == 1
        and final is not None and final.counted and final.count == 1
        and greedy.rounds_completed == 6 and greedy.counts_performed == 24
        and greedy.candidates_enumerated == 24 and greedy.boolean_vectors_scanned == 96
        and greedy.matrix_allocated and final.boolean_vectors == 4
    )
    rows.append(_row(
        'lemma-6.2-greedy-tiny', tiny.n, tiny.m, greedy_ok, 'COUNTED',
        brute=brute, count=greedy.count, rounds=greedy.rounds_completed,
        list_length=1, counts_performed=greedy.counts_performed,
        boolean_vectors_scanned=greedy.boolean_vectors_scanned,
        common_zero_count=None if final is None else final.count,
        common_zero_vectors=0 if final is None else final.boolean_vectors,
        appendix_d_implemented=False, rho_bound_claimed=False,
        theorem_1_1_complete=False, lean_proof_complete=False,
        note='List length 1 fits. Final N(A) is 1, matching brute #F and count_common_zeros.',
    ))
    print(
        f'  [U] greedy n=2: N(A)={greedy.count}; counts={greedy.counts_performed}; '
        f'{"COUNTED" if greedy_ok else "FAIL"}',
        flush=True,
    )

    candidate = greedy_select_equations(
        tiny, list_length=1, B=1, budgets=Section6Budgets(max_candidate_lists=1),
    )
    candidate_ok = (
        candidate.refused and candidate.count is None and candidate.candidates_enumerated == 0
        and candidate.counts_performed == 0 and candidate.boolean_vectors_scanned == 0
        and candidate.matrix_allocated is False and candidate.equations is None
        and 'no candidates enumerated' in candidate.reason
    )
    rows.append(_row(
        'lemma-6.2-candidate-budget', tiny.n, tiny.m, candidate_ok, 'REFUSED',
        count=None, candidate_budget=1, candidates_enumerated=0, boolean_vectors_scanned=0,
        matrix_allocated=False, appendix_d_implemented=False,
        theorem_1_1_complete=False, lean_proof_complete=False,
        reason=candidate.reason,
    ))
    print(f'  [U] candidate budget: {"REFUSED before enumeration" if candidate_ok else "FAIL"}', flush=True)

    cube = greedy_select_equations(
        tiny, list_length=1, B=1, budgets=Section6Budgets(max_boolean_vectors=1),
    )
    cube_ok = (
        cube.refused and cube.count is None and cube.candidates_enumerated == 0
        and cube.boolean_vectors_scanned == 0 and cube.matrix_allocated is False
        and 'no assignments scanned' in cube.reason
    )
    rows.append(_row(
        'lemma-6.2-count-budget', tiny.n, tiny.m, cube_ok, 'REFUSED',
        count=None, boolean_vector_budget=1, candidates_enumerated=0,
        boolean_vectors_scanned=0, matrix_allocated=False, appendix_d_implemented=False,
        theorem_1_1_complete=False, lean_proof_complete=False,
        reason=cube.reason,
    ))
    print(f'  [U] count budget: {"REFUSED before the cube" if cube_ok else "FAIL"}', flush=True)

    paper_P = lemma_6_1_parameters(tiny.n, tiny.m, 1).P
    paper = greedy_select_equations(tiny, list_length=paper_P, B=1)
    paper_ok = (
        paper_P == 80_000 and paper.refused and paper.count is None
        and paper.candidates_enumerated == 0 and paper.boolean_vectors_scanned == 0
        and paper.matrix_allocated is False and paper.equations is None
    )
    rows.append(_row(
        'lemma-6.2-paper-scale-P', tiny.n, tiny.m, paper_ok, 'REFUSED',
        count=None, P=paper_P, candidates_enumerated=0, boolean_vectors_scanned=0,
        matrix_allocated=False, appendix_d_implemented=False, rho_bound_claimed=False,
        theorem_1_1_complete=False, lean_proof_complete=False,
        reason=paper.reason,
    ))
    print(f'  [U] paper-scale P={paper_P}: {"REFUSED" if paper_ok else "FAIL"}', flush=True)

    sparse = CNF(3, [
        [(0, False), (1, False)],
        [(0, False), (2, False)],
        [(1, False), (2, False)],
    ])
    partitioned = count_via_disjoint_sparsification(
        sparse, k=2, M=2,
        leaf_counter=count_violations_common_zeros,
        verify_disjoint=True,
    )
    sparse_brute = brute_count(sparse)
    sparse_ok = (
        partitioned.counted and partitioned.count == sparse_brute == 4
        and partitioned.leaf_counts == (3, 1) and partitioned.leaves == 2
        and partitioned.nodes == 3 and partitioned.depth == 1
        and partitioned.same_declared_n and partitioned.max_occurrence == 1
        and partitioned.occurrence_bound == 6 and partitioned.max_occurrence <= partitioned.occurrence_bound
        and partitioned.tree is not None
        and all(leaf.n == sparse.n for leaf in partitioned.tree.leaves)
    )
    rows.append(_row(
        'disjoint-sparsification-n3', sparse.n, sparse.m, sparse_ok, 'COUNTED',
        brute=sparse_brute, count=partitioned.count, leaf_counts=list(partitioned.leaf_counts),
        leaves=partitioned.leaves, nodes=partitioned.nodes, depth=partitioned.depth,
        occurrence_bound=partitioned.occurrence_bound, max_occurrence=partitioned.max_occurrence,
        same_declared_n=partitioned.same_declared_n, boolean_vectors_scanned=partitioned.boolean_vectors_scanned,
        step5_copying=False, appendix_d_implemented=False, rho_bound_claimed=False,
        theorem_1_1_complete=False, lean_proof_complete=False,
        note='Leaf common-zero counts sum to #F. Declared n stays 3. Step 5 copying is not this tree.',
    ))
    print(
        f'  [U] disjoint n=3: #F={partitioned.count}; leaves={list(partitioned.leaf_counts)}; '
        f'{"COUNTED" if sparse_ok else "FAIL"}',
        flush=True,
    )

    tree_budget = count_via_disjoint_sparsification(
        sparse, k=2, M=2, budgets=Section6Budgets(max_tree_nodes=1),
    )
    tree_ok = (
        tree_budget.refused and tree_budget.count is None and tree_budget.leaves == 0
        and tree_budget.boolean_vectors_scanned == 0 and tree_budget.tree is None
        and tree_budget.leaf_counts == () and 'no leaves returned' in tree_budget.reason
    )
    rows.append(_row(
        'disjoint-tree-budget', sparse.n, sparse.m, tree_ok, 'REFUSED',
        count=None, node_budget=1, leaves=0, boolean_vectors_scanned=0,
        appendix_d_implemented=False, theorem_1_1_complete=False, lean_proof_complete=False,
        reason=tree_budget.reason,
    ))
    print(f'  [U] tree budget: {"REFUSED before the leaves" if tree_ok else "FAIL"}', flush=True)

    schedule = proposition_6_4_schedule(20)
    branching = choose_branching_M(2, Fraction(1))
    witness = largest_s(521_284, 19)
    schedule_ok = (
        schedule.lambda_n == 5 and schedule.nu_n == 3 and schedule.s == 0
        and schedule.exhaustive_search and schedule.sigma is None
        and schedule.rho_bound_claimed is False
        and schedule.theorem_1_1_complete is False and schedule.lean_proof_complete is False
        and witness == 2 and branching.chosen and branching.M == 64
        and branching.C == 4 and branching.h == 10 and branching.gamma == 4
    )
    rows.append(_row(
        'proposition-6.4-schedule', 20, 0, schedule_ok, 'CHECKED',
        lambda_n=schedule.lambda_n, nu_n=schedule.nu_n, s=schedule.s,
        sigma=schedule.sigma, exhaustive_search=True, witness_lambda=521_284,
        witness_nu=19, witness_s=witness, branching_M_at_sigma_1=branching.M,
        count_checked=False, rho_bound_claimed=False, appendix_d_implemented=False,
        theorem_1_1_complete=False, lean_proof_complete=False,
        note=(
            'n=20 has s=0, so the schedule is the exhaustive-search case. '
            's=2 is only an integer witness at lambda=521284. No rho bound is claimed.'
        ),
    ))
    print(
        f'  [U] schedule n=20: lambda={schedule.lambda_n}, nu={schedule.nu_n}, s={schedule.s}; '
        f'{"CHECKED" if schedule_ok else "FAIL"}',
        flush=True,
    )
    if STEP19_SCOPE and 'Theorem 1.1' not in STEP19_SCOPE:
        rows[-1]['ok'] = False
        rows[-1]['label'] = 'FAIL'
    return rows


def run_appendix_d_suite() -> list[dict]:
    print('=== (V) Appendix D finite modular networks and local recovery ===', flush=True)
    rows = []
    ring = extension_ring_for_size(2, 5, 2)
    unit = ring.element((1, 1, 0, 0, 2))
    V = ring.from_base((0, 1))
    points = tuple(ring.residue_point(j) for j in range(7))
    coefficients = tuple(ring.element((j, j + 1, 0, 0, 2 * j)) for j in range(7))
    values = tuple(evaluate_polynomial(ring, coefficients, point) for point in points)
    ring_ok = (
        ring.extension_degree == 4 and ring.polynomial == 0b10011
        and ring.mul(unit, ring.inverse(unit)) == ring.one and ring.mul(V, V) == ring.zero
        and interpolate_polynomial(ring, points, values) == coefficients
    )
    rows.append(_row(
        'appendix-d-coefficient-ring', 0, 0, ring_ok, 'CHECKED',
        s=2, precision=5, ell=2, extension_degree=4, polynomial_bits=bin(ring.polynomial),
        interpolation_points=7, all_coordinates_checked=True,
    ))

    demo = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])
    counted = count_cnf_modular(demo)
    rows.append(_row(
        'appendix-d-cnf-network', demo.n, demo.m,
        counted.counted and counted.count == brute_count(demo) == 196
        and counted.edges == 12 and counted.vertices == 14 and counted.max_degree == 3
        and counted.free_variables == 2 and counted.boolean_vectors_scanned == 4096
        and counted.tensor_evaluations == 57344,
        'COUNTED', **asdict(counted), clause_ports=6, backend='direct modular leaf',
    ))
    tree = CNF(3, [[(0, False)]] * 4)
    tree_count = count_cnf_modular(tree)
    rows.append(_row(
        'appendix-d-equality-tree', tree.n, tree.m,
        tree_count.counted and tree_count.count == brute_count(tree) == 4
        and tree_count.edges == 10 and tree_count.vertices == 11 and tree_count.max_degree == 3,
        'COUNTED', **asdict(tree_count), clause_ports=4, trivalent_equality_vertices=2,
    ))
    empty, unsat = count_cnf_modular(CNF(5, [])), count_cnf_modular(CNF(5, [[]]))
    rows.append(_row(
        'appendix-d-trivial-formulas', 5, 0,
        empty.counted and empty.count == 32 and unsat.counted and unsat.count == 0
        and empty.edges == unsat.edges == 0 and unsat.network_allocated is False,
        'COUNTED', counts=[empty.count, unsat.count], edges=0, tensor_evaluations=0,
        empty_network_assignments=empty.boolean_vectors_scanned,
    ))

    creation, later = ModularRing(2), ModularRing(8)
    lifted = BooleanTensorNetwork(later, 1, (
        TensorVertex('row', 0, (0,), creation, lambda bits: (-1,)),
        TensorVertex('column', 1, (0,), creation, lambda bits: (6,)),
    ))
    lifted_result = contract_network(lifted)
    rows.append(_row(
        'appendix-d-creation-precision', 0, 1, lifted_result.value == (12,), 'CHECKED',
        creation_precision=2, call_precision=8, value=list(lifted_result.value or ()),
        lifted_entry_values=[3, 2], boolean_vectors_scanned=lifted_result.boolean_vectors_scanned,
        note='Each entry is reduced at precision 2 and lifted separately; the two edge terms sum to 12.',
    ))

    local_ring = ring.at_precision(9)
    kernel = cauchy_local_polynomial(local_ring, 1, 1, 0, 1)
    expected_kernel = tuple(local_ring.scalar(value) for value in (-3, -22, -48, -32))
    identity_ok = all(
        cauchy_local_polynomial(local_ring, 1, 0, x, y)[0] == local_ring.scalar(int(x == y))
        for x in range(2) for y in range(2)
    )
    rows.append(_row(
        'appendix-d-integral-cauchy-coordinates', 0, 0,
        kernel == expected_kernel and identity_ok, 'CHECKED', a=1, L=2, alpha=4,
        t=1, row=0, column=1, exact_scalar_coefficients=[-3, -22, -48, -32],
        precision=9, identity_at_zero_checked=True,
    ))

    base = ModularRing(5, 2)
    network = BooleanTensorNetwork(base, 2, (
        TensorVertex('row', 0, (0, 1), base, lambda bits: (1 + bits[0], 2 - bits[1])),
        TensorVertex('column', 1, (0, 1), base, lambda bits: (1 + bits[1], -bits[0])),
    ))
    extension = embed_network(network, ring.polynomial)
    layout = PacketLayout(1, ((0, 1),))
    recovered = recover_network_modular(extension, layout, 1)
    at_one = evaluate_deformation(extension, layout, 1, ring.zero)
    recovered_ok = (
        recovered.recovered and recovered.value == contract_network(extension).value == ring.from_base((9, 4))
        and at_one.value == ring.from_base((26, 13))
        and recovered.port_pairs_scanned == 224 and recovered.tensor_evaluations == 448
        and recovered.plan is not None and recovered.plan.degree_bound == 6
        and [node.precision for node in recovered.plan.nodes] == [133, 133]
    )
    rows.append(_row(
        'appendix-d-local-recovery', 0, 2, recovered_ok, 'CHECKED',
        coefficients=list(ring.base_coefficients(recovered.value)) if recovered.value is not None else None,
        value_at_one=[26, 13], precision=5, ell=2,
        plan=asdict(recovered.plan) if recovered.plan is not None else None,
        port_pairs_scanned=recovered.port_pairs_scanned, tensor_evaluations=recovered.tensor_evaluations,
        recursive_contraction=False,
    ))

    nonlinear = BooleanTensorNetwork(base, 2, (
        TensorVertex('row', 0, (0, 1), base, lambda bits: (1 + bits[0] * bits[1], 2 - bits[1])),
        network.vertices[1],
    ))
    nonlinear = embed_network(nonlinear, ring.polynomial)
    full_degree = recover_network_modular(nonlinear, layout, 1)
    rows.append(_row(
        'appendix-d-nonlinear-row-recovery', 0, 2,
        full_degree.recovered and full_degree.value == contract_network(nonlinear).value == ring.from_base((8, 5)),
        'CHECKED', coefficients=[8, 5], precision=5, ell=2,
        port_pairs_scanned=full_degree.port_pairs_scanned, tensor_evaluations=full_degree.tensor_evaluations,
        note='The degree-6 rational deformation is checked coefficientwise by the independent pytest oracle.',
    ))

    tapered = recover_network_modular(extension, PacketLayout(2, ((0, None, 1, None),)), 2)
    tapered_ok = (
        tapered.recovered and tapered.value == ring.from_base((9, 4)) and tapered.plan is not None
        and [node.multiplicity for node in tapered.plan.nodes] == [6, 3, 6, 6, 3, 6]
        and [node.precision for node in tapered.plan.nodes] == [149, 101, 149, 149, 101, 149]
        and tapered.port_pairs_scanned == 1440 and tapered.tensor_evaluations == 2880
    )
    rows.append(_row(
        'appendix-d-tapered-precision-levels', 0, 2, tapered_ok, 'CHECKED',
        coefficients=[9, 4], precision=5, ell=2,
        plan=asdict(tapered.plan) if tapered.plan is not None else None,
        port_pairs_scanned=tapered.port_pairs_scanned, tensor_evaluations=tapered.tensor_evaluations,
        recursive_contraction=False,
    ))

    refused = count_cnf_modular(demo, budgets=ModularBudgets(max_boolean_vectors=4095))
    rows.append(_row(
        'appendix-d-network-vector-budget', demo.n, demo.m,
        refused.refused and refused.count is None and not refused.network_allocated
        and refused.boolean_vectors_scanned == refused.tensor_evaluations == 0,
        'REFUSED', **asdict(refused),
    ))
    for name, budgets in (
        ('precision', ModularBudgets(max_precision_bits=132)),
        ('port-pairs', ModularBudgets(max_boolean_vectors=223)),
        ('jet-storage', ModularBudgets(max_jet_slots=751)),
    ):
        refused_recovery = recover_network_modular(extension, layout, 1, budgets=budgets)
        rows.append(_row(
            f'appendix-d-recovery-{name}-budget', 0, 2,
            refused_recovery.refused and refused_recovery.value is None and not refused_recovery.jets_allocated
            and refused_recovery.port_pairs_scanned == refused_recovery.tensor_evaluations == 0,
            'REFUSED', value=None, reason=refused_recovery.reason,
            port_pairs_scanned=0, tensor_evaluations=0, jets_allocated=False,
        ))
    short_field = recover_network_modular(network, layout, 1)
    rows.append(_row(
        'appendix-d-insufficient-residue-points', 0, 2,
        short_field.refused and short_field.value is None and 'too few distinct' in short_field.reason
        and short_field.port_pairs_scanned == 0 and not short_field.jets_allocated,
        'REFUSED', value=None, reason=short_field.reason, field_size=2, required_points=7,
        port_pairs_scanned=0, jets_allocated=False,
    ))
    for row in rows:
        row.update(recursive_contraction=False, theorem_1_1_complete=False, lean_proof_complete=False)
        print(f"  [V] {row['name']}: {row['label']}", flush=True)
    print(f'  {STEP20_SCOPE}', flush=True)
    return rows


def run_packet_packing_suite() -> list[dict]:
    print('=== (W) Appendix D collision labels, packets and child topology ===', flush=True)
    rows = []
    proposal = appendix_d_proposal(16, 20)
    integers = proposal.integers
    rows.append(_row(
        'appendix-d-active-schedule-guard', 0, 20,
        (integers.lambda_s, integers.u, integers.a, proposal.H, proposal.delta) == (4, 0, 0, 8, 13)
        and not proposal.leaf and proposal.hypotheses is None,
        'CHECKED', s=16, lambda_s=4, u=0, a=0, H=8, delta=13,
        L_formed=False, D_arr_formed=False, reason=proposal.reason,
    ))
    witness = appendix_d_integers(131072)
    deeper = appendix_d_integers(262144, 1)
    rows.append(_row(
        'appendix-d-logarithmic-schedule-witness', 0, 0,
        (witness.u, witness.a, witness.previous_a) == (1, 2, 0)
        and (deeper.u, deeper.a, deeper.previous_a) == (1, 6, 4)
        and witness.r == 64 and witness.epsilon == Fraction(1, 131072),
        'CHECKED', lambda_s=131072, depth=0, u=1, a=2, previous_a=0,
        second_lambda=262144, second_depth=1, second_u=1, second_a=6, second_previous_a=4,
        r=64, epsilon='1/131072', size_s_constructed=False,
    ))
    admitted = packing_hypotheses(128, 2, 1, 4)
    too_small = packing_hypotheses(32, 2, 1, 4)
    paper_r = packing_hypotheses(8192, 64, 1, 4)
    rows.append(_row(
        'appendix-d-squared-collision-guard', 0, 0,
        admitted.admitted and not too_small.collision
        and not paper_r.collision and paper_r.density and paper_r.packet_degree,
        'CHECKED', finite_r2_threshold=128, r64_D_arr=4096, r64_collision_threshold=134217728,
        r64_M=8192, r64_D6_admitted=False,
        note='The collision condition contains D_arr squared. These finite tuples are not D.4 schedule members.',
    ))

    def no_entry(bits):
        raise AssertionError('collision packing must not evaluate tensor entries')

    ring = ModularRing(1)

    def matching(M):
        return BooleanTensorNetwork(ring, M, tuple(
            TensorVertex(f'{side}-{edge}', side, (edge,), ring, no_entry)
            for side in (0, 1) for edge in range(M)
        ))

    ends = ((0, 2), (0, 2), (0, 3), (1, 2), (1, 3), (1, 3))
    small = BooleanTensorNetwork(ring, len(ends), tuple(
        TensorVertex(str(v), int(v >= 2), tuple(i for i, edge in enumerate(ends) if v in edge), ring, no_entry)
        for v in range(4)
    ))
    labeling = choose_collision_labels(small, 2, 1, trace=True)

    def pair_oracle(assigned):
        value = Fraction(-len(ends) ** 2, 4)
        for first in ends:
            for second in ends:
                probability = Fraction(1)
                for left, right in zip(first, second):
                    if left == right:
                        factor = Fraction(1)
                    elif left in assigned and right in assigned:
                        factor = Fraction(int(assigned[left] == assigned[right]))
                    else:
                        factor = Fraction(1, 2)
                    probability *= factor
                value += probability
        return value

    assigned = {}
    ce_ok = labeling.initial_potential == Fraction(23, 2) and labeling.final_potential == 1
    for choice in labeling.choices:
        ce_ok &= choice.before == pair_oracle(assigned)
        ce_ok &= all(value == pair_oracle({**assigned, choice.vertex: label})
                     for label, value in enumerate(choice.candidates))
        ce_ok &= choice.after == min(choice.candidates) <= choice.before
        assigned[choice.vertex] = choice.label
    rows.append(_row(
        'appendix-d-shared-parallel-collision-expectation', 0, small.edges, ce_ok, 'CHECKED',
        initial_potential=str(labeling.initial_potential), final_potential=str(labeling.final_potential),
        labels=list(labeling.labels), pair_loads=[list(row) for row in labeling.pair_loads],
        candidate_expectations_checked=True, label_trials=labeling.label_trials,
        score_updates=labeling.score_updates, tensor_evaluations=0,
    ))

    finite = pack_network(matching(128), 2, 1, delta=4)
    finite_ok = (
        finite.packed and finite.check.ok and finite.plan.hypotheses.admitted
        and finite.packing.H_pad == 32 and len(finite.packing.layout.packets) == 64
        and finite.packing.labeling.final_potential == 0 and finite.check.max_child_edges == 128
        and finite.check.max_packet_degree == 4 and not finite.check.child_shrink
    )
    rows.append(_row(
        'appendix-d-finite-D6-packing', 0, 128, finite_ok, 'CHECKED',
        r=2, a=1, delta=4, D6_admitted=True, D4_schedule_member=False,
        H_pad=finite.packing.H_pad if finite.packing else None,
        packets=len(finite.packing.layout.packets) if finite.packing else None,
        max_group_ports=64, max_child_edges=128, max_packet_degree=4, child_shrink=False,
        initial_potential='96', final_potential='0', label_trials=finite.label_trials,
        score_updates=finite.score_updates, tensor_evaluations=finite.tensor_evaluations,
    ))
    multi = pack_network(matching(32), 2, 2, require_hypotheses=False)
    multi_ok = (
        multi.packed and multi.check.ok and multi.packing.H_pad == 4
        and len(multi.packing.layout.packets) == 16 and multi.check.max_packet_degree == 8
        and {sum(low == high for low, high in symbols) for symbols in multi.packing.packet_symbols} == {0, 1, 2}
    )
    rows.append(_row(
        'appendix-d-multicoordinate-packet-orbits', 0, 32, multi_ok, 'CHECKED',
        r=2, a=2, H_pad=4, packets=16, max_packet_degree=8, equal_coordinate_counts=[0, 1, 2],
        D6_admitted=False, tensor_evaluations=multi.tensor_evaluations,
    ))

    large = matching(8192)
    strict = pack_network(large, 64, 1, delta=4)
    rows.append(_row(
        'appendix-d-r64-D6-refusal', 0, large.edges,
        strict.refused and strict.packing is None and not strict.labels_allocated and strict.label_trials == 0,
        'REFUSED', r=64, a=1, delta=4, D6_admitted=False, packing=None,
        labels_allocated=False, packets_allocated=False, label_trials=0, reason=strict.reason,
    ))
    actual = pack_network(large, 64, 1, delta=4, require_hypotheses=False)
    child = child_topology(actual.packing, 0, 0) if actual.packed else None
    actual_ok = (
        actual.packed and actual.check.ok and not actual.plan.hypotheses.admitted
        and actual.packing.H_pad == 2 and len(actual.packing.layout.packets) == 4096
        and actual.packing.labeling.initial_potential == 8190 and actual.packing.labeling.final_potential == 0
        and all(load == 2 for row in actual.packing.labeling.pair_loads for load in row)
        and actual.check.row_group_ports == actual.check.column_group_ports == ((128,) * 64,)
        and child is not None and child.edges == 256 and child.max_degree == 4
        and actual.check.child_shrink and actual.label_trials == 1048576 and actual.score_updates == 1572864
    )
    rows.append(_row(
        'appendix-d-r64-actual-packet-and-child-bounds', 0, large.edges, actual_ok, 'CHECKED',
        r=64, a=1, delta=4, D6_admitted=False, D4_schedule_member=False,
        H_pad=2, packets=4096, initial_potential='8190', final_potential='0', pair_class_load=2,
        max_group_ports=128, max_child_edges=256, max_packet_degree=4, child_shrink=True,
        label_trials=actual.label_trials, score_updates=actual.score_updates,
        tensor_evaluations=actual.tensor_evaluations, child_topology_only=True,
    ))

    cnf = CNF(2, [[(0, False)], [(0, False)]])
    from .modular_network import encode_cnf_network
    encoded = encode_cnf_network(cnf)
    field = extension_ring_for_size(2, encoded.network.ring.precision)
    network = embed_network(encoded.network, field.polynomial)
    automatic = pack_network(network, 2, 1, require_hypotheses=False)
    recovered = recover_network_modular(network, automatic.packing.layout, 1) if automatic.packed else None
    count = (recovered.value[0] << encoded.free_variables) if recovered is not None and recovered.recovered else None
    rows.append(_row(
        'appendix-d-automatic-layout-cnf-recovery', cnf.n, cnf.m,
        count == brute_count(cnf) == 2 and recovered.port_pairs_scanned == 6656
        and recovered.tensor_evaluations == 33280,
        'COUNTED', count=count, brute=2, network_edges=4, free_variables=1,
        port_pairs_scanned=recovered.port_pairs_scanned if recovered else 0,
        tensor_evaluations=recovered.tensor_evaluations if recovered else 0,
        backend='Step 20 local recovery with automatically constructed packets',
        D6_admitted=False,
    ))

    for name, limits, trace in (
        ('label-trials', PackingBudgets(max_label_trials=1), False),
        ('trace', PackingBudgets(max_trace_entries=1), True),
        ('padding', PackingBudgets(max_padded_slots=167), False),
        ('score-updates', PackingBudgets(max_score_updates=1023), False),
    ):
        refused = pack_network(matching(128), 2, 1, delta=4, budgets=limits, trace=trace)
        rows.append(_row(
            f'appendix-d-packing-{name}-budget', 0, 128,
            refused.refused and refused.packing is None and not refused.labels_allocated
            and not refused.packets_allocated and refused.label_trials == refused.score_updates == 0,
            'REFUSED', packing=None, label_trials=0, score_updates=0, tensor_evaluations=0,
            labels_allocated=False, packets_allocated=False, reason=refused.reason,
        ))
    child_refused = False
    child_reason = 'packing did not complete'
    if actual.packed:
        try:
            child_topology(actual.packing, 0, 0, budgets=PackingBudgets(max_child_edges=255))
        except ModularBudgetRefusal as error:
            child_refused = True
            child_reason = str(error)
    rows.append(_row(
        'appendix-d-child-topology-budget', 0, 8192, child_refused, 'REFUSED',
        child_edges=256, child_edge_budget=255, child_topology_allocated=False,
        packing_already_constructed=True, tensor_evaluations=0, reason=child_reason,
    ))
    for row in rows:
        row.update(recursive_contraction=False, theorem_1_1_complete=False, lean_proof_complete=False)
        print(f"  [W] {row['name']}: {row['label']}", flush=True)
    print(f'  {STEP21_SCOPE}', flush=True)
    return rows


REFERENCE_SUITES = (
    (9, run_reference_design_suite),
    (10, run_common_count_suite),
    (11, run_structured_auxiliary_suite),
    (12, run_sparse_placement_suite),
    (13, run_component_count_suite),
    (14, run_support_count_suite),
    (15, run_tensor_count_suite),
    (16, run_asymptotic_family_suite),
    (17, run_section5_suite),
    (18, run_section5_replacement_suite),
    (19, run_section6_suite),
    (20, run_appendix_d_suite),
    (21, run_packet_packing_suite),
)


def write_reference_report(
    path: str | Path, rows: list[dict], steps: list[int], *, failed_step: int | None = None,
) -> None:
    """Save the checks actually run, including failures, with their finite scope."""
    from . import __version__

    report = dict(
        schema_version=1,
        version=__version__,
        goal='Reproduce the local reference exact-counting algorithm and verify its time and space bounds.',
        theorem_1_1_complete=False,
        lean_proof_complete=False,
        scope=('Finite reference Sections 3 and 4 checks, sparse admitted placements, independent '
               'component counts, disjoint replacement supports, finite connected Taylor contraction, '
               'finite members of the Lemma 4.8 parameter schedule, finite Section 5.1–5.2 '
               'sparse-system checks, a finite Section 5.3–5.5 replacement sum in A_ell, and '
               'finite Section 6 equation selection with disjoint sparsification, and finite Appendix D '
               'coefficient rings, CNF networks, Cauchy deformation, tapered Hermite recovery, '
               'collision parameter schedules, packet layouts and prospective child topologies. '
               'Proposition 6.4 does not claim an O_k(rho(n)) bound. '
               'Appendix D recursive contraction and Theorem 1.1 are not implemented; '
               'no general-input asymptotic claim.'),
        wall_times_recorded=False,
        steps=steps,
        failed_step=failed_step,
        all_checks_passed=failed_step is None and bool(rows) and all(row['ok'] for row in rows),
        rows=rows,
    )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
