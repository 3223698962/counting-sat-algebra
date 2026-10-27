"""Executable alignment checks for reference Sections 3–5.

Section 5 checks are the finite monomial, Lemma 5.1, and truncated-product
objects. They do not implement Section 6 or Theorem 1.1.
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
    Section5Budgets, beigel_tarui_complement, brute_common_zeros, build_truncated_product,
    construct_parity_polynomials, count_common_zeros, evaluate_truncated_product, integer_lift,
    monomial_allowance, sparse_polynomial, sparse_system, truncated_indicator_residue,
)


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
               'finite members of the Lemma 4.8 parameter schedule, and finite Section 5.1–5.2 '
               'sparse-system checks. Section 6, Appendix D, and Theorem 1.1 are not implemented; '
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
