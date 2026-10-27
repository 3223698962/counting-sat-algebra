"""Executable alignment checks for reference Sections 3 and 4."""

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
from .reference_parameters import check_reference_parameters
from .scheduler import schedule_placement
from .structured_auxiliary import build_structured_auxiliary
from .sparse_placement import place_sparse_formula, schedule_sparse_placement
from .support_count import _support_components, count_placed_formula_by_support
from .tensor_count import ContractionBudgets, _contract_component, _plan_components, count_placed_formula_by_contraction


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


REFERENCE_SUITES = (
    (9, run_reference_design_suite),
    (10, run_common_count_suite),
    (11, run_structured_auxiliary_suite),
    (12, run_sparse_placement_suite),
    (13, run_component_count_suite),
    (14, run_support_count_suite),
    (15, run_tensor_count_suite),
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
               'component counts, disjoint replacement supports and finite connected Taylor contraction; '
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
