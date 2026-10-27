"""Executable alignment checks for reference Sections 3 and 4."""

from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction

from .cnf import CNF, brute_count
from .component_count import count_cnf_components
from .common_count import count_placed_formula
from .hermite import eval_H_hasse_jet
from .replacement_matrix import build_auxiliary_slice, step3_examples
from .reference_parameters import check_reference_parameters
from .scheduler import schedule_placement
from .structured_auxiliary import build_structured_auxiliary
from .sparse_placement import schedule_sparse_placement


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
                     boolean_vector_bound_bits=counted.plan.boolean_vector_bound.bit_length(),
                     boolean_vectors=counted.boolean_vectors, reason=counted.reason))
    print(f"  [N] global common sum: {'REFUSED before enumeration' if refused else 'FAIL'}", flush=True)
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
