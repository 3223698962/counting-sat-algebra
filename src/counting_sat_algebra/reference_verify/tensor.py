"""Step 15 finite checks: budgeted Taylor elimination on connected supports. Not Appendix D recursion.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction
from dataclasses import asdict

from .protocol import _row, clean_refusal
from ..cnf import CNF, brute_count
from ..component_count import count_cnf_components
from ..hermite import eval_H_hasse_jet
from ..replacement_matrix import build_auxiliary_slice, place_formula
from ..structured_auxiliary import build_structured_auxiliary
from ..sparse_placement import place_sparse_formula
from ..support_count import _support_components, count_placed_formula_by_support
from ..tensor_count import (
    ContractionBudgets,
    _contract_component,
    _plan_components,
    count_placed_formula_by_contraction,
)


def run_tensor_count_suite() -> list[dict]:
    print('=== (Q) Connected replacement supports by finite Taylor contraction ===', flush=True)

    def details(result):
        return dict(
            plan=asdict(result.plan) if result.plan else None,
            components_counted=result.components_counted,
            boolean_vectors=result.boolean_vectors,
            factor_evaluations=result.factor_evaluations,
            jet_products=result.jet_products,
            peak_table_entries=result.peak_table_entries,
            peak_table_coefficients=result.peak_table_coefficients,
        )

    def path_placement(n, free=0):
        cnf = CNF(n + free, [[(v, False), (v + 1, False)] for v in range(n - 1)])
        placed = place_sparse_formula(
            cnf, build_structured_auxiliary(2, 2, 1), ((0, 0),) * cnf.n, ((0, 0),) * cnf.m
        )
        return cnf, placed

    cnf, placed = path_placement(13)
    legacy = count_placed_formula_by_support(placed)
    result = count_placed_formula_by_contraction(placed)
    expected = brute_count(cnf)
    ok = (
        result.counted
        and result.count == expected == 610
        and legacy.status == 'refused'
        and result.plan.max_bucket_variables == 2
        and result.boolean_vectors == 98
    )
    rows = [
        _row(
            'connected-support-n13',
            cnf.n,
            cnf.m,
            ok,
            'COUNTED',
            brute=expected,
            count=result.count,
            previous_support_status=legacy.status,
            same_replacement_polynomials=True,
            **details(result),
        )
    ]
    print(
        f'  [Q] connected n=13: #F={result.count}, local vectors={result.boolean_vectors}; {"COUNTED" if ok else "FAIL"}',
        flush=True,
    )

    cnf, placed = path_placement(1024, free=2)
    a, b = 0, 1
    for _ in range(1026):
        a, b = b, a + b
    result = count_placed_formula_by_contraction(placed)
    supplied = count_cnf_components(cnf, backend='tensor', max_component_variables=1024)
    ok = (
        result.counted
        and result.count == 4 * a
        and result.plan.largest_component == 1024
        and result.plan.free_variables == 2
        and result.plan.max_bucket_variables == 2
        and result.boolean_vectors == 8186
        and result.peak_table_entries == 4094
        and supplied.counted
        and supplied.count == result.count
        and supplied.boolean_vectors == 8186
    )
    rows.append(
        _row(
            'connected-path-n1024-free2',
            cnf.n,
            cnf.m,
            ok,
            'COUNTED',
            count_hex=hex(result.count) if result.counted else None,
            expected_count='4*Fibonacci(1026)',
            count_validation='independent recurrence for bit strings with no adjacent zeros',
            component_backend_count_matches=supplied.counted and supplied.count == result.count,
            component_backend_boolean_vectors=supplied.boolean_vectors,
            same_replacement_polynomials=True,
            **details(result),
        )
    )
    print(
        f'  [Q] connected n=1024 plus two free bits: #F=4*Fibonacci(1026), local vectors={result.boolean_vectors}; {"COUNTED" if ok else "FAIL"}',
        flush=True,
    )

    cnf = CNF(2, [[(0, False)], [(1, True)]])
    dense = place_formula(cnf, build_auxiliary_slice(2, 2, 2), ((0, 0), (1, 0)), ((1, 1), (0, 1)))
    placed = place_sparse_formula(
        cnf, build_structured_auxiliary(2, 2, 2), dense.variable_schedules, dense.clause_schedules
    )
    result = count_placed_formula_by_contraction(placed)
    components, free = _support_components(placed)
    plan, local_plans = _plan_components(components, free, ContractionBudgets())
    stats = dict(
        boolean_vectors=0,
        factor_evaluations=0,
        jet_products=0,
        peak_table_entries=0,
        peak_table_coefficients=0,
    )
    jet = _contract_component(components[0], local_plans[0], stats)
    expected_jet = eval_H_hasse_jet(cnf, dense.family, Fraction(1), plan.max_component_degree)
    ok = (
        result.counted
        and result.count == 1
        and len(components) == 1
        and any(jet[1:])
        and jet == expected_jet
    )
    rows.append(
        _row(
            'coupled-support-full-taylor',
            cnf.n,
            cnf.m,
            bool(ok),
            'COUNTED',
            brute=1,
            count=result.count,
            jet=list(map(str, jet)),
            independent_jet=list(map(str, expected_jet)),
            coefficient_validation_vectors=stats['boolean_vectors'],
            **details(result),
        )
    )
    print(
        f'  [Q] coupled replacements: nonconstant Taylor coefficients {"CHECKED" if ok else "FAIL"}',
        flush=True,
    )

    cnf = CNF(13, [[(v, False) for v in range(13)]])
    placed = place_sparse_formula(
        cnf, build_structured_auxiliary(2, 2, 1), ((0, 0),) * cnf.n, ((0, 0),)
    )
    result = count_placed_formula_by_contraction(placed)
    refused = clean_refusal(result.status, result.count, result.boolean_vectors)
    rows.append(
        _row(
            'wide-tensor-budget',
            cnf.n,
            cnf.m,
            refused,
            'REFUSED',
            count=None,
            boolean_vectors=result.boolean_vectors,
            reason=result.reason,
        )
    )
    print(
        f'  [Q] width-13 factor: {"REFUSED before numerical tables" if refused else "FAIL"}',
        flush=True,
    )
    return rows
