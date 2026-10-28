"""Step 14 finite checks: disjoint replacement-support counts.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from dataclasses import asdict

from .protocol import _row, clean_refusal
from ..cnf import CNF
from ..structured_auxiliary import build_structured_auxiliary
from ..sparse_placement import place_sparse_formula, schedule_sparse_placement
from ..support_count import count_placed_formula_by_support


def run_support_count_suite() -> list[dict]:
    print('=== (P) Global placed polynomial by disjoint replacement supports ===', flush=True)
    cnf = CNF(16384, [[(v, bool(v % 2))] for v in range(16384)])
    outcome = schedule_sparse_placement(cnf, r=2, Q=3, R=2)
    if not outcome.accepted:
        return [
            _row(
                'n16384-support-contraction', cnf.n, cnf.m, False, 'COUNTED', reason=outcome.reason
            )
        ]
    result = count_placed_formula_by_support(outcome.placed, vector_budget=49152)
    ok = (
        result.counted
        and result.count == 1
        and result.plan.components == 10240
        and result.plan.largest_component == 4
        and result.plan.max_component_degree == 4
        and result.plan.boolean_vector_bound == 49152
        and result.boolean_vectors == 96
        and result.components_counted == result.plan.components
        and result.local_evaluations + result.memo_hits == result.components_counted
    )
    rows = [
        _row(
            'n16384-support-contraction',
            cnf.n,
            cnf.m,
            ok,
            'COUNTED',
            count=result.count,
            expected_count=1,
            count_validation='one assignment for independent alternating-polarity unit clauses',
            global_placement_checked=True,
            numerical_hypotheses_admitted=outcome.parameters.admitted,
            same_replacement_polynomials=True,
            memo_enabled=True,
            plan=asdict(result.plan) if result.plan else None,
            components_counted=result.components_counted,
            local_evaluations=result.local_evaluations,
            memo_hits=result.memo_hits,
            boolean_vectors=result.boolean_vectors,
            factor_evaluations=result.factor_evaluations,
        )
    ]
    print(
        f'  [P] n=16384: #F={result.count}; support size <= 4; vectors={result.boolean_vectors}/49152; {"COUNTED" if ok else "FAIL"}',
        flush=True,
    )

    # Reuse the exact same placement, but evaluate every support component.
    # This checks the advertised bound without relying on repeated local factors.
    uncached = count_placed_formula_by_support(outcome.placed, vector_budget=49152, memo_entries=0)
    uncached_ok = (
        uncached.counted
        and uncached.count == result.count == 1
        and uncached.plan == result.plan
        and uncached.boolean_vectors == 49152
        and uncached.components_counted == uncached.local_evaluations == 10240
        and uncached.memo_hits == 0
    )
    rows.append(
        _row(
            'n16384-support-no-memo',
            cnf.n,
            cnf.m,
            uncached_ok,
            'COUNTED',
            count=uncached.count,
            expected_count=1,
            count_validation='one assignment for independent alternating-polarity unit clauses',
            global_placement_checked=True,
            numerical_hypotheses_admitted=outcome.parameters.admitted,
            same_replacement_polynomials=True,
            memo_enabled=False,
            count_matches_memoized=uncached.counted
            and result.counted
            and uncached.count == result.count,
            plan=asdict(uncached.plan) if uncached.plan else None,
            components_counted=uncached.components_counted,
            local_evaluations=uncached.local_evaluations,
            memo_hits=uncached.memo_hits,
            boolean_vectors=uncached.boolean_vectors,
            factor_evaluations=uncached.factor_evaluations,
        )
    )
    print(
        f'  [P] same placement, no memo: #F={uncached.count}; vectors={uncached.boolean_vectors}/49152; {"COUNTED" if uncached_ok else "FAIL"}',
        flush=True,
    )

    cnf = CNF(13, [[(v, False), (v + 1, False)] for v in range(12)])
    placed = place_sparse_formula(
        cnf, build_structured_auxiliary(2, 2, 1), ((0, 0),) * 13, ((0, 0),) * 12
    )
    result = count_placed_formula_by_support(placed)
    refused = clean_refusal(result.status, result.count, result.boolean_vectors)
    rows.append(
        _row(
            'connected-support-budget',
            cnf.n,
            cnf.m,
            refused,
            'REFUSED',
            count=None,
            plan=asdict(result.plan),
            boolean_vectors=result.boolean_vectors,
            reason=result.reason,
        )
    )
    print(
        f'  [P] connected support n=13: {"REFUSED before enumeration" if refused else "FAIL"}',
        flush=True,
    )
    return rows
