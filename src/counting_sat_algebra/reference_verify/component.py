"""Step 13 finite checks: exact counts of supplied CNFs by independent components.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction

from .protocol import _row
from ..cnf import CNF, brute_count
from ..component_count import count_cnf_components
from ..reference_parameters import check_reference_parameters


def run_component_count_suite() -> list[dict]:
    print('=== (O) Exact component counts on supplied CNFs ===', flush=True)
    cnf = CNF(8, [[(0, False), (1, False), (2, True)], [(3, False), (4, True), (5, False)]])
    result = count_cnf_components(cnf)
    ok = result.counted and result.count == brute_count(cnf) == 196
    rows = [
        _row(
            'component-demo-n8',
            8,
            2,
            ok,
            'COUNTED',
            brute=196,
            count=result.count,
            components=result.components,
            free_variables=result.free_variables,
            boolean_vectors=result.boolean_vectors,
        )
    ]
    print(
        f'  [O] n=8: #F={result.count}, declared free variables={result.free_variables}; {"COUNTED" if ok else "FAIL"}',
        flush=True,
    )
    cnf = CNF(
        32769, [[(3 * i, False), (3 * i + 1, True), (3 * i + 2, False)] for i in range(10923)]
    )
    parameters = check_reference_parameters(
        cnf.n,
        tuple(tuple(v for v, _neg in clause) for clause in cnf.clauses),
        r=2,
        Q=3,
        R=2,
        tau=Fraction(1, 2),
    )
    result = count_cnf_components(cnf)
    ok = result.counted and result.count == 7**10923 and parameters.admitted
    rows.append(
        _row(
            'disjoint-3cnf-n32769',
            cnf.n,
            cnf.m,
            ok,
            'COUNTED',
            count_hex=hex(result.count) if result.counted else None,
            expected_count='7**10923',
            count_validation='closed form for disjoint three-literal clauses',
            numerical_hypotheses_admitted=parameters.admitted,
            global_placement_checked=False,
            components=result.components,
            local_evaluations=result.local_evaluations,
            memo_hits=result.memo_hits,
            boolean_vectors=result.boolean_vectors,
        )
    )
    print(
        f'  [O] n=32769: #F=7**10923; local evaluations={result.local_evaluations}; memo hits={result.memo_hits}; {"COUNTED" if ok else "FAIL"}',
        flush=True,
    )
    return rows
