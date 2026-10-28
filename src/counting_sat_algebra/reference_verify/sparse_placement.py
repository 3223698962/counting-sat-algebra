"""Step 12 finite checks: sparse admitted placement and the common-count budget refusal.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from dataclasses import asdict

from .protocol import _row, clean_refusal
from ..cnf import CNF
from ..common_count import count_placed_formula
from ..sparse_placement import schedule_sparse_placement


def run_sparse_placement_suite() -> list[dict]:
    print('=== (N) Sparse placement satisfying all Lemma 4.7 hypotheses ===', flush=True)
    cnf = CNF(16384, [[(v, bool(v % 2))] for v in range(16384)])
    outcome = schedule_sparse_placement(cnf, r=2, Q=3, R=2)
    ok = outcome.accepted and outcome.placed.statistics.polynomial_entries == 59392
    rows = [
        _row(
            'admit-n16384-sparse',
            cnf.n,
            cnf.m,
            ok,
            'CHECKED',
            numerical_hypotheses_admitted=outcome.parameters.admitted,
            tail=str(outcome.parameters.tail_probability),
            storage=asdict(outcome.placed.statistics) if outcome.placed else None,
            multiplicity=outcome.placed.copies if outcome.placed else None,
            group_ranks=outcome.design.group_ranks if outcome.design else (),
            count_checked=False,
        )
    ]
    print(f'  [N] n=16384: hypotheses and sparse design {"CHECKED" if ok else "FAIL"}', flush=True)
    if not ok:
        print(outcome.reason, flush=True)
        return rows
    storage = outcome.placed.statistics
    print(
        f'  [N] polynomial entries={storage.polynomial_entries}; dense slots={storage.dense_polynomial_slots}; measurement entries={storage.measurement_entries}',
        flush=True,
    )
    counted = count_placed_formula(outcome.placed)
    refused = clean_refusal(counted.status, counted.count, counted.boolean_vectors)
    rows.append(
        _row(
            'n16384-common-budget',
            cnf.n,
            cnf.m,
            refused,
            'REFUSED',
            count=None,
            q=counted.plan.q,
            J=counted.plan.J,
            substitution_pair_bound_bits=counted.plan.substitution_pairs.bit_length(),
            measurement_tuple_bound_bits=[
                value.bit_length() for value in counted.plan.measurement_tuples
            ],
            local_cube_bound_bits=[
                [value.bit_length() for value in sizes] for sizes in counted.plan.local_cube_sizes
            ],
            boolean_vector_bound_bits=counted.plan.boolean_vector_bound.bit_length(),
            boolean_vectors=counted.boolean_vectors,
            reason=counted.reason,
        )
    )
    print(
        f'  [N] global common sum: {"REFUSED before enumeration" if refused else "FAIL"}',
        flush=True,
    )
    return rows
