"""Step 10 finite checks: shared-measurement counting and Taylor recovery.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction
from dataclasses import asdict

from .protocol import _row
from ..cnf import brute_count
from ..common_count import count_placed_formula
from ..hermite import eval_H_hasse_jet
from ..replacement_matrix import step3_examples


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
        rows.append(
            _row(
                example.name,
                placed.cnf.n,
                placed.cnf.m,
                bool(ok),
                'COUNTED',
                brute=expected,
                count=result.count,
                jets=[list(map(str, jet)) for jet in result.jets],
                plan=asdict(result.plan),
                boolean_vectors=result.boolean_vectors,
                measurement_tuples=result.measurement_tuples,
                factor_evaluations=result.factor_evaluations,
            )
        )
        print(
            f'  [L] {example.name}: #F={result.count}, pairs={result.plan.substitution_pairs}, local vectors={result.boolean_vectors}/{result.plan.boolean_vector_bound}; {"COUNTED" if ok else "FAIL"}',
            flush=True,
        )
    return rows
