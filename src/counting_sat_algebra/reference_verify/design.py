"""Step 9 finite checks: conditional expectation and the rational design.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction

from .protocol import _row
from ..cnf import CNF
from ..reference_parameters import check_reference_parameters
from ..scheduler import schedule_placement


def run_reference_design_suite() -> list[dict]:
    print('=== (K) Reference conditional expectation and rational design ===', flush=True)
    cnf = CNF(1024, [[(v, bool(v % 2))] for v in range(1024)])
    outcome = schedule_placement(cnf, r=2, Q=2, R=1)
    ok = (
        outcome.accepted
        and outcome.design_method == 'linear'
        and outcome.balance.admitted
        and outcome.multiplicity_ok
        and outcome.group_size_ok
        and outcome.final_potential <= outcome.initial_potential <= 2
    )
    rows = [
        _row(
            'admit-n1024-linear',
            cnf.n,
            cnf.m,
            bool(ok),
            'CHECKED',
            multiplicity=outcome.multiplicity,
            multiplicity_bound=str(outcome.multiplicity_bound),
            max_group_size=outcome.max_group_size,
            trials=outcome.trials,
            initial_potential=str(outcome.initial_potential),
            final_potential=str(outcome.final_potential),
            count_checked=False,
        )
    ]
    print(
        f'  [K] n=1024: placement {"CHECKED" if ok else "FAIL"}; multiplicity={outcome.multiplicity}, largest group={outcome.max_group_size}',
        flush=True,
    )
    for n, Q, R, admitted in ((1024, 2, 1, False), (16384, 3, 2, True)):
        parameters = check_reference_parameters(
            n, tuple((v,) for v in range(n)), r=2, Q=Q, R=R, tau=Fraction(1, 2)
        )
        correct = parameters.admitted == admitted
        rows.append(
            _row(
                f'lemma4.7-n{n}',
                n,
                n,
                correct,
                'CHECKED',
                numerical_hypotheses_admitted=parameters.admitted,
                tail=str(parameters.tail_probability),
                tail_limit=str(parameters.tail_limit),
                placement_checked=False,
                count_checked=False,
                reasons=list(parameters.reasons),
            )
        )
        print(
            f'  [K] n={n}, Q={Q}, R={R}: tail={parameters.tail_probability}, limit={parameters.tail_limit}; hypotheses {"admitted" if parameters.admitted else "rejected"}',
            flush=True,
        )
    return rows
