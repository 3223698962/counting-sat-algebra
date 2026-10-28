"""Step 22 finite checks: weighted modes before packet tensors.

Labels stay COUNTED, CHECKED, or REFUSED. A check that does not hold is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import prod

from .protocol import _row, stamp_completion
from ..modular_modes import (
    STEP22_SCOPE,
    ModeBudgets,
    ModeRequest,
    SupportRecord,
    admit_exponent_tuple,
    admit_support,
    build_mode_table,
    recover_crossing_residue,
    recover_network_with_modes,
    retained_degree_ok,
)
from ..modular_network import BooleanTensorNetwork, TensorVertex
from ..modular_recovery import PacketLayout, embed_network
from ..modular_ring import ModularRing


def _oracle_residue(a, t, row, column, zeta, precision):
    """Independent rational value of W_xy(t + alpha*zeta), reduced mod 2^precision."""
    length = 1 << a
    denominator = prod(-point for point in range(-length, length) if point)
    polynomial = [Fraction(1, denominator)]
    for point in range(-length, length):
        if point == column - row:
            continue
        updated = [Fraction(0)] * (len(polynomial) + 1)
        for index, coefficient in enumerate(polynomial):
            updated[index] += coefficient * -point
            updated[index + 1] += coefficient
        polynomial = updated
    value = Fraction(t + (1 << (a + 1)) * zeta)
    total = Fraction(0)
    for coefficient in reversed(polynomial):
        total = total * value + coefficient
    modulus = 1 << precision
    return total.numerator * pow(total.denominator, -1, modulus) % modulus


def _brute_tuples(weights, precision):
    spans = [range((precision - 1) // weight + 1) for weight in weights]
    return tuple(
        exponents
        for exponents in product(*spans)
        if sum(weight * exponent for weight, exponent in zip(weights, exponents)) < precision
    )


def _request(**overrides):
    payload = dict(
        a=2,
        t=1,
        edges=2,
        precision=4,
        row_symbol=0,
        column_symbol=1,
        real_columns=(0, 1, 2, 3),
    )
    payload.update(overrides)
    return ModeRequest(**payload)


def _parallel_network():
    ring = ModularRing(5, 2)

    def row(bits):
        return ring.from_base((1 + bits[0], 2 - bits[1]))

    def column(bits):
        return ring.from_base((1 + bits[1], -bits[0]))

    return BooleanTensorNetwork(
        ring,
        2,
        (
            TensorVertex('row', 0, (0, 1), ring, row),
            TensorVertex('column', 1, (0, 1), ring, column),
        ),
    )


def run_weighted_modes_suite() -> list[dict]:
    print('=== (X) Appendix D finite weighted modes ===', flush=True)
    rows = []
    table = build_mode_table(_request())
    modes = table.description.modes if table.description is not None else ()
    structural = (
        not table.refused
        and table.count is None
        and table.description is not None
        and table.description.tensor_evaluated is False
        and len(modes) == 6
        and [mode.weight for mode in modes] == [1, 2, 3, 1, 2, 3]
        and [mode.tile_origin for mode in modes] == [0, 0, 0, 1, 1, 1]
        and [mode.source_group for mode in modes] == [1, 1, 1, 0, 0, 0]
        and [mode.target_group for mode in modes] == [0, 0, 0, 1, 1, 1]
        and modes[0].real_columns == (0, 2)
        and modes[3].real_columns == (1, 3)
        and table.tuples is not None
        and len(table.tuples) == 18
        and table.partial_choices == 67
        and table.partial_choices <= 6 * 18 + 1
        and table.tensor_evaluations == 0
    )
    rows.append(
        _row(
            'appendix-d-weighted-mode-table',
            0,
            2,
            structural,
            'CHECKED',
            modes=len(modes),
            tuples=0 if table.tuples is None else len(table.tuples),
            partial_choices=table.partial_choices,
            tensor_evaluated=False,
            packet_cache_rows_origins_built=False,
            tensor_evaluations=0,
            count=None,
        )
    )

    local = build_mode_table(_request(row_symbol=0, column_symbol=0))
    rows.append(
        _row(
            'appendix-d-local-packet-empty-modes',
            0,
            2,
            not local.refused
            and local.description is not None
            and local.description.modes == ()
            and local.tuples == ((),)
            and local.count is None
            and local.tensor_evaluations == 0,
            'CHECKED',
            modes=0,
            tuples=1,
            partial_choices=local.partial_choices,
            tensor_evaluations=0,
            count=None,
        )
    )

    weights = tuple(mode.weight for mode in modes)
    brute = _brute_tuples(weights, 4)
    rows.append(
        _row(
            'appendix-d-exponent-tuple-enumeration',
            0,
            2,
            table.tuples == brute and table.partial_choices <= len(brute) * len(weights) + 1,
            'CHECKED',
            tuples=len(brute),
            partial_choices=table.partial_choices,
            partial_bound=len(brute) * len(weights) + 1,
            count=None,
        )
    )

    zero = admit_support(SupportRecord(0, 0, 0))
    linear = admit_support(SupportRecord(1, 0, 1))
    rectangle = all(
        retained_degree_ok(SupportRecord(deg_u, deg_z, deg_x), 4, 5)
        for deg_u in range(4)
        for deg_z in range(5)
        for deg_x in range(deg_u + deg_z + 1)
    )
    rows.append(
        _row(
            'appendix-d-weighted-support-admissibility',
            0,
            0,
            zero.support_stored
            and linear.support_stored
            and rectangle
            and zero.count is None
            and linear.count is None,
            'CHECKED',
            zero_admitted=True,
            linear_admitted=True,
            retained_rectangle=True,
            sqrtlog_incidence_polynomial=False,
            count=None,
        )
    )

    residue = recover_crossing_residue(_request(), 0, 2, 1)
    oracle = _oracle_residue(2, 1, 0, 2, 1, 4)
    rows.append(
        _row(
            'appendix-d-crossing-cauchy-residue',
            0,
            2,
            residue.count == oracle == 8
            and residue.coefficient_valuations == (3, 4, 5)
            and residue.tensor_evaluations == 0
            and residue.recovery_tensor_evaluations == 0,
            'COUNTED',
            count=residue.count,
            oracle=oracle,
            precision=4,
            coefficient_valuations=list(residue.coefficient_valuations),
            tensor_evaluations=0,
            packet_tensor_evaluated=False,
        )
    )

    network = embed_network(_parallel_network(), 0b10011)
    counted = recover_network_with_modes(
        network,
        PacketLayout(1, ((0, 1),)),
        1,
        _request(a=1, real_columns=(0, 1)),
    )
    rows.append(
        _row(
            'appendix-d-modes-then-step20-recovery',
            0,
            2,
            counted.recovery_called
            and counted.base_value == (9, 4)
            and counted.count == 9
            and counted.port_pairs_scanned == 224
            and counted.recovery_tensor_evaluations == 448
            and counted.tensor_evaluations == 0
            and counted.description is not None
            and counted.description.tensor_evaluated is False,
            'COUNTED',
            count=counted.count,
            base_value=list(counted.base_value or ()),
            port_pairs_scanned=counted.port_pairs_scanned,
            recovery_tensor_evaluations=counted.recovery_tensor_evaluations,
            tensor_evaluations=0,
            packet_tensor_evaluated=False,
        )
    )

    mode_budget = build_mode_table(_request(), budgets=ModeBudgets(max_modes=5))
    rows.append(
        _row(
            'appendix-d-mode-count-budget',
            0,
            2,
            mode_budget.refused
            and mode_budget.count is None
            and mode_budget.description is None
            and mode_budget.tuples is None
            and not mode_budget.modes_allocated
            and not mode_budget.tuples_allocated
            and mode_budget.tensor_evaluations == 0,
            'REFUSED',
            count=None,
            modes_allocated=False,
            tuples_allocated=False,
            tensor_evaluations=0,
            reason=mode_budget.reason,
        )
    )

    weight_class = build_mode_table(_request(edges=1))
    rows.append(
        _row(
            'appendix-d-weight-class-budget',
            0,
            1,
            weight_class.refused
            and weight_class.count is None
            and not weight_class.modes_allocated
            and weight_class.tuples is None,
            'REFUSED',
            count=None,
            modes_allocated=False,
            tuples_allocated=False,
            reason=weight_class.reason,
        )
    )

    tuple_budget = build_mode_table(_request(), budgets=ModeBudgets(max_tuples=255))
    rows.append(
        _row(
            'appendix-d-tuple-budget',
            0,
            2,
            tuple_budget.refused
            and tuple_budget.count is None
            and tuple_budget.tuples is None
            and not tuple_budget.tuples_allocated
            and not tuple_budget.modes_allocated,
            'REFUSED',
            count=None,
            modes_allocated=False,
            tuples_allocated=False,
            reason=tuple_budget.reason,
        )
    )

    partial_budget = build_mode_table(_request(), budgets=ModeBudgets(max_partial_choices=1536))
    rows.append(
        _row(
            'appendix-d-partial-choice-budget',
            0,
            2,
            partial_budget.refused
            and partial_budget.count is None
            and partial_budget.partial_choices == 0
            and not partial_budget.modes_allocated,
            'REFUSED',
            count=None,
            partial_choices=0,
            modes_allocated=False,
            reason=partial_budget.reason,
        )
    )

    inadmissible = admit_support(SupportRecord(0, 0, 1))
    rows.append(
        _row(
            'appendix-d-inadmissible-support',
            0,
            0,
            inadmissible.refused and inadmissible.count is None and not inadmissible.support_stored,
            'REFUSED',
            count=None,
            support_stored=False,
            reason=inadmissible.reason,
        )
    )

    support_budget = admit_support(
        SupportRecord(0, 0, 0), budgets=ModeBudgets(max_support_records=0)
    )
    rows.append(
        _row(
            'appendix-d-support-record-budget',
            0,
            0,
            support_budget.refused
            and support_budget.count is None
            and not support_budget.support_stored,
            'REFUSED',
            count=None,
            support_stored=False,
            reason=support_budget.reason,
        )
    )

    over_budget = admit_exponent_tuple((1, 2), (1, 1), 3)
    rows.append(
        _row(
            'appendix-d-exponent-over-budget',
            0,
            0,
            over_budget.refused
            and over_budget.count is None
            and over_budget.tuples is None
            and not over_budget.tuples_allocated,
            'REFUSED',
            count=None,
            tuples_allocated=False,
            weighted_degree=3,
            precision=3,
            reason=over_budget.reason,
        )
    )

    for row in rows:
        stamp_completion(row)
        print(f'  [X] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP22_SCOPE}', flush=True)
    return rows
