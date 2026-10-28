"""Finite D.4.5 mode tables, crossing residues, and budget refusals."""

import json
from fractions import Fraction
from itertools import product
from math import prod

import pytest

from counting_sat_algebra.modular_modes import (
    ModeBudgets, ModeRequest, SupportRecord, admit_exponent_tuple, admit_support,
    build_mode_table, recover_crossing_residue, recover_network_with_modes,
    retained_degree_ok, support_admissible,
)
from counting_sat_algebra.modular_network import BooleanTensorNetwork, TensorVertex
from counting_sat_algebra.modular_recovery import (
    PacketLayout, cauchy_local_polynomial, embed_network,
)
from counting_sat_algebra.modular_ring import ModularRing, evaluate_polynomial


def _request(**overrides):
    payload = dict(
        a=2, t=1, edges=2, precision=4, row_symbol=0, column_symbol=1,
        real_columns=(0, 1, 2, 3),
    )
    payload.update(overrides)
    return ModeRequest(**payload)


def _rational_residue(a, t, row, column, zeta, precision):
    length = 1 << a
    polynomial = [Fraction(1, prod(-point for point in range(-length, length) if point))]
    for point in range(-length, length):
        if point == column - row:
            continue
        updated = [Fraction(0)] * (len(polynomial) + 1)
        for index, coefficient in enumerate(polynomial):
            updated[index] += -point * coefficient
            updated[index + 1] += coefficient
        polynomial = updated
    value = Fraction(t + (1 << (a + 1)) * zeta)
    total = Fraction(0)
    for coefficient in reversed(polynomial):
        total = total * value + coefficient
    modulus = 1 << precision
    return total.numerator * pow(total.denominator, -1, modulus) % modulus


def _brute(weights, precision):
    spans = [range((precision - 1) // weight + 1) for weight in weights]
    return [
        exponents for exponents in product(*spans)
        if sum(weight * exponent for weight, exponent in zip(weights, exponents)) < precision
    ]


def test_mode_table_matches_tile_weights_and_independent_tuples():
    witness = build_mode_table(_request())
    assert witness.modes_allocated and witness.count is None
    description = witness.description
    assert description is not None and description.tensor_evaluated is False
    assert [mode.weight for mode in description.modes] == [1, 2, 3, 1, 2, 3]
    assert [mode.source_group for mode in description.modes] == [1, 1, 1, 0, 0, 0]
    assert [mode.target_group for mode in description.modes] == [0, 0, 0, 1, 1, 1]
    assert description.modes[0].real_columns == (0, 2)
    weights = tuple(mode.weight for mode in description.modes)
    assert list(witness.tuples) == _brute(weights, 4)
    assert witness.partial_choices <= len(witness.tuples) * len(weights) + 1


def test_crossing_residue_matches_rational_oracle_and_cauchy_polynomial():
    request = _request()
    witness = recover_crossing_residue(request, 0, 2, 1)
    ring = ModularRing(4)
    direct = evaluate_polynomial(
        ring, cauchy_local_polynomial(ring, 2, 1, 0, 2), ring.scalar(1),
    )
    assert witness.count == direct[0] == _rational_residue(2, 1, 0, 2, 1, 4) == 8
    assert witness.coefficient_valuations == (3, 4, 5)
    assert witness.tensor_evaluations == 0
    assert witness.theorem_1_1_complete is witness.lean_proof_complete is False
    assert witness.recursive_contraction is False


@pytest.mark.parametrize('a,t,row,column,zeta,precision', [
    (1, 1, 0, 0, 0, 5),
    (1, 1, 0, 0, 1, 6),
    (1, -1, 1, 1, 2, 5),
    (2, 1, 0, 2, 1, 4),
    (2, 1, 1, 3, 3, 5),
    (2, 2, 0, 0, 1, 3),
    (2, -2, 1, 1, 1, 3),
])
def test_series_matches_rational_product_on_crossing_slots(a, t, row, column, zeta, precision):
    length = 1 << a
    coordinate = (t & -t).bit_length() - 1
    assert ((row >> coordinate) & 1) == ((column >> coordinate) & 1)
    witness = recover_crossing_residue(
        _request(
            a=a, t=t, edges=length, precision=precision,
            real_columns=tuple(range(length)),
        ),
        row, column, zeta,
    )
    assert witness.count == _rational_residue(a, t, row, column, zeta, precision)
    assert all(value >= index + 1 for index, value in enumerate(witness.coefficient_valuations))


def test_local_slot_and_equal_symbols_are_not_expanded():
    with pytest.raises(ValueError, match='local slots'):
        recover_crossing_residue(_request(), 0, 1, 1)
    local = build_mode_table(_request(row_symbol=3, column_symbol=3))
    assert local.description is not None and local.description.modes == ()
    assert local.tuples == ((),) and local.count is None
    with pytest.raises(ValueError, match='equal packet symbols'):
        recover_crossing_residue(_request(row_symbol=3, column_symbol=3), 0, 0, 1)


def test_support_gate_refuses_before_storing_an_inadmissible_degree():
    assert support_admissible(SupportRecord(0, 0, 0))
    assert support_admissible(SupportRecord(1, 0, 1))
    assert not support_admissible(SupportRecord(0, 0, 1))
    refused = admit_support(SupportRecord(2, 1, 4))
    assert refused.refused and refused.count is None and not refused.support_stored
    stored = admit_support(SupportRecord(2, 1, 3))
    assert stored.support_stored and stored.count is None
    assert all(
        retained_degree_ok(SupportRecord(deg_u, deg_z, deg_x), 4, 5)
        for deg_u in range(4) for deg_z in range(5) for deg_x in range(deg_u + deg_z + 1)
    )
    budget = admit_support(SupportRecord(0, 0, 0), budgets=ModeBudgets(max_support_records=0))
    assert budget.refused and budget.count is None and not budget.support_stored


def test_over_budget_exponent_returns_no_tuple():
    refused = admit_exponent_tuple((2, 1), (2, 0), 4)
    assert refused.refused and refused.tuples is None and refused.count is None
    assert not refused.tuples_allocated
    kept = admit_exponent_tuple((2, 1), (1, 1), 4)
    assert kept.tuples == ((1, 1),) and kept.count is None


@pytest.mark.parametrize('budgets', [
    ModeBudgets(max_modes=5),
    ModeBudgets(max_tuples=255),
    ModeBudgets(max_partial_choices=1536),
    ModeBudgets(max_precision_bits=3),
    ModeBudgets(max_packet_slots=3),
    ModeBudgets(max_real_columns=3),
])
def test_table_budgets_allocate_nothing_and_return_no_count(budgets):
    witness = build_mode_table(_request(), budgets=budgets)
    assert witness.refused and witness.count is None
    assert witness.description is None and witness.tuples is None
    assert not witness.modes_allocated and not witness.tuples_allocated
    assert witness.partial_choices == witness.tensor_evaluations == 0


def test_weight_class_refusal_precedes_the_mode_table():
    witness = build_mode_table(_request(edges=1))
    assert witness.refused and '6*M/D' in witness.reason
    assert witness.count is None and not witness.modes_allocated


def test_mode_refusal_does_not_call_step20_recovery():
    def forbidden(bits):
        raise AssertionError('recovery entry was called after a mode refusal')

    ring = ModularRing(5, 2)
    network = embed_network(BooleanTensorNetwork(ring, 2, (
        TensorVertex('row', 0, (0, 1), ring, forbidden),
        TensorVertex('column', 1, (0, 1), ring, forbidden),
    )), 0b10011)
    witness = recover_network_with_modes(
        network, PacketLayout(1, ((0, 1),)), 1,
        _request(a=1, edges=2, real_columns=(0, 1)),
        mode_budgets=ModeBudgets(max_modes=1),
    )
    assert witness.refused and witness.count is None and witness.base_value is None
    assert not witness.recovery_called and not witness.modes_allocated
    assert witness.port_pairs_scanned == witness.tensor_evaluations == 0


def test_admitted_modes_keep_the_step20_two_edge_value():
    ring = ModularRing(5, 2)

    def row(bits):
        return ring.from_base((1 + bits[0], 2 - bits[1]))

    def column(bits):
        return ring.from_base((1 + bits[1], -bits[0]))

    network = embed_network(BooleanTensorNetwork(ring, 2, (
        TensorVertex('row', 0, (0, 1), ring, row),
        TensorVertex('column', 1, (0, 1), ring, column),
    )), 0b10011)
    witness = recover_network_with_modes(
        network, PacketLayout(1, ((0, 1),)), 1,
        _request(a=1, real_columns=(0, 1)),
    )
    assert witness.base_value == (9, 4) and witness.count == 9
    assert witness.recovery_called and witness.tensor_evaluations == 0
    assert witness.recovery_tensor_evaluations == 448
    assert witness.port_pairs_scanned == 224
    assert witness.description is not None and witness.description.tensor_evaluated is False
    assert witness.recursive_contraction is False


def test_step22_report_keeps_theorem_flags_false(tmp_path, capsys):
    from counting_sat_algebra import __version__
    from counting_sat_algebra.__main__ import main

    output = tmp_path / 'step22.json'
    assert main(['verify', '--step22', '--output', str(output)]) == 0
    report = json.loads(output.read_text(encoding='utf-8'))
    assert report['version'] == __version__ and report['steps'] == [22]
    assert report['all_checks_passed'] and report['failed_step'] is None
    assert report['theorem_1_1_complete'] is report['lean_proof_complete'] is False
    assert len(report['rows']) == 13 and all(row['ok'] for row in report['rows'])
    assert {row['label'] for row in report['rows']} == {'COUNTED', 'CHECKED', 'REFUSED'}
    rows = {row['name']: row for row in report['rows']}
    assert rows['appendix-d-crossing-cauchy-residue']['count'] == 8
    assert rows['appendix-d-modes-then-step20-recovery']['base_value'] == [9, 4]
    assert rows['appendix-d-mode-count-budget']['count'] is None
    assert rows['appendix-d-exponent-over-budget']['tuples_allocated'] is False
    assert all(row['recursive_contraction'] is False for row in report['rows'])
    assert '13/13 checks passed' in capsys.readouterr().out
