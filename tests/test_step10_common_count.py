"""Section 3's streamed restricted sums against independent cube jets."""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction

import pytest

from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.common_count import count_placed_formula
from counting_sat_algebra.hermite import eval_H_hasse_jet
from counting_sat_algebra.replacement_matrix import (
    build_auxiliary_slice, place_formula, step3_examples,
)
from counting_sat_algebra.reference_parameters import check_reference_parameters


@pytest.mark.parametrize("example", step3_examples(), ids=lambda example: example.name)
def test_shared_measurement_jets_match_the_full_cube(example):
    placed = example.placed
    result = count_placed_formula(placed)
    assert result.counted and result.count == brute_count(placed.cnf)
    assert placed.auxiliary.Q * result.plan.J > result.plan.degree_cap
    assert result.boolean_vectors <= result.plan.boolean_vector_bound
    for t, jet in enumerate(result.jets, 1):
        assert jet == eval_H_hasse_jet(placed.cnf, placed.family, Fraction(t), result.plan.q)


def test_second_order_mobius_recovery_negative_literals_and_fixed_bit_restrictions():
    sl = build_auxiliary_slice(2, 2, 2)
    cnf = CNF(2, [[(0, False)], [(1, True)], [(0, True), (1, False)]])
    placed = place_formula(cnf, sl, ((0, 0), (1, 0)), ((1, 1),) * 3)
    result = count_placed_formula(placed)
    assert result.counted and result.plan.q == 2
    assert result.count == brute_count(cnf) == 0
    for t, jet in enumerate(result.jets, 1):
        assert jet == eval_H_hasse_jet(cnf, placed.family, Fraction(t), 2)


@pytest.mark.parametrize("cnf", [CNF(0, []), CNF(0, [[]]), CNF(3, []), CNF(3, [[]])])
def test_empty_formula_empty_clause_and_free_variables(cnf):
    sl = build_auxiliary_slice(2, 2, 1)
    placed = place_formula(cnf, sl, ((0, 0),) * cnf.n, ((0, 0),) * cnf.m)
    result = count_placed_formula(placed)
    assert result.counted and result.count == brute_count(cnf)


def test_budget_refuses_before_any_substitution_or_local_cube(monkeypatch):
    placed = step3_examples()[1].placed

    def boom(*args, **kwargs):
        raise AssertionError("an enumeration started")

    monkeypatch.setattr('counting_sat_algebra.common_count._substitutions', boom)
    for kwargs in ({'vector_budget': 0}, {'pair_budget': 0}, {'jet_budget': 0}):
        result = count_placed_formula(placed, **kwargs)
        assert result.status == 'refused' and result.count is None
        assert result.jets == () and result.boolean_vectors == result.measurement_tuples == 0
    with pytest.raises(ValueError, match='nonnegative'):
        count_placed_formula(placed, vector_budget=-1)


def test_invalid_design_returns_no_count_and_does_not_enumerate(monkeypatch):
    placed = replace(step3_examples()[0].placed, measurement_rows=())

    def boom(*args, **kwargs):
        raise AssertionError("an enumeration started")

    monkeypatch.setattr('counting_sat_algebra.common_count._substitutions', boom)
    result = count_placed_formula(placed)
    assert result.status == 'invalid-design' and result.count is None
    assert 'outside the measurement span' in result.reason


def test_reference_tail_is_checked_exactly_without_matrix_allocation(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError("a matrix was allocated")

    monkeypatch.setattr('counting_sat_algebra.replacement_matrix.build_auxiliary_slice', boom)
    # Q=2, R=1 meets the balance gate at n=1024, but not the tail gate.
    first = check_reference_parameters(1024, tuple((v,) for v in range(1024)), r=2, Q=2, R=1, tau=Fraction(1, 2))
    assert first.balance.admitted and first.balance.degree_gap_ok
    assert first.tail_probability == Fraction(7, 16) and first.tail_limit == Fraction(1, 4)
    assert not first.admitted and any('tail' in reason for reason in first.reasons)
    # All numerical Lemma 4.7 hypotheses hold at this finite parameter point.
    second = check_reference_parameters(16384, tuple((v,) for v in range(16384)), r=2, Q=3, R=2, tau=Fraction(1, 2))
    assert second.tail_probability == Fraction(5, 32)
    assert second.admitted and second.reasons == ()
    invalid = check_reference_parameters(1024, ((0,),) * 1024, r=2, Q=2, R=1, tau=Fraction(0), A0=-1)
    assert not invalid.admitted
    assert any('interval' in reason for reason in invalid.reasons)


def test_counter_has_no_global_cube_oracle_dependency(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('the full-cube oracle was called')

    monkeypatch.setattr('counting_sat_algebra.cnf.brute_count', boom)
    monkeypatch.setattr('counting_sat_algebra.replacements.eval_H_at', boom)
    monkeypatch.setattr('counting_sat_algebra.hermite.eval_H_hasse_jet', boom)
    result = count_placed_formula(step3_examples()[1].placed)
    assert result.counted and result.count == 1
