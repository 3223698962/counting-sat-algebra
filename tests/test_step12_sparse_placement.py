"""Sparse Section 4.5 placements versus the dense oracle and admitted inputs."""

from dataclasses import replace
from fractions import Fraction
import random

import pytest

from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.common_count import count_placed_formula
from counting_sat_algebra.hermite import eval_H_hasse_jet
from counting_sat_algebra.linear_algebra import (
    SparseEliminationBudgetExceeded, SparseRationalRowSpan,
)
from counting_sat_algebra.replacement_matrix import (
    build_auxiliary_slice, place_formula, step3_examples, vector_in_row_span,
)
from counting_sat_algebra.replacements import ReplacementFamily
from counting_sat_algebra.sparse_placement import (
    SparsePlacementBudgetExceeded, SparsePlacementBudgets, check_sparse_design,
    place_sparse_formula, schedule_sparse_placement,
)
from counting_sat_algebra.structured_auxiliary import build_structured_auxiliary


def _dense_family(placed):
    return ReplacementFamily(
        placed.cnf.n, placed.auxiliary.R - 1,
        [[[list(row.get(v, (Fraction(0),))) for v in range(placed.cnf.n)] for row in clause] for clause in placed.rows],
        [[v for v, _neg in clause] for clause in placed.cnf.clauses],
        [[neg for _v, neg in clause] for clause in placed.cnf.clauses],
    )


@pytest.mark.parametrize('example', step3_examples(), ids=lambda example: example.name)
def test_sparse_rows_measurements_and_jets_match_dense_placements(example):
    dense = example.placed
    auxiliary = build_structured_auxiliary(dense.auxiliary.r, dense.auxiliary.Q, dense.auxiliary.R, dense.auxiliary.A0)
    sparse = place_sparse_formula(dense.cnf, auxiliary, dense.variable_schedules, dense.clause_schedules)
    for e, clause in enumerate(dense.family.rows):
        for occ, row in enumerate(clause):
            assert sparse.rows[e][occ] == {v: tuple(poly) for v, poly in enumerate(row) if any(poly)}
    expected = {tuple(row) for row in dense.measurement_rows if any(row)}
    actual = {tuple(row.get(v, Fraction(0)) for v in range(dense.cnf.n)) for row in sparse.measurement_rows}
    assert actual == expected
    assert sparse.copies == dense.copies and check_sparse_design(sparse).ok
    result = count_placed_formula(sparse)
    assert result.counted and result.count == brute_count(dense.cnf)
    for t, jet in enumerate(result.jets, 1):
        assert jet == eval_H_hasse_jet(dense.cnf, dense.family, Fraction(t), result.plan.q)


def test_sparse_row_span_matches_dense_oracle_and_budgets():
    rng = random.Random(12)
    for width in range(7):
        rows = [[Fraction(rng.randrange(-2, 3), rng.randrange(1, 4)) for _ in range(width)] for _ in range(5)]
        span = SparseRationalRowSpan([{i: value for i, value in enumerate(row) if value} for row in rows], width)
        for _ in range(10):
            target = [Fraction(rng.randrange(-2, 3)) for _ in range(width)]
            assert span.contains({i: value for i, value in enumerate(target) if value}) == vector_in_row_span(rows, target)
    # An enormous declared width does not allocate a dense row.
    span = SparseRationalRowSpan([{10**9: Fraction(1), 10**9 + 1: Fraction(2)}], 10**9 + 2)
    assert span.rank == 1 and span.stored_entries == 2
    assert span.contains({10**9: Fraction(3), 10**9 + 1: Fraction(6)})
    with pytest.raises(SparseEliminationBudgetExceeded):
        SparseRationalRowSpan([{0: Fraction(1)}], 1, max_entries=0)
    with pytest.raises(ValueError):
        SparseRationalRowSpan([{1: Fraction(1)}], 1)


def test_repeated_variables_in_different_occurrence_coordinates_coalesce_exactly():
    cnf = CNF(3, [[(0, False), (1, False)], [(0, False), (2, True)], [(0, True), (1, False)]])
    variables, clauses = ((0, 0), (1, 0), (1, 1)), ((1, 1), (1, 0), (0, 1))
    dense = place_formula(cnf, build_auxiliary_slice(2, 2, 2), variables, clauses)
    sparse = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 2), variables, clauses)
    for e, clause in enumerate(dense.family.rows):
        for occ, row in enumerate(clause):
            assert sparse.rows[e][occ] == {v: tuple(poly) for v, poly in enumerate(row) if any(poly)}
    result = count_placed_formula(sparse)
    assert result.counted and result.count == brute_count(cnf)
    for t, jet in enumerate(result.jets, 1):
        assert jet == eval_H_hasse_jet(cnf, dense.family, Fraction(t), result.plan.q)


def test_beyond_dense_slice_cap_counts_a_crossing_clause():
    cnf = CNF(3, [[(0, False), (1, True), (2, False)]])
    auxiliary = build_structured_auxiliary(2, 5, 3)
    sparse = place_sparse_formula(cnf, auxiliary,
                                 ((0, 1, 0, 0, 1), (1, 0, 1, 0, 0), (0, 0, 1, 1, 0)),
                                 ((1, 1, 1, 0, 1),))
    assert auxiliary.dimension == 1024 and check_sparse_design(sparse).ok
    result = count_placed_formula(sparse)
    assert result.counted and result.count == brute_count(cnf) == 7
    family = _dense_family(sparse)
    for t, jet in enumerate(result.jets, 1):
        assert jet == eval_H_hasse_jet(cnf, family, Fraction(t), result.plan.q)


def test_admitted_16384_variables_use_sparse_storage_and_refuse_count_before_enumeration(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('a dense matrix, dense placement or Boolean enumeration was used')

    monkeypatch.setattr('counting_sat_algebra.replacement_matrix.build_auxiliary_slice', boom)
    monkeypatch.setattr('counting_sat_algebra.replacement_matrix.place_formula', boom)
    monkeypatch.setattr('counting_sat_algebra.common_count._substitutions', boom)
    cnf = CNF(16384, [[(v, bool(v % 2))] for v in range(16384)])
    outcome = schedule_sparse_placement(cnf, r=2, Q=3, R=2)
    assert outcome.accepted and outcome.parameters.admitted
    assert outcome.parameters.tail_probability == Fraction(5, 32)
    assert outcome.placed.copies == 256 and outcome.design.max_group_size == 8192
    assert outcome.design.group_ranks == ((2048, 2048),) * 3
    statistics = outcome.placed.statistics
    assert statistics.dense_polynomial_slots == 536870912
    assert statistics.polynomial_entries == 59392 and statistics.measurement_entries == 14336
    assert statistics.pullback_candidates <= 65536
    assert outcome.schedules.final_potential <= outcome.schedules.initial_potential <= 2
    result = count_placed_formula(outcome.placed)
    assert result.status == 'refused' and result.count is None and result.boolean_vectors == 0
    assert result.plan.q == 5461 and result.plan.J == 5462
    assert result.plan.boolean_vector_bound.bit_length() > 8192
    assert 'no enumeration started' in result.reason


@pytest.mark.parametrize('budget', ['max_modes', 'max_mode_coordinates', 'max_pullback_candidates',
                                   'max_coefficient_entries', 'max_measurement_entries'])
def test_construction_budgets_refuse_before_modes_or_support(monkeypatch, budget):
    def boom(*args, **kwargs):
        raise AssertionError('an enumeration started before its budget check')

    auxiliary = build_structured_auxiliary(2, 1, 1)
    monkeypatch.setattr(type(auxiliary), 'measurement_modes', boom)
    monkeypatch.setattr(type(auxiliary), 'replacement_support', boom)
    with pytest.raises(SparsePlacementBudgetExceeded, match='construction not started'):
        place_sparse_formula(CNF(1, [[(0, False)]]), auxiliary, ((0,),), ((0,),),
                             budgets=replace(SparsePlacementBudgets(), **{budget: 0}))


def test_reference_and_trial_gates_refuse_without_search(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('schedule search was run')

    monkeypatch.setattr('counting_sat_algebra.sparse_placement.choose_schedules', boom)
    refused = schedule_sparse_placement(CNF(1024, [[(v, False)] for v in range(1024)]), r=2, Q=2, R=1)
    assert not refused.accepted and refused.placed is None
    assert 'tail' in refused.reason and refused.schedules is None
    cnf = CNF(16384, [[(v, False)] for v in range(16384)])
    refused = schedule_sparse_placement(cnf, r=2, Q=3, R=2, budgets=replace(SparsePlacementBudgets(), max_trials=1))
    assert refused.status == 'refused' and refused.placed is None and refused.parameters.admitted
    assert refused.schedules is None and 'candidate trial' in refused.reason


def test_sparse_validation_rejects_wrong_restoration_metadata_and_measurements():
    example = step3_examples()[0].placed
    auxiliary = build_structured_auxiliary(2, 1, 1)
    sparse = place_sparse_formula(example.cnf, auxiliary, example.variable_schedules, example.clause_schedules)
    assert not check_sparse_design(replace(sparse, measurement_rows=())).ok
    assert not check_sparse_design(replace(sparse, rows=(({},), sparse.rows[1]))).ok
    assert not check_sparse_design(replace(sparse, occurrence_mass=0)).ok
    assert not check_sparse_design(replace(sparse, variable_schedules=())).ok
    fractional = replace(sparse, measurement_rows=({0: Fraction(1, 2)},))
    assert any('integral' in problem for problem in check_sparse_design(fractional).problems)
    with pytest.raises(ValueError, match='nonnegative'):
        SparsePlacementBudgets(max_modes=-1)


@pytest.mark.parametrize('cnf', [CNF(0, []), CNF(0, [[]]), CNF(3, []), CNF(3, [[]])])
def test_empty_sparse_formulas_and_free_variables(cnf):
    auxiliary = build_structured_auxiliary(2, 2, 1)
    placed = place_sparse_formula(cnf, auxiliary, ((0, 0),) * cnf.n, ((0, 0),) * cnf.m)
    result = count_placed_formula(placed)
    assert result.counted and result.count == brute_count(cnf)
