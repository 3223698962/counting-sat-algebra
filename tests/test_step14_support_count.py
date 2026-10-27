"""Actual placed-factor supports and local Taylor contraction, with independent oracles."""

from dataclasses import replace
from fractions import Fraction

import pytest

from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.common_count import count_placed_formula
from counting_sat_algebra.hermite import eval_H_hasse_jet
from counting_sat_algebra.replacement_matrix import build_auxiliary_slice, place_formula, step3_examples
from counting_sat_algebra.sparse_placement import place_sparse_formula, schedule_sparse_placement
from counting_sat_algebra.structured_auxiliary import build_structured_auxiliary
from counting_sat_algebra.support_count import _component_jet, count_placed_formula_by_support


@pytest.mark.parametrize('example', step3_examples(), ids=lambda example: example.name)
def test_same_placement_counts_match_shared_measurements_and_cube(example):
    result = count_placed_formula_by_support(example.placed, memo_entries=0)
    assert result.counted and result.count == brute_count(example.placed.cnf)
    assert result.count == count_placed_formula(example.placed).count
    assert result.boolean_vectors <= result.plan.boolean_vector_bound


def test_support_factors_merge_originally_independent_cnf_components_and_preserve_full_jets():
    # The original CNF has two independent units; the replacements couple them.
    cnf = CNF(2, [[(0, False)], [(1, True)]])
    dense = place_formula(cnf, build_auxiliary_slice(2, 2, 2),
                          ((0, 0), (1, 0)), ((1, 1), (0, 1)))
    sparse = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 2),
                                  dense.variable_schedules, dense.clause_schedules)
    clauses = tuple(tuple((neg, tuple(sorted(row.items())))
                          for (_v, neg), row in zip(clause, rows))
                    for clause, rows in zip(cnf.clauses, sparse.rows))
    stats = dict(boolean_vectors=0, factor_evaluations=0)
    jet = _component_jet(clauses, 2, 2, stats)
    assert jet == eval_H_hasse_jet(cnf, dense.family, Fraction(1), 2)
    assert any(jet[1:])  # This really checks a nonconstant H(T).
    result = count_placed_formula_by_support(sparse)
    assert result.counted and result.count == brute_count(cnf) == 1
    assert result.plan.components == 1 and result.plan.largest_component == 2


def test_factor_product_free_variables_and_disabled_memoization():
    cnf = CNF(8, [[(0, False), (1, True)], [(2, False), (3, True)]])
    sparse = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 1),
                                  ((0, 0),) * 8, ((0, 0),) * 2)
    cached = count_placed_formula_by_support(sparse)
    uncached = count_placed_formula_by_support(sparse, memo_entries=0)
    assert cached.count == uncached.count == brute_count(cnf) == 144
    assert cached.plan.components == 2 and cached.plan.free_variables == 4
    assert cached.local_evaluations == 1 and cached.memo_hits == 1
    assert uncached.local_evaluations == 2 and uncached.memo_hits == 0
    assert uncached.boolean_vectors == uncached.plan.boolean_vector_bound == 8


@pytest.mark.parametrize('cnf', [CNF(0, []), CNF(0, [[]]), CNF(3, []), CNF(3, [[]])])
def test_empty_formulas_and_declared_free_variables(cnf):
    sparse = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 1),
                                  ((0, 0),) * cnf.n, ((0, 0),) * cnf.m)
    result = count_placed_formula_by_support(sparse)
    assert result.counted and result.count == brute_count(cnf)


def test_refusals_happen_before_local_contraction_and_return_no_partial_count(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('enumeration started before the budget check')
    monkeypatch.setattr('counting_sat_algebra.support_count._component_jet', boom)
    placed = step3_examples()[1].placed
    for kwargs in ({'max_variables': 0}, {'max_component_variables': 0},
                   {'max_component_degree': 0}, {'vector_budget': 0}):
        result = count_placed_formula_by_support(placed, **kwargs)
        assert result.status == 'refused' and result.count is None and result.boolean_vectors == 0
    invalid = replace(step3_examples()[0].placed, measurement_rows=())
    result = count_placed_formula_by_support(invalid)
    assert result.status == 'invalid-design' and result.count is None
    with pytest.raises(ValueError, match='nonnegative'):
        count_placed_formula_by_support(placed, memo_entries=-1)


def test_admitted_global_placement_is_counted_without_common_sum_or_cube_oracles(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('an original-CNF backend or global-sum oracle was called')
    monkeypatch.setattr('counting_sat_algebra.cnf.brute_count', boom)
    monkeypatch.setattr('counting_sat_algebra.replacements.eval_H_at', boom)
    monkeypatch.setattr('counting_sat_algebra.hermite.eval_H_hasse_jet', boom)
    monkeypatch.setattr('counting_sat_algebra.component_count.count_cnf_components', boom)
    monkeypatch.setattr('counting_sat_algebra.common_count._substitutions', boom)
    cnf = CNF(16384, [[(v, bool(v % 2))] for v in range(16384)])
    outcome = schedule_sparse_placement(cnf, r=2, Q=3, R=2)
    result = count_placed_formula_by_support(outcome.placed)
    assert outcome.accepted and result.counted and result.count == 1
    assert result.plan.components == 10240 and result.plan.largest_component == 4
    assert result.plan.max_component_degree == 4 and result.plan.boolean_vector_bound == 49152
    assert result.components_counted == result.plan.components
    assert result.boolean_vectors == 96
    uncached = count_placed_formula_by_support(outcome.placed, memo_entries=0, vector_budget=49152)
    assert uncached.counted and uncached.count == result.count == 1
    assert uncached.plan == result.plan
    assert uncached.boolean_vectors == uncached.plan.boolean_vector_bound
    assert uncached.components_counted == uncached.local_evaluations == 10240
    assert uncached.memo_hits == 0

    # The advertised budget is independent of cache hits, even for this input.
    monkeypatch.setattr('counting_sat_algebra.support_count._component_jet', boom)
    refused = count_placed_formula_by_support(outcome.placed, vector_budget=49151)
    assert refused.status == 'refused' and refused.count is None
    assert refused.boolean_vectors == refused.local_evaluations == 0


def test_large_connected_replacement_support_remains_an_explicit_refusal():
    cnf = CNF(13, [[(v, False), (v + 1, False)] for v in range(12)])
    sparse = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 1),
                                  ((0, 0),) * 13, ((0, 0),) * 12)
    result = count_placed_formula_by_support(sparse)
    assert result.status == 'refused' and result.count is None and result.boolean_vectors == 0
    assert result.plan.largest_component == 13
