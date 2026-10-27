"""Exact finite contraction, including long connected supports and full jets."""

import random
from dataclasses import replace
from fractions import Fraction

import pytest

from counting_sat_algebra.cnf import CNF, brute_count, random_3cnf
from counting_sat_algebra.hermite import eval_H_hasse_jet
from counting_sat_algebra.replacement_matrix import build_auxiliary_slice, place_formula, step3_examples
from counting_sat_algebra.sparse_placement import place_sparse_formula
from counting_sat_algebra.structured_auxiliary import build_structured_auxiliary
from counting_sat_algebra.support_count import count_placed_formula_by_support
from counting_sat_algebra.tensor_count import ContractionBudgets, count_placed_formula_by_contraction


def _identity_placement(cnf):
    return place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 1),
                                ((0, 0),) * cnf.n, ((0, 0),) * cnf.m)


def _fib(index):
    a, b = 0, 1
    for _ in range(index):
        a, b = b, a + b
    return a


@pytest.mark.parametrize('example', step3_examples(), ids=lambda example: example.name)
def test_contraction_matches_independent_counts_on_the_same_dense_placement(example):
    result = count_placed_formula_by_contraction(example.placed)
    assert result.counted and result.count == brute_count(example.placed.cnf)
    assert result.count == count_placed_formula_by_support(example.placed).count
    assert result.boolean_vectors == result.plan.boolean_vector_bound
    assert result.jet_products == result.plan.jet_product_bound
    assert result.peak_table_entries == result.plan.peak_table_entries


def test_coupled_replacements_preserve_every_nonconstant_taylor_coefficient(monkeypatch):
    from counting_sat_algebra import tensor_count
    cnf = CNF(2, [[(0, False)], [(1, True)]])
    dense = place_formula(cnf, build_auxiliary_slice(2, 2, 2),
                          ((0, 0), (1, 0)), ((1, 1), (0, 1)))
    sparse = place_sparse_formula(cnf, build_structured_auxiliary(2, 2, 2),
                                  dense.variable_schedules, dense.clause_schedules)
    real = tensor_count._contract_component
    observed = []

    def capture(*args):
        jet = real(*args)
        observed.append(jet)
        return jet

    monkeypatch.setattr(tensor_count, '_contract_component', capture)
    result = count_placed_formula_by_contraction(sparse)
    assert result.counted and result.count == 1 and result.plan.components == 1
    assert len(observed) == 1 and any(observed[0][1:])
    assert observed[0] == eval_H_hasse_jet(cnf, dense.family, Fraction(1), result.plan.max_component_degree)


@pytest.mark.parametrize('cnf', [CNF(0, []), CNF(0, [[]]), CNF(5, []), CNF(3, [[]])])
def test_empty_formulas_empty_factors_and_declared_free_variables(cnf):
    result = count_placed_formula_by_contraction(_identity_placement(cnf))
    assert result.counted and result.count == brute_count(cnf)
    assert result.boolean_vectors == 0


def test_negative_repeated_and_tautological_literals_match_random_cube_oracles():
    rng = random.Random(15)
    for n in range(3, 8):
        cnf = random_3cnf(n, n + 2, rng=rng)
        cnf.clauses.extend([[(0, False), (0, True)], [(1, True), (1, True)]])
        cnf.n += 2
        result = count_placed_formula_by_contraction(_identity_placement(cnf))
        assert result.counted and result.count == brute_count(cnf)


def test_long_connected_support_is_counted_without_cube_or_component_backends(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('a global cube or earlier count backend was used')

    monkeypatch.setattr('counting_sat_algebra.cnf.brute_count', boom)
    monkeypatch.setattr('counting_sat_algebra.hermite.eval_H_hasse_jet', boom)
    monkeypatch.setattr('counting_sat_algebra.replacements.eval_H_at', boom)
    monkeypatch.setattr('counting_sat_algebra.common_count._substitutions', boom)
    monkeypatch.setattr('counting_sat_algebra.component_count.count_cnf_components', boom)
    monkeypatch.setattr('counting_sat_algebra.support_count._component_jet', boom)
    n = 1024
    cnf = CNF(n + 2, [[(v, False), (v + 1, False)] for v in range(n - 1)])
    placed = _identity_placement(cnf)
    assert count_placed_formula_by_support(placed).status == 'refused'
    result = count_placed_formula_by_contraction(placed)
    assert result.counted and result.count == 4 * _fib(n + 2)
    assert result.plan.components == result.components_counted == 1
    assert result.plan.largest_component == n and result.plan.free_variables == 2
    assert result.plan.max_bucket_variables == 2
    assert result.boolean_vectors == result.plan.boolean_vector_bound == 8186
    assert result.peak_table_entries == result.plan.peak_table_entries == 4094


@pytest.mark.parametrize('field', [
    'max_variables', 'max_component_degree', 'max_bucket_variables', 'max_boolean_vectors',
    'max_jet_products', 'max_table_entries', 'max_table_coefficients', 'max_planning_steps',
])
def test_every_budget_refuses_before_numerical_tables_are_created(field, monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('numerical tables created before all budgets were checked')
    monkeypatch.setattr('counting_sat_algebra.tensor_count._factor_tensor', boom)
    placed = step3_examples()[1].placed
    result = count_placed_formula_by_contraction(placed, budgets=replace(ContractionBudgets(), **{field: 0}))
    assert result.status == 'refused' and result.count is None and result.plan is None
    assert result.boolean_vectors == result.components_counted == result.factor_evaluations == 0


def test_invalid_designs_and_invalid_budget_types_are_rejected():
    placed = replace(step3_examples()[0].placed, measurement_rows=())
    result = count_placed_formula_by_contraction(placed)
    assert result.status == 'invalid-design' and result.count is None
    with pytest.raises(ValueError, match='nonnegative'):
        ContractionBudgets(max_table_entries=-1)
    with pytest.raises(TypeError):
        ContractionBudgets(max_boolean_vectors=0.5)


@pytest.mark.parametrize('wide_factor', [False, True])
def test_large_initial_or_induced_support_width_refuses_without_enumeration(wide_factor, monkeypatch):
    cnf = CNF(13, [[(v, False) for v in range(13)]]) if wide_factor else CNF(
        13, [[(u, False), (v, False)] for u in range(13) for v in range(u + 1, 13)],
    )
    placed = _identity_placement(cnf)

    def boom(*args, **kwargs):
        raise AssertionError('wide support was enumerated')
    monkeypatch.setattr('counting_sat_algebra.tensor_count._factor_tensor', boom)
    result = count_placed_formula_by_contraction(placed)
    assert result.status == 'refused' and result.count is None and result.boolean_vectors == 0
    assert ('factor support' if wide_factor else 'elimination bucket') in result.reason


def test_a_later_oversized_component_prevents_earlier_numeric_work(monkeypatch):
    cnf = CNF(14, [[(0, False)], [(v, False) for v in range(1, 14)]])
    placed = _identity_placement(cnf)

    def boom(*args, **kwargs):
        raise AssertionError('an early component was counted before global preflight')
    monkeypatch.setattr('counting_sat_algebra.tensor_count._factor_tensor', boom)
    result = count_placed_formula_by_contraction(placed)
    assert result.status == 'refused' and result.count is None
    assert result.components_counted == result.boolean_vectors == 0


def test_tight_operation_and_live_table_budgets_include_the_pending_output(monkeypatch):
    placed = _identity_placement(CNF(13, [[(v, False), (v + 1, False)] for v in range(12)]))
    first = count_placed_formula_by_contraction(placed)
    assert first.count == 610
    plan = first.plan
    tight = ContractionBudgets(
        max_variables=13, max_component_degree=0, max_bucket_variables=2,
        max_boolean_vectors=plan.boolean_vector_bound, max_jet_products=plan.jet_product_bound,
        max_table_entries=plan.peak_table_entries, max_table_coefficients=plan.peak_table_coefficients,
        max_planning_steps=plan.planning_steps,
    )
    assert count_placed_formula_by_contraction(placed, budgets=tight).count == 610

    def boom(*args, **kwargs):
        raise AssertionError('a below-bound budget reached numerical contraction')
    monkeypatch.setattr('counting_sat_algebra.tensor_count._factor_tensor', boom)
    for field in ('max_boolean_vectors', 'max_jet_products', 'max_table_entries',
                  'max_table_coefficients', 'max_planning_steps'):
        refused = count_placed_formula_by_contraction(placed, budgets=replace(tight, **{field: getattr(tight, field) - 1}))
        assert refused.status == 'refused' and refused.count is None and refused.boolean_vectors == 0
