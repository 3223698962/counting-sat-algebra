"""Finite Section 5.3–5.5 replacement sum. Not Section 6, Appendix D, or Theorem 1.1."""

import json

import pytest

from counting_sat_algebra import __version__
from counting_sat_algebra.__main__ import main
from counting_sat_algebra.algebra_count import (
    STEP18_SCOPE, AlgebraBudgets, count_replacement_sum, grouped_unit_design,
    h_degree_bound, incidence_cut, input_degree_bound, restrict_system, sum_supplied_design,
)
from counting_sat_algebra.reference_verify import run_section5_replacement_suite
from counting_sat_algebra.sparse_systems import (
    brute_common_zeros, construct_parity_polynomials, integer_lift, sparse_system,
)


def _oracle(system, ell):
    """Boolean cube sum in Z[V]/(V^ell). No replacement row and no Hermite step."""
    parity = construct_parity_polynomials(ell)
    active = [equation for equation in system.equations if not equation.is_zero]
    total = [0] * ell
    for mask in range(1 << system.n):
        bits = [(mask >> variable) & 1 for variable in range(system.n)]
        accumulator = [0] * ell
        accumulator[0] = 1
        for equation in active:
            lift = integer_lift(equation, bits)
            height = lift * (lift - 1)
            column = []
            height_power = 1
            for constant, slope in parity.coefficients:
                column.append((constant + slope * lift) * height_power)
                height_power *= height
            combined = [0] * ell
            for left, left_value in enumerate(accumulator):
                if left_value == 0:
                    continue
                for right, right_value in enumerate(column):
                    if left + right < ell and right_value:
                        combined[left + right] += left_value * right_value
            accumulator = combined
        for slot, value in enumerate(accumulator):
            total[slot] += value
    return tuple(total)


def test_degree_bounds_match_equations_5_2_and_5_5():
    assert input_degree_bound(2, 2, 5) == 20
    assert h_degree_bound(2, 2, 5, 1) == 20
    assert input_degree_bound(1, 1, 3) == 5
    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)
    cut = incidence_cut(system)
    assert cut.incidence_total == 8 and cut.occurrence_cap == 4 and cut.heavy == ()
    assert cut.neutral_scopes == ((0,), (1,), (2,), (3,))


def test_perturbed_replacement_residue_matches_the_cube_oracle():
    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)
    result = count_replacement_sum(system)
    oracle = _oracle(system, 5)
    assert result.counted and result.count == brute_common_zeros(system) == 4
    assert result.h0_coefficients == oracle
    assert sum(oracle) == 2308 and 2308 % 32 == 4
    assert result.sample_at_one != str(result.count)
    assert result.boolean_vectors == 399 and result.vector_bound == 714
    assert result.q == 0 and result.jet_length == 1 and result.design_Q == 21
    assert result.slope == '1/21' and result.neutral_factors == 4
    assert result.lemma_48_guard_passed is False and result.lemma_48_dimension_formed is False
    assert 'Theorem 1.1' in result.note


def test_positive_jet_order_and_a_split_partition_agree_with_the_oracle():
    system = sparse_system(2, [[[0]]], k=1, c=1)
    jets = count_replacement_sum(system, Q=3, R=1)
    assert jets.counted and jets.count == brute_common_zeros(system) == 2
    assert jets.q == 1 and jets.jet_length == 2 and jets.boolean_vectors == 105
    assert jets.h0_coefficients == _oracle(system, 3) == (2, 0, 0)
    assert jets.sample_at_one == '68/81'

    design = grouped_unit_design(system, ((0,), (1,)), (0,), R=1, Q=3)
    split = sum_supplied_design(system, design, ell=3)
    assert split.counted and split.count == 2 and split.h0_coefficients == (2, 0, 0)
    assert split.boolean_vectors == 90 and split.q == 1


def test_high_incidence_branches_add_to_the_unrestricted_count():
    system = sparse_system(4, [[[0]], [[0]], [[0]], [[0]], [[0]], [[0], [1]]], k=1, c=1)
    cut = incidence_cut(system)
    assert cut.heavy == (0,) and cut.occurrence_cap == 6 and cut.incidence_total == 11
    fixed = restrict_system(system, {0: 0})
    assert fixed.n == 3
    assert [equation.monomials for equation in fixed.equations if not equation.is_zero] == [(frozenset({0}),)]
    unsat = restrict_system(system, {0: 1})
    assert any(equation.is_one for equation in unsat.equations)
    result = count_replacement_sum(system)
    assert result.counted and result.count == brute_common_zeros(system) == 4
    assert result.branch_counts == (4, 0) and result.boolean_vectors == 110
    assert result.lemma_48_r == 1 and result.lemma_48_Q == 0 and result.lemma_48_R == 8
    assert result.lemma_48_dimension_formed is False


def test_constant_and_empty_systems_do_not_scan():
    empty = count_replacement_sum(sparse_system(3, [], k=1, c=1))
    assert empty.counted and empty.count == 8 and empty.boolean_vectors == 0
    assert empty.jets_allocated is False and empty.enumerated is False
    constant = count_replacement_sum(sparse_system(3, [[[]]], k=1, c=1))
    assert constant.counted and constant.count == 0 and constant.boolean_vectors == 0
    redundant = sparse_system(2, [[], [[0]]], k=1, c=1)
    assert count_replacement_sum(redundant).count == brute_common_zeros(redundant) == 2


def test_kernel_failure_returns_no_count_and_does_not_enumerate(monkeypatch):
    system = sparse_system(2, [[[0]]], k=1, c=1)
    design = grouped_unit_design(system, ((1,), (0,)), (0,), R=1, Q=2)

    def boom(*_args, **_kwargs):
        raise AssertionError('enumeration started after an invalid design')

    monkeypatch.setattr('counting_sat_algebra.algebra_count._enumerate_design', boom)
    result = sum_supplied_design(system, design, ell=3)
    assert result.status == 'invalid-design' and result.count is None
    assert result.boolean_vectors == 0 and result.enumerated is False and result.jets_allocated is False


def test_budget_refusals_happen_before_enumeration(monkeypatch):
    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)
    branched = sparse_system(4, [[[0]], [[0]], [[0]], [[0]], [[0]], [[0], [1]]], k=1, c=1)

    def boom(*_args, **_kwargs):
        raise AssertionError('enumeration started after a budget refusal')

    monkeypatch.setattr('counting_sat_algebra.algebra_count._enumerate_design', boom)
    for budgets in (
        AlgebraBudgets(max_boolean_vectors=1),
        AlgebraBudgets(max_coefficient_bits=1),
        AlgebraBudgets(max_jet_length=0),
        AlgebraBudgets(max_algebra_slots=1),
        AlgebraBudgets(max_accumulator_bits=1),
    ):
        result = count_replacement_sum(system, budgets)
        assert result.refused and result.count is None and result.boolean_vectors == 0
        assert result.jets_allocated is False and result.enumerated is False
    bits = count_replacement_sum(system, AlgebraBudgets(max_coefficient_bits=1))
    assert bits.coefficients_constructed is False
    accumulator = count_replacement_sum(system, AlgebraBudgets(max_accumulator_bits=1))
    assert accumulator.coefficients_constructed is True and accumulator.jets_allocated is False
    branches = count_replacement_sum(branched, AlgebraBudgets(max_branches=1))
    assert branches.refused and branches.branches == 2 and branches.branches_enumerated == 0
    assert 'no branch enumerated' in branches.reason


def test_small_ell_does_not_pretend_to_be_the_count():
    system = sparse_system(2, [[[0]]], k=1, c=1)
    with pytest.raises(ValueError, match='n\\+1'):
        count_replacement_sum(system, ell=2)


def test_scope_text_keeps_theorem_1_1_open():
    assert 'does not implement Theorem 1.1' in STEP18_SCOPE
    assert 'Section 6' in STEP18_SCOPE and 'Appendix D' in STEP18_SCOPE
    assert '2^{o(n)}' in STEP18_SCOPE


def test_step18_reference_rows_and_cli(tmp_path, capsys):
    rows = run_section5_replacement_suite()
    by_name = {row['name']: row for row in rows}
    assert by_name['replacement-degree-and-incidence']['label'] == 'CHECKED'
    assert by_name['replacement-degree-and-incidence']['input_degree_bound'] == 20
    counted = by_name['replacement-sum-n4']
    assert counted['label'] == 'COUNTED' and counted['count'] == counted['brute'] == 4
    assert counted['h0_coefficients'] == [-2, 6, 0, -192, 2496]
    assert counted['boolean_vectors'] == 399
    branched = by_name['high-incidence-branches-n4']
    assert branched['label'] == 'COUNTED' and branched['count'] == 4
    assert branched['branch_counts'] == [4, 0]
    guard = by_name['lemma-4.8-guard-on-branch']
    assert guard['label'] == 'CHECKED' and guard['dimension_formed'] is False
    assert guard['r'] == 1 and guard['Q'] == 0 and guard['R'] == 8
    for name in (
        'replacement-vector-budget', 'replacement-branch-budget', 'replacement-coefficient-budget',
    ):
        row = by_name[name]
        assert row['label'] == 'REFUSED' and row['count'] is None and row['boolean_vectors'] == 0
        assert row['jets_allocated'] is False
    assert all(row['ok'] for row in rows)
    assert {row['label'] for row in rows} == {'COUNTED', 'CHECKED', 'REFUSED'}

    output = tmp_path / 'step18.json'
    assert main(['verify', '--step18', '--output', str(output)]) == 0
    report = json.loads(output.read_text(encoding='utf-8'))
    text = capsys.readouterr().out
    assert report['version'] == __version__ == '0.19.0'
    assert report['steps'] == [18] and report['all_checks_passed'] is True
    assert report['theorem_1_1_complete'] is report['lean_proof_complete'] is False
    assert report['failed_step'] is None and report['wall_times_recorded'] is False
    assert 'Theorem 1.1' in report['scope'] and 'not implemented' in report['scope']
    assert 'Section 6' in report['scope'] and 'Appendix D' in report['scope']
    assert 'ALL TESTS PASSED' in text
    assert 'Theorem 1.1 is not implemented. This program does not refute ETH.' in text
    assert 'No Theorem 1.1, no 2^{o(n)} claim.' in text
    assert len(report['rows']) == len(rows)
