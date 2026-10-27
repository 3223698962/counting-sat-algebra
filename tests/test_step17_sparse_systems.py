"""Finite Section 5.1–5.2 checks. These are not Theorem 2.2 or Theorem 1.1."""

import decimal
import json

import pytest

from counting_sat_algebra import __version__
from counting_sat_algebra.__main__ import main
from counting_sat_algebra.reference_verify import run_section5_suite
from counting_sat_algebra.sparse_systems import (
    STEP17_SCOPE, BudgetRefusal, Section5Budgets, beigel_tarui_complement,
    brute_common_zeros, build_truncated_product, construct_parity_polynomials,
    count_common_zeros, count_retained_terms, evaluate_truncated_product, integer_lift,
    monomial_allowance, sparse_polynomial, sparse_system, truncated_indicator_residue,
)


def _decimal_allowance(n, c):
    if n < 2:
        return 0
    decimal.getcontext().prec = 80
    log2 = decimal.Decimal(n).ln() / decimal.Decimal(2).ln()
    value = log2 ** c
    nearest = value.to_integral_value(rounding=decimal.ROUND_HALF_EVEN)
    if abs(value - nearest) < decimal.Decimal('1e-10'):
        value = nearest
    return int(value.to_integral_value(rounding=decimal.ROUND_FLOOR))


def test_monomial_allowance_matches_exact_powers_and_a_decimal_cross_check():
    assert monomial_allowance(1, 1) == 0
    assert monomial_allowance(2, 4) == 1
    assert monomial_allowance(1 << 10, 3) == 1000
    assert monomial_allowance(16, 2) == 16
    for n in range(1, 65):
        for c in (1, 2, 3):
            assert monomial_allowance(n, c) == _decimal_allowance(n, c)


def test_normalization_cancels_duplicates_and_collapses_powers():
    polynomial = sparse_polynomial(4, [[0, 0, 1], [0, 1], [2], [2]], degree_bound=2, monomial_bound=4)
    assert polynomial.monomials == ()
    assert polynomial.is_zero and integer_lift(polynomial, [1, 1, 1, 1]) == 0
    one = sparse_polynomial(2, [[]], degree_bound=0, monomial_bound=1)
    assert one.is_one and one.degree == 0
    assert integer_lift(one, [0, 0]) == 1
    kept = sparse_polynomial(4, [[1], [0], [0], [0]], degree_bound=1, monomial_bound=4)
    assert kept.monomials == (frozenset({0}), frozenset({1}))


def test_sparsity_and_degree_hypotheses_reject_the_input_list():
    with pytest.raises(ValueError, match='3 monomials'):
        sparse_system(4, [[[0], [1], [2]]], k=1, c=1)
    with pytest.raises(ValueError, match='degree 3'):
        sparse_system(4, [[[0, 1, 2]]], k=2, c=1)
    with pytest.raises(ValueError, match='outside'):
        sparse_polynomial(3, [[3]])
    system = sparse_system(4, [[[0], [1]]], k=1, c=1)
    assert system.meets_theorem_size is False
    wide = sparse_system(20, [[[0, 1], [2]]], k=2, c=1)
    assert wide.meets_theorem_size is True
    assert wide.m_test == 1


def test_boolean_lift_has_the_same_parity_as_the_f2_equation():
    system = sparse_system(4, [[[0, 1], [2], []], [[3]]], k=2, c=2)
    for mask in range(16):
        assignment = [(mask >> variable) & 1 for variable in range(4)]
        for equation in system.equations:
            lift = integer_lift(equation, assignment)
            direct = 0
            for monomial in equation.monomials:
                bit = 1
                for variable in monomial:
                    bit &= assignment[variable]
                direct ^= bit
            assert lift % 2 == direct


def test_lemma_5_1_agrees_on_integers_and_reproduces_the_binomial_form():
    for ell in range(1, 9):
        first = construct_parity_polynomials(ell)
        second = construct_parity_polynomials(ell)
        raw = beigel_tarui_complement(ell)
        assert first.coefficients == second.coefficients
        assert raw[0] == 1 and all(coefficient == 0 for coefficient in raw[1:ell])
        assert len(raw) == 2 * ell
        for y in range(-8, 17):
            assert first.indicator_agrees(y)
            assert first.evaluate(y) % (1 << ell) == (0 if y % 2 else 1)


def test_ell_2_matches_the_displayed_polynomial_and_drops_the_double_height_term():
    parity = construct_parity_polynomials(2)
    assert parity.coefficients == ((1, -1), (-1, 2))
    assert beigel_tarui_complement(2) == (1, 0, -3, 2)
    product = build_truncated_product(2, 2)
    assert [term.powers for term in product.terms] == [(0, 0), (0, 1), (1, 0)]
    assert product.shape_ok()
    assert all(term.linear_factors <= 2 and term.quadratic_factors <= 1 for term in product.terms)

    def height(y):
        return y * (y - 1)

    def displayed(y1, y2):
        return (
            (1 - y1) * (1 - y2)
            + (2 * y1 - 1) * height(y1) * (1 - y2)
            + (1 - y1) * (2 * y2 - 1) * height(y2)
        )

    def omitted(y1, y2):
        return (2 * y1 - 1) * (2 * y2 - 1) * height(y1) * height(y2)

    for y1 in range(-5, 8):
        for y2 in range(-5, 8):
            got = evaluate_truncated_product(product, parity, (y1, y2))
            full = parity.evaluate(y1) * parity.evaluate(y2)
            assert got == displayed(y1, y2)
            assert got + omitted(y1, y2) == full
            assert omitted(y1, y2) % 4 == 0
            assert truncated_indicator_residue(parity, (y1, y2)) == got % 4


def test_retained_term_count_is_the_bounded_binomial_and_refuses_without_tuples():
    assert count_retained_terms(0, 5, 1) == 1
    assert count_retained_terms(2, 2, 3) == 3
    assert count_retained_terms(2, 5, 100) == 15
    assert count_retained_terms(8, 8, 100) is None
    product = build_truncated_product(3, 3, Section5Budgets(max_terms=10, max_ell=3))
    assert len(product.terms) == 10 and product.shape_ok()
    with pytest.raises(BudgetRefusal, match='term table not allocated'):
        build_truncated_product(3, 3, Section5Budgets(max_terms=9, max_ell=3))
    with pytest.raises(BudgetRefusal, match='parity coefficients not constructed'):
        construct_parity_polynomials(30)


def test_truncated_count_matches_the_f2_oracle_on_tiny_systems():
    cases = [
        (4, [[[0, 1], [2]], [[3], []]], 2, 1, 4),
        (4, [[[0], [1]]], 1, 1, 8),
        (4, [], 2, 1, 16),
        (4, [[[]]], 2, 1, 0),
        (4, [[], [[0]]], 1, 1, 8),
        (4, [[[0], [0]]], 1, 1, 16),
        (3, [[[0, 1]]], 2, 1, 6),
        (8, [[[0], [1], [2]]], 2, 1, 128),
        (5, [[[0, 1], [2]], [[3], [4]]], 2, 1, 8),
        (4, [[[0, 1], [2]], [[3], []], []], 2, 1, 4),
    ]
    for n, equations, k, c, expected in cases:
        system = sparse_system(n, equations, k=k, c=c)
        result = count_common_zeros(system)
        assert result.counted and result.count == expected == brute_common_zeros(system)
        if equations and any(equation and equation != [[]] and equation != [[0], [0]] for equation in equations):
            assert result.enumerated and result.boolean_vectors == 1 << n
            assert result.term_table_allocated
        else:
            assert result.boolean_vectors == 0 and result.term_table_allocated is False


def test_pointwise_residue_is_the_common_zero_indicator():
    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)
    parity = construct_parity_polynomials(5)
    product = build_truncated_product(2, 5)
    assert all(term.lift_degree <= 2 + 2 * 4 for term in product.terms)
    total = 0
    for mask in range(16):
        assignment = [(mask >> variable) & 1 for variable in range(4)]
        lifts = [integer_lift(equation, assignment) for equation in system.equations]
        residue = evaluate_truncated_product(product, parity, lifts) % 32
        vector = truncated_indicator_residue(parity, lifts)
        indicator = 0 if any(lift % 2 for lift in lifts) else 1
        assert residue == vector == indicator
        total += indicator
    assert total == 4


def test_budget_refusals_allocate_nothing_and_return_no_count(monkeypatch):
    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)

    def boom(*args, **kwargs):
        raise AssertionError('a coefficient table or term table was allocated')

    monkeypatch.setattr('counting_sat_algebra.sparse_systems.construct_parity_polynomials', boom)
    monkeypatch.setattr('counting_sat_algebra.sparse_systems.build_truncated_product', boom)
    for budgets in (
        Section5Budgets(max_terms=2),
        Section5Budgets(max_ell=4),
        Section5Budgets(max_boolean_vectors=1),
        Section5Budgets(max_equations=1),
    ):
        result = count_common_zeros(system, budgets)
        assert result.refused and result.count is None
        assert result.boolean_vectors == 0 and result.enumerated is False
        assert result.term_table_allocated is False
        assert result.coefficients_constructed is False

    monkeypatch.undo()
    monkeypatch.setattr('counting_sat_algebra.sparse_systems.build_truncated_product', boom)
    for budgets in (
        Section5Budgets(max_coefficient_bits=1),
        Section5Budgets(max_accumulator_bits=1),
    ):
        result = count_common_zeros(system, budgets)
        assert result.refused and result.count is None and result.boolean_vectors == 0
        assert result.term_table_allocated is False
    bits = count_common_zeros(system, Section5Budgets(max_coefficient_bits=1))
    assert bits.coefficients_constructed is False
    accumulator = count_common_zeros(system, Section5Budgets(max_accumulator_bits=1))
    assert accumulator.coefficients_constructed is True
    assert 'term table not allocated' in accumulator.reason


def test_small_ell_does_not_pretend_to_be_the_count():
    system = sparse_system(4, [[[0], [1]]], k=1, c=1)
    with pytest.raises(ValueError, match='n\\+1'):
        count_common_zeros(system, ell=4)
    with pytest.raises(ValueError, match='positive'):
        count_common_zeros(system, ell=0)


def test_scope_text_keeps_theorem_1_1_open():
    assert 'does not implement Theorem 1.1' in STEP17_SCOPE
    assert '2^{o(n)}' in STEP17_SCOPE
    result = count_common_zeros(sparse_system(3, [[[0]]], k=1, c=1))
    assert result.counted and result.count == 4
    assert 'Theorem 1.1' in result.note


def test_step17_reference_rows_and_cli(tmp_path, capsys):
    rows = run_section5_suite()
    by_name = {row['name']: row for row in rows}
    assert by_name['lemma-5.1-integers']['label'] == 'CHECKED'
    assert by_name['lemma-5.1-ell-2-product']['label'] == 'CHECKED'
    assert by_name['lemma-5.1-ell-2-product']['p_coefficients'] == [1, 0, -3, 2]
    assert by_name['sparse-f2-lift']['label'] == 'CHECKED'
    counted = by_name['tiny-common-zeros-n4']
    assert counted['label'] == 'COUNTED' and counted['count'] == counted['brute'] == 4
    assert counted['boolean_vectors'] == 16 and counted['terms'] == 15
    assert by_name['empty-system-n5']['count'] == 32 and by_name['empty-system-n5']['boolean_vectors'] == 0
    assert by_name['constant-one-equation']['count'] == 0
    for name in ('term-table-budget', 'coefficient-bit-budget', 'cube-budget'):
        row = by_name[name]
        assert row['label'] == 'REFUSED' and row['count'] is None
        assert row['boolean_vectors'] == 0 and row['term_table_allocated'] is False
    assert all(row['ok'] for row in rows)
    assert {row['label'] for row in rows} == {'COUNTED', 'CHECKED', 'REFUSED'}

    output = tmp_path / 'step17.json'
    assert main(['verify', '--step17', '--output', str(output)]) == 0
    report = json.loads(output.read_text(encoding='utf-8'))
    text = capsys.readouterr().out
    assert report['version'] == __version__
    assert report['steps'] == [17] and report['all_checks_passed'] is True
    assert report['theorem_1_1_complete'] is report['lean_proof_complete'] is False
    assert report['failed_step'] is None and report['wall_times_recorded'] is False
    assert 'Theorem 1.1' in report['scope'] and 'not implemented' in report['scope']
    assert 'ALL TESTS PASSED' in text
    assert 'Theorem 1.1 is not implemented. This program does not refute ETH.' in text
    assert 'No Theorem 1.1, no 2^{o(n)} claim.' in text
    assert len(report['rows']) == len(rows)
