"""Independent algebra, query semantics and global resource contracts."""

import json
import random
from dataclasses import replace
from itertools import product
from math import comb

import pytest

from counting_sat_algebra import __main__, lean_sparsify
from counting_sat_algebra.cnf import CNF, brute_count, clause_poly_bool
from counting_sat_algebra.lean_query import TernaryHashRow, count_cnf_ternary, count_ternary_query
from counting_sat_algebra.lean_reference import DECLARATIONS, audit_lean_reference, lean_code
from counting_sat_algebra.lean_ring import (
    AlgebraMeter,
    LeanAlgebraBudgets,
    LeanAlgebraRefusal,
    TruncatedRing,
    construct_ternary_indicator,
)
from counting_sat_algebra.lean_squarefree import SquarefreeAlgebra


def _dense_product(left, right, t, J, modulus):
    result = {}
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            u, z = i // J + j // J, i % J + j % J
            if u < t and z < J:
                result[u, z] = result.get((u, z), 0) + a * b
    return tuple(result.get((u, z), 0) % modulus for u in range(t) for z in range(J))


@pytest.mark.parametrize('t,J', [(1, 1), (2, 3), (3, 4), (0, 3), (3, 0)])
def test_dense_ring_matches_independent_bivariate_convolution(t, J):
    rng = random.Random(240 + t + J)
    ring = TruncatedRing(t, J, 81)
    for _ in range(5):
        a, b, c = (
            ring.element(tuple(rng.randrange(-100, 100) for _ in range(t * J))) for _ in range(3)
        )
        expected = () if not t * J else _dense_product(a.coefficients, b.coefficients, t, J, 81)
        assert (a * b).coefficients == expected
        assert a * (b * c) == (a * b) * c
        assert a * (b + c) == a * b + a * c
        assert a - a == ring.scalar(0)
        assert a**0 == ring.scalar(1)


def test_marker_readout_must_follow_the_complete_truncated_product():
    ring = TruncatedRing(2, 3, 27)
    U, Z = ring.monomial(1, 0), ring.monomial(0, 1)
    assert (U * U).readout() == (0, 0, 0)
    assert U.readout() == (1, 0, 0)
    assert (U * Z).readout() == (0, 1, 0)
    assert (U + Z).readout() == (1, 1, 0)
    assert (Z**3).readout() == (0, 0, 0)
    with pytest.raises(ValueError, match='same truncated ring'):
        _ = U * TruncatedRing(3, 3, 27).scalar(1)


def test_dense_power_budget_covers_all_binary_exponentiation_products():
    ring = TruncatedRing(1, 1, 27, LeanAlgebraBudgets(max_scalar_operations=2))
    assert (ring.scalar(2) ** 1).coefficients == (2,)
    with pytest.raises(LeanAlgebraRefusal, match='power operation'):
        _ = ring.scalar(2) ** 4


@pytest.mark.parametrize(
    'budgets',
    [
        LeanAlgebraBudgets(max_precision=1),
        LeanAlgebraBudgets(max_taylor_length=1),
        LeanAlgebraBudgets(max_dense_slots=1),
        LeanAlgebraBudgets(max_live_coefficients=1),
        LeanAlgebraBudgets(max_algebra_bits=1),
    ],
)
def test_dense_dimensions_and_storage_refuse_before_coefficient_allocation(budgets):
    with pytest.raises(LeanAlgebraRefusal):
        TruncatedRing(2, 2, 81, budgets)


@pytest.mark.parametrize('t', range(1, 9))
def test_constructed_indicator_has_the_independent_hermite_characterization(t):
    indicator = construct_ternary_indicator(t)
    # Degree <3*t and these 3*t Hermite conditions uniquely determine the
    # zero-residue idempotent at the three pairwise unit-separated roots.
    assert len(indicator.polynomial) == 3 * t
    for root in range(3):
        for derivative in range(t):
            hasse = (
                sum(
                    comb(j, derivative) * coefficient * root ** (j - derivative)
                    for j, coefficient in enumerate(indicator.polynomial)
                    if j >= derivative
                )
                % indicator.modulus
            )
            assert hasse == int(root == 0 and derivative == 0)
    samples = range(-81, 82) if t > 5 else range(-(3**t), 3**t)
    assert all(indicator.evaluate(value) == int(value % 3 == 0) for value in samples)
    assert all(len(row) == 3 for row in indicator.digits)
    assert indicator.newton_steps == (t - 1).bit_length()


def test_ternary_coefficient_words_match_lean_fixed_width_order():
    indicator = construct_ternary_indicator(2)
    assert indicator.digits == ((1, 3, 5), (3, 7, 3))
    encoded, width = indicator.encode(), indicator.coefficient_word_bits
    assert len(encoded) == 3 * indicator.t * width == 24
    decoded = tuple(
        sum(encoded[start + bit] << bit for bit in range(width))
        for start in range(0, len(encoded), width)
    )
    assert decoded == tuple(c for row in indicator.digits for c in row)


def test_squarefree_multiplication_and_sum_match_boolean_semantics():
    algebra = SquarefreeAlgebra(4, 81, AlgebraMeter(LeanAlgebraBudgets()))
    X = algebra.variable(1)
    assert (X * X).terms == X.terms
    assert not (X * (algebra.one - X)).terms
    a = algebra.polynomial({frozenset(): -3, frozenset((0, 2)): 4, frozenset((1,)): 2})
    b = algebra.polynomial({frozenset((2,)): 7, frozenset((0, 3)): -5})
    multiplied = a * b
    for bits in product((0, 1), repeat=4):
        assert multiplied.evaluate(bits) == a.evaluate(bits) * b.evaluate(bits) % 81
    expected = sum(multiplied.evaluate(bits) for bits in product((0, 1), repeat=4)) % 81
    assert algebra.sum_boolean(multiplied) == expected
    with pytest.raises(TypeError):
        multiplied.terms[frozenset()] = 1


def _query_oracle(cnf, rows, mark=None, fixed=None):
    total = 0
    for bits in product((0, 1), repeat=cnf.n):
        if any(bits[v] != bit for v, bit in (fixed or {}).items()):
            continue
        violations = [1 - clause_poly_bool(bits, clause) for clause in cnf.clauses]
        accepted = all(
            sum(weight * violations[index] for index, weight in row.terms) % 3 == 0 for row in rows
        )
        total += accepted and (mark is None or violations[mark])
    return total


def test_shared_marker_high_coefficients_and_marked_violation_semantics():
    formula = CNF(5, [[(v, False)] for v in range(3)])
    rows = (TernaryHashRow(((0, 1), (1, 1), (2, 1))),)
    accepted = count_ternary_query(formula, rows)
    marked = count_ternary_query(formula, rows, marked_clause=0)
    assert accepted.count == _query_oracle(formula, rows) == 8
    assert marked.count == _query_oracle(formula, rows, 0) == 4
    assert accepted.count != brute_count(formula)
    assert accepted.marker_sums == (8, 36, 45, 0)
    assert marked.marker_sums == (4, 36, 45, 0)
    assert sum(accepted.marker_sums) % 81 == accepted.count
    assert accepted.same_declared_n and accepted.precision == 4
    assert accepted.boolean_vectors == marked.boolean_vectors == 0


@pytest.mark.parametrize('fixed', [{0: 1}, {0: 0, 4: 1}, {3: 0}, {0: 0, 1: 1, 2: 1}])
def test_high_conditioning_retains_original_source_ids_and_declared_n(fixed):
    formula = CNF(5, [[(v, False)] for v in range(3)])
    rows = (TernaryHashRow(((0, 1), (1, 1), (2, 1))),)
    for mark in (None, 0, 2):
        result = count_ternary_query(formula, rows, marked_clause=mark, fixed=fixed)
        assert result.count == _query_oracle(formula, rows, mark, fixed)
        assert result.n == 5 and result.fixed_variables == tuple(sorted(fixed))
        assert result.precision == 4


def test_random_history_and_marked_queries_match_an_independent_assignment_oracle():
    rng = random.Random(241)
    for _ in range(12):
        clauses = [
            [(rng.randrange(4), bool(rng.randrange(2))) for _ in range(rng.randrange(1, 4))]
            for _ in range(5)
        ]
        formula = CNF(4, clauses)
        rows = tuple(
            TernaryHashRow(tuple((rng.randrange(5), rng.randrange(3)) for _ in range(3)))
            for _ in range(2)
        )
        for mark in (None, 0, 4):
            result = count_ternary_query(formula, rows, marked_clause=mark)
            assert result.counted and result.count == _query_oracle(formula, rows, mark)
            assert result.boolean_vectors == 0 and result.count_word_bits <= 5


def test_source_repetitions_tautologies_empty_clauses_and_zero_rows():
    formula = CNF(4, [[(0, False), (0, False)], [(1, False), (1, True)], []])
    rows = (TernaryHashRow(((0, 1), (0, 2))), TernaryHashRow(()))
    assert rows[0].terms == ((0, 1), (0, 2))
    assert count_ternary_query(formula, rows).count == 16
    assert count_ternary_query(formula, rows, marked_clause=0).count == 8
    assert count_ternary_query(formula, rows, marked_clause=1).count == 0
    assert count_ternary_query(formula, rows, marked_clause=2).count == 16
    assert count_cnf_ternary(formula).count == 0
    assert count_cnf_ternary(CNF(8, [])).count == 256
    assert count_ternary_query(CNF(8, []), fixed={2: 1, 6: 0}).count == 64


@pytest.mark.parametrize(
    'budgets',
    [
        LeanAlgebraBudgets(max_variables=3),
        LeanAlgebraBudgets(max_clauses=1),
        LeanAlgebraBudgets(max_literals=1),
        LeanAlgebraBudgets(max_formula_bits=1),
        LeanAlgebraBudgets(max_query_rows=1),
        LeanAlgebraBudgets(max_row_terms=1),
        LeanAlgebraBudgets(max_precision=1),
        LeanAlgebraBudgets(max_dense_slots=1),
        LeanAlgebraBudgets(max_live_coefficients=1),
        LeanAlgebraBudgets(max_algebra_bits=1),
        LeanAlgebraBudgets(max_scalar_operations=0),
    ],
)
def test_query_preflight_refusals_do_not_construct_a_source_polynomial(budgets):
    formula = CNF(4, [[(0, False), (1, False)], [(2, True)]])
    result = count_cnf_ternary(formula, budgets=budgets)
    assert result.status == 'refused' and result.count is None
    assert result.source_polynomials == result.monomial_products == result.boolean_vectors == 0
    assert not result.indicator_constructed and result.count_word_bits == 0


@pytest.mark.parametrize(
    'budgets',
    [
        LeanAlgebraBudgets(max_monomial_products=0),
        LeanAlgebraBudgets(max_polynomial_terms=1),
    ],
)
def test_late_algebra_refusal_publishes_no_residual_or_partial_count(budgets):
    formula = CNF(4, [[(0, False), (1, False)], [(2, True)]])
    result = count_cnf_ternary(formula, budgets=budgets)
    assert result.status == 'refused' and result.count is None
    assert result.indicator_constructed and result.scalar_operations > 0
    assert result.count_word_bits == 0 and result.marker_sums == ()


def test_malformed_query_parameters_do_not_change_the_requested_semantics():
    formula = CNF(4, [[(0, False)]])
    with pytest.raises(ValueError, match='weight'):
        TernaryHashRow(((0, 3),))
    with pytest.raises(ValueError, match='clause index'):
        count_ternary_query(formula, (TernaryHashRow(((2, 1),)),))
    with pytest.raises(ValueError, match='marked clause'):
        count_ternary_query(formula, marked_clause=4)
    with pytest.raises(ValueError, match=r'3\^precision'):
        count_cnf_ternary(formula, precision=2)
    with pytest.raises(ValueError, match='Boolean values'):
        count_ternary_query(formula, fixed={0: 2})
    with pytest.raises(ValueError, match='history rows'):
        count_ternary_query(formula, ([],))
    with pytest.raises(ValueError, match='nonnegative'):
        LeanAlgebraBudgets(max_precision=-1)


def _star():
    return CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])


def test_native_frontier_calls_no_existing_counter_or_truth_table_backend(monkeypatch):
    expected = brute_count(_star())

    def forbidden(*_args, **_kwargs):
        raise AssertionError('a truth table or previous algebraic leaf was called')

    from counting_sat_algebra import cnf, equation_selection

    monkeypatch.setattr(cnf, 'brute_count', forbidden)
    monkeypatch.setattr(equation_selection, 'brute_count', forbidden)
    monkeypatch.setattr(equation_selection, 'build_disjoint_sparsification', forbidden)
    monkeypatch.setattr(lean_sparsify, 'count_cnf_components', forbidden)
    result = lean_sparsify.count_cnf_streaming(
        _star(), M=2, leaf_backend='native', total_vector_budget=0
    )
    assert result.count == expected == 9 and result.native_queries == result.leaves_counted == 2
    assert result.boolean_vectors == result.memo_hits == result.completed_leaves_retained == 0
    assert result.monomial_products > 0 and result.scalar_operations > 0
    assert result.max_pending_nodes <= result.depth + 1
    assert (
        result.theorem_1_1_complete
        is result.lean_proof_complete
        is result.recursive_contraction
        is False
    )


def test_native_operation_budget_is_global_across_the_counting_frontier():
    result = lean_sparsify.count_cnf_streaming(
        _star(),
        M=2,
        leaf_backend='native',
        native_budgets=LeanAlgebraBudgets(max_scalar_operations=1000),
    )
    assert result.status == 'refused' and result.count is None
    assert result.native_queries == 2 and result.leaves_counted == 1
    assert 0 < result.scalar_operations <= 1000
    assert result.completed_leaves_retained == result.count_word_bits == 0


def test_many_native_branches_accumulate_counts_without_saving_completed_queries():
    formula = CNF(9, [[(3 * i, False), (3 * i + j, False)] for i in range(3) for j in (1, 2)])
    result = lean_sparsify.count_cnf_streaming(formula, M=2, leaf_backend='native')
    assert result.count == brute_count(formula) == 125
    assert result.native_queries == result.leaves == 8
    assert result.max_pending_nodes == 4 and result.depth == 3
    assert result.completed_leaves_retained == result.boolean_vectors == 0


def test_native_backend_and_refusals_are_available_from_dimacs(tmp_path, capsys):
    path = tmp_path / 'star.cnf'
    path.write_text('p cnf 4 3\n1 2 0\n1 3 0\n1 4 0\n', encoding='utf-8')
    args = [
        'count',
        str(path),
        '--backend',
        'lean',
        '--leaf-backend',
        'native',
        '--sparsify-M',
        '2',
        '--total-vector-budget',
        '0',
        '--json',
    ]
    assert __main__.main(args) == 0
    report = json.loads(capsys.readouterr().out)
    assert report['count'] == '9' and report['satisfiable'] is True
    assert report['native_queries'] == 2 and report['boolean_vectors'] == 0
    assert __main__.main(args + ['--max-algebra-operations', '1000']) == 2
    report = json.loads(capsys.readouterr().out)
    assert report['count'] is report['satisfiable'] is report['count_encoding'] is None


def test_lean_lexer_distinguishes_code_from_nested_comments_and_escaped_strings():
    text = (
        '/- sorry /- axiom hidden -/ admit -/\n-- sorry\n'
        + 'def s := "sorry \\" axiom"\n'
        + 'theorem actual : True := by sorry\n'
    )
    code = lean_code(text)
    assert code.count('sorry') == 1 and 'axiom' not in code and 'admit' not in code
    assert code.count('\n') == text.count('\n')


def _fake_lean(root):
    main = root / 'OAI' / 'Algorithms' / 'SubexponentialSAT' / 'Main.lean'
    main.parent.mkdir(parents=True)
    imports = sorted({file[:-5].replace('/', '.') for file, _ in DECLARATIONS})
    main.write_text(
        'import Mathlib.Data.Nat.Basic\n' + ''.join(f'import {name}\n' for name in imports),
        encoding='utf-8',
    )
    for file in {file for file, _ in DECLARATIONS}:
        target = root / file
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            '/- sorry /- nested admit -/ -/\n'
            + ''.join(
                f'theorem {name} : True := by trivial\n'
                for path, name in DECLARATIONS
                if path == file
            ),
            encoding='utf-8',
        )
    (root / 'lean-toolchain').write_text('leanprover/lean4:v4.34.1\n', encoding='utf-8')
    return main


def test_reference_source_audit_hashes_imports_but_never_claims_kernel_verification(tmp_path):
    main = _fake_lean(tmp_path)
    report = audit_lean_reference(tmp_path)
    assert report['static_checks_passed'] and report['local_modules'] == 3
    assert report['proof_placeholders'] == report['missing_local_imports'] == []
    assert report['kernel_verified'] is report['python_algorithm_complete'] is False
    original = report['source_sha256']
    main.write_text(main.read_text(encoding='utf-8') + '-- a source change\n', encoding='utf-8')
    assert audit_lean_reference(tmp_path)['source_sha256'] != original


@pytest.mark.parametrize(
    'suffix',
    ['import OAI.Missing\n', 'theorem hole : True := by sorry\n', 'axiom unchecked : True\n'],
)
def test_reference_audit_refuses_missing_imports_and_code_level_placeholders(tmp_path, suffix):
    main = _fake_lean(tmp_path)
    main.write_text(main.read_text(encoding='utf-8') + suffix, encoding='utf-8')
    report = audit_lean_reference(tmp_path)
    assert not report['static_checks_passed'] and report['kernel_verified'] is False


def test_native_count_respects_explicit_resource_receipt_caps():
    limits = LeanAlgebraBudgets(max_algebra_bits=10_000, max_live_coefficients=1000)
    result = count_cnf_ternary(_star(), budgets=limits)
    assert result.count == 9
    assert result.peak_algebra_bits <= limits.max_algebra_bits
    assert result.peak_live_coefficients <= limits.max_live_coefficients
    refused = count_cnf_ternary(_star(), budgets=replace(limits, max_live_coefficients=1))
    assert refused.status == 'refused' and refused.count is None
