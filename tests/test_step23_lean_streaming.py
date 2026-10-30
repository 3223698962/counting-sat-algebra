"""Semantic, streaming and resource contracts transferred from SqrtLogSpace."""

import json
import random
from itertools import combinations, product

import pytest

from counting_sat_algebra import __main__, lean_sparsify
from counting_sat_algebra.cnf import CNF, brute_count, is_satisfied
from counting_sat_algebra.lean_encoding import (
    cnf_numbers, count_word, decode_cnf, decode_count_word, encode_cnf,
    formula_bit_length, literal_code,
)
from counting_sat_algebra.lean_parameters import choose_sparse_parameter, parameter_test
from counting_sat_algebra.lean_sparsify import (
    StreamingBudgets, StreamingRefusal, count_cnf_streaming, normalize_lean_cnf,
    select_priority, split_priority,
)


def _star(n=4):
    return CNF(n, [[(0, False), (var, False)] for var in range(1, n)])


def _models(cnf):
    return {bits for bits in product((0, 1), repeat=cnf.n) if is_satisfied(bits, cnf)}


def _frame(bits):
    return ''.join('1' + bit for bit in bits) + '0'


def _list(items):
    return ''.join('1' + _frame(item) for item in items) + '0'


def _oracle_encoding(cnf, canonical=True):
    rows = [[2 * var + int(not neg) for var, neg in clause] for clause in cnf.clauses]
    if canonical:
        rows = sorted({tuple(sorted(set(row))) for row in rows})
    return '1' * cnf.n + '0' + _list([_list(['1' * code + '0' for code in row]) for row in rows])


@pytest.mark.parametrize('canonical', [False, True])
def test_lean_encoding_matches_explicit_framing_and_preserves_free_bits(canonical):
    formula = CNF(5, [[(0, False), (1, True), (0, False)], [], [(0, False), (1, True)]])
    word = encode_cnf(formula, canonical=canonical)
    assert ''.join(str(bit) for bit in word) == _oracle_encoding(formula, canonical)
    assert len(word) == formula_bit_length(formula, canonical=canonical)
    decoded = decode_cnf(word)
    assert decoded.n == 5 and _models(decoded) == _models(formula)
    assert literal_code((0, True)) == 0 and literal_code((0, False)) == 1
    assert len(encode_cnf(CNF(0, []))) == 2
    assert len(encode_cnf(CNF(0, [[]]))) == 6


def test_canonical_encoding_distinguishes_syntax_and_declared_n():
    left = CNF(2, [[(1, True), (0, False), (0, False)], [(0, False)]])
    right = CNF(2, [[(0, False)], [(0, False), (1, True)], [(0, False)]])
    assert encode_cnf(left) == encode_cnf(right)
    assert encode_cnf(left, canonical=False) != encode_cnf(right, canonical=False)
    assert encode_cnf(left) != encode_cnf(CNF(3, left.clauses))
    tautology = CNF(1, [[(0, False), (0, True)]])
    assert cnf_numbers(tautology) == ((0, 1),)
    assert encode_cnf(tautology) != encode_cnf(CNF(1, []))
    with pytest.raises(ValueError, match='above budget'):
        encode_cnf(left, max_bits=0)


@pytest.mark.parametrize('word', [b'', b'\x01', b'\x02', b'\x00', b'\x00\x01', b'\x00\x00\x00'])
def test_malformed_cnf_words_are_rejected(word):
    with pytest.raises(ValueError):
        decode_cnf(word)


@pytest.mark.parametrize('count', [0, 1, 2, 3, 196, 1 << 500, (1 << 500) - 1])
def test_binary_count_word_has_lean_zero_and_little_endian_semantics(count):
    word = count_word(count)
    assert decode_count_word(word) == count
    assert len(word) == max(1, count.bit_length())
    assert count_word(0) == b'\x00' and count_word(6) == b'\x00\x01\x01'


def test_sparse_parameter_is_the_first_passing_dyadic_candidate():
    for k in range(2, 6):
        for numerator, denominator in ((1, 1), (1, 7), (3, 11)):
            parameter = choose_sparse_parameter(k, numerator, denominator)
            M = parameter.M
            depth = 4 * sum((3 * M) ** j for j in range(k - 1))
            right = 4 * sum(3 ** j for j in range(k - 1))
            assert M >= max(2, right)
            assert right * (4 * depth * M - 1).bit_length() * denominator <= numerator * M
            assert all(not parameter_test(k, numerator, denominator, e) for e in range(parameter.exponent))
            assert parameter.occurrence_bound == 2 * sum(M ** j for j in range(k))
    with pytest.raises(ValueError, match='exponent budget'):
        choose_sparse_parameter(3, max_exponent=0)


def test_normalization_preserves_models_and_does_not_propagate_units():
    formula = CNF(5, [[(0, False)], [(0, False), (1, True)], [(1, False), (1, True)], [(0, False)]])
    normalized = normalize_lean_cnf(formula)
    assert normalized.n == 5 and normalized.clauses == [[(0, False)]]
    assert _models(normalized) == _models(formula)
    contradictory = normalize_lean_cnf(CNF(3, [[(0, False)], [(0, True)]]))
    assert len(contradictory.clauses) == 2 and contradictory.n == 3


def test_priority_and_all_selected_children_partition_every_assignment():
    formula = _star()
    clauses = tuple(frozenset(clause) for clause in formula.clauses)
    choice = select_priority(clauses, 2)
    assert choice is not None and choice.width == 2 and choice.core == frozenset(((0, False),))
    assert choice.degree == 3 and choice.threshold == 2
    left, right = split_priority(clauses, choice)
    assert left == (frozenset(((0, False),)),)
    assert all(len(clause) == 1 for clause in right) and len(right) == 4
    a, b = (CNF(formula.n, [list(clause) for clause in child]) for child in (left, right))
    assert not (_models(a) & _models(b))
    assert _models(a) | _models(b) == _models(formula)


def test_priority_uses_smallest_width_then_largest_core():
    formula = tuple(frozenset(clause) for clause in [
        [(0, False), (1, False), (2, False)],
        [(0, False), (1, False), (3, False)],
        [(0, False), (4, False), (5, False)],
        [(0, False), (4, False), (6, False)],
    ])
    choice = select_priority(formula, 2)
    assert choice is not None and len(choice.core) == 2 and choice.width == 3
    narrow = formula + (frozenset(((7, False), (8, False))), frozenset(((7, False), (9, False))))
    choice = select_priority(narrow, 2)
    assert choice is not None and choice.width == 2
    with pytest.raises(StreamingRefusal):
        select_priority(narrow, 2, max_candidates=0)


@pytest.mark.parametrize('backend', ['common', 'tensor'])
def test_streaming_counts_with_algebraic_leaves_and_no_truth_table_backend(backend, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError('the existing materialized tree or cube oracle was called')

    from counting_sat_algebra import cnf, equation_selection

    monkeypatch.setattr(cnf, 'brute_count', forbidden)
    monkeypatch.setattr(equation_selection, 'brute_count', forbidden)
    monkeypatch.setattr(equation_selection, 'build_disjoint_sparsification', forbidden)
    result = count_cnf_streaming(_star(), M=2, leaf_backend=backend)
    assert result.counted and result.count == 9
    assert result.leaves == result.leaves_counted == 2
    assert result.nodes == result.planning_nodes == 3
    assert result.max_pending_nodes <= result.depth + 1
    assert result.completed_leaves_retained == 0 and result.same_declared_n
    assert result.accumulator_bits <= 5 and result.count_word_bits == 4
    assert result.boolean_vectors > 0 and result.memo_hits == 0
    assert result.theorem_1_1_complete is result.lean_proof_complete is result.recursive_contraction is False


def test_random_small_formulas_match_independent_cube_counts():
    rng = random.Random(230)
    for _ in range(14):
        clauses = []
        for _ in range(rng.randrange(1, 9)):
            variables = rng.sample(range(4), rng.choice((1, 2, 3)))
            clauses.append([(var, rng.choice((False, True))) for var in variables])
        formula = CNF(4, clauses)
        expected = brute_count(formula)
        result = count_cnf_streaming(formula, M=2, leaf_backend='tensor')
        assert result.counted and result.count == expected
        assert result.boolean_vectors <= 5_000_000


@pytest.mark.parametrize('budgets', [
    StreamingBudgets(max_nodes=1), StreamingBudgets(max_depth=0),
    StreamingBudgets(max_frontier_bits=1), StreamingBudgets(max_formula_bits=1),
    StreamingBudgets(max_search_candidates=0), StreamingBudgets(max_width=1),
    StreamingBudgets(max_clauses=2), StreamingBudgets(max_variables=3),
    StreamingBudgets(max_parameter_exponent=0),
])
def test_structural_preflight_refuses_before_any_leaf_call(budgets, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError('a leaf evaluator ran after a structural refusal')

    monkeypatch.setattr(lean_sparsify, 'count_cnf_components', forbidden)
    result = count_cnf_streaming(_star(), budgets=budgets, M=None if budgets.max_parameter_exponent == 0 else 2)
    assert result.status == 'refused' and result.count is None
    assert result.nodes == result.boolean_vectors == result.leaves_counted == 0


def test_global_leaf_budget_refuses_without_publishing_the_partial_sum():
    result = count_cnf_streaming(_star(), M=2, total_vector_budget=12)
    assert result.status == 'refused' and result.count is None
    assert result.leaves_counted == 1 and result.boolean_vectors == 10
    assert result.accumulator_bits == 4 and result.count_word_bits == 0
    assert result.completed_leaves_retained == 0


@pytest.mark.parametrize('formula,expected', [(CNF(8, []), 256), (CNF(8, [[]]), 0)])
def test_empty_cases_preserve_declared_variables_without_enumeration(formula, expected):
    result = count_cnf_streaming(formula)
    assert result.count == expected and result.boolean_vectors == 0
    assert result.count_word_bits <= formula.n + 1 and result.same_declared_n


def test_many_independent_branches_keep_only_the_active_frontier():
    formula = CNF(9, [[(3 * i, False), (3 * i + offset, False)]
                      for i in range(3) for offset in (1, 2)])
    result = count_cnf_streaming(formula, M=2, leaf_backend='tensor')
    assert result.counted and result.count == brute_count(formula)
    assert result.branches > 1 and result.max_pending_nodes <= result.depth + 1
    assert result.completed_leaves_retained == 0


def test_lean_backend_is_available_to_dimacs_and_refusals_have_no_answer(tmp_path, capsys):
    path = tmp_path / 'star.cnf'
    path.write_text('p cnf 4 3\n1 2 0\n1 3 0\n1 4 0\n', encoding='utf-8')
    args = ['count', str(path), '--backend', 'lean', '--sparsify-M', '2', '--json']
    assert __main__.main(args) == 0
    report = json.loads(capsys.readouterr().out)
    assert report['count'] == '9' and report['satisfiable'] is True
    assert report['leaves_counted'] == 2 and report['completed_leaves_retained'] == 0
    assert __main__.main(args + ['--max-sparsify-nodes', '1']) == 2
    report = json.loads(capsys.readouterr().out)
    assert report['count'] is report['satisfiable'] is report['count_encoding'] is None


def test_priority_key_matches_an_independent_finite_degree_search():
    rng = random.Random(231)
    for _ in range(20):
        formula = tuple(frozenset((v, False) for v in row) for row in
                        rng.sample(list(combinations(range(5), 3)), 6))
        viable = []
        for width in range(2, 4):
            for size in range(1, width):
                for variables in combinations(range(5), size):
                    core = frozenset((v, False) for v in variables)
                    degree = sum(len(clause) == width and core <= clause for clause in formula)
                    if degree >= 2 ** (width - size):
                        viable.append((width, -size))
        choice = select_priority(formula, 2)
        assert choice is not None and (choice.width, -len(choice.core)) == min(viable)
