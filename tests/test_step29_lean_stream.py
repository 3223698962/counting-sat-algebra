"""Step 29 independent checks for the streamed pointwise correction-pair leaf."""

import json
import random
from itertools import product

import pytest

from counting_sat_algebra.__main__ import main
from counting_sat_algebra.cnf import CNF, brute_count, clause_poly_bool
from counting_sat_algebra.lean_query import TernaryHashRow, count_cnf_ternary, count_ternary_query
from counting_sat_algebra.lean_ring import LeanAlgebraBudgets
from counting_sat_algebra.lean_sparsify import count_cnf_streaming
from counting_sat_algebra.lean_stream import count_cnf_streamed, streamed_query_count
from counting_sat_algebra.lean_taylor import TaylorBudgets

FIXTURE = CNF(5, [[(v, False)] for v in range(3)])
HISTORY = (TernaryHashRow(((0, 1), (1, 1), (2, 1))),)
STAR = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
HUB = CNF(8, [[(0, False), (i, False)] for i in range(1, 8)])


def _random_pair(seed):
    rng = random.Random(seed)
    n = rng.randrange(3, 7)
    clauses = []
    for _ in range(rng.randrange(1, 5)):
        width = rng.randrange(1, 4)
        clause = []
        for _ in range(width):
            variable = rng.randrange(n)
            clause.append((variable, bool(rng.randrange(2))))
        clauses.append(clause)
    rows = []
    for _ in range(rng.randrange(0, 3)):
        row = []
        for _ in range(rng.randrange(1, 4)):
            row.append((rng.randrange(len(clauses)), rng.randrange(3)))
        rows.append(TernaryHashRow(tuple(row)))
    mark = rng.randrange(len(clauses)) if clauses and rng.randrange(2) else None
    fixed = {}
    if n > 2 and rng.randrange(2):
        fixed = {rng.randrange(n): rng.randrange(2)}
    return CNF(n, clauses), tuple(rows), mark, fixed


@pytest.mark.parametrize('seed', range(291, 309))
def test_streamed_queries_match_the_native_algebra(seed):
    cnf, rows, mark, fixed = _random_pair(seed)
    native = count_ternary_query(cnf, rows, marked_clause=mark, fixed=fixed)
    stream = streamed_query_count(cnf, rows, marked_clause=mark, fixed=fixed)
    assert stream.counted and native.counted
    assert stream.count == native.count
    assert stream.marker_sums == native.marker_sums
    assert stream.precision == native.precision and stream.modulus == native.modulus
    assert stream.boolean_vectors == 0 and stream.polynomial_terms == 0
    assert stream.pairs == stream.point_evaluations
    if rows or mark is not None:
        assert stream.pairs > 0


def test_streamed_counts_match_brute_on_small_cnfs():
    for cnf, expected in ((STAR, 9), (HUB, 129)):
        stream = count_cnf_streamed(cnf)
        assert stream.count == count_cnf_ternary(cnf).count == brute_count(cnf) == expected


def test_streamed_degree_bound_and_order_formula():
    stream = streamed_query_count(FIXTURE, HISTORY)
    assert stream.degree_bound == 1
    assert stream.order == 1 * (3 * (stream.precision - 1) + 2) + 1
    assert stream.precision == 4 and stream.modulus == 81
    hub_fixed = streamed_query_count(HUB, None, fixed={7: 1})
    assert hub_fixed.degree_bound == 12


def test_streamed_path_needs_no_support_table():
    tight = LeanAlgebraBudgets(max_polynomial_terms=1)
    native = count_cnf_ternary(STAR, budgets=tight)
    stream = count_cnf_streamed(STAR, algebra_budgets=tight)
    assert native.status == 'refused' and native.count is None
    assert stream.counted and stream.count == 9 and stream.polynomial_terms == 0


def test_streamed_budget_refusals_without_partial_counts():
    refused = streamed_query_count(STAR, None, budgets=TaylorBudgets(max_pairs=40))
    assert refused.status == 'refused' and refused.count is None and refused.pairs == 0
    refused = streamed_query_count(
        STAR, None, budgets=TaylorBudgets(max_pairs=1 << 20, max_assignments=40)
    )
    assert refused.status == 'refused' and refused.count is None
    refused = streamed_query_count(
        STAR, None,
        budgets=TaylorBudgets(max_pairs=1 << 30, max_assignments=1 << 30),
        algebra_budgets=LeanAlgebraBudgets(max_precision=1),
    )
    assert refused.status == 'refused' and refused.count_word_bits == 0


def test_streamed_precision_rejects_a_too_small_explicit_modulus():
    with pytest.raises(ValueError, match='exceed'):
        streamed_query_count(STAR, None, precision=1)


def test_streamed_terminal_cases():
    assert count_cnf_streamed(CNF(8, [])).count == 256
    assert count_cnf_streamed(CNF(8, [[]])).count == 0
    tautology = count_cnf_streamed(CNF(6, [[(0, False), (0, True)]]))
    assert tautology.count == 64 and tautology.n == 6
    empty_history_marked = streamed_query_count(FIXTURE, (), marked_clause=0)
    assert empty_history_marked.count == sum(
        1 for bits in product((0, 1), repeat=5)
        if not clause_poly_bool(bits, FIXTURE.clauses[0])
    ) == 16


def test_streamed_leaf_backend_matches_native():
    forest = CNF(6, [[(0, False), (1, False)], [(0, False), (2, False)],
                     [(3, False), (4, False)], [(3, False), (5, False)]])
    for cnf, expected in ((STAR, 9), (forest, 25)):
        stream = count_cnf_streaming(cnf, M=2, leaf_backend='stream')
        native = count_cnf_streaming(cnf, M=2, leaf_backend='native')
        assert stream.counted and stream.count == native.count == expected
        assert stream.stream_pairs > 0 and stream.stream_point_evaluations == stream.stream_pairs
    shared = count_cnf_streaming(
        STAR, M=2, leaf_backend='stream', taylor_budgets=TaylorBudgets(max_pairs=40)
    )
    assert shared.status == 'refused' and shared.count is None and shared.stream_pairs <= 40


def test_count_command_stream_leaf(tmp_path, capsys):
    source = tmp_path / 'demo.cnf'
    source.write_text('p cnf 8 2\n1 2 3 0\n4 5 6 0\n', encoding='utf-8')
    assert main(['count', str(source), '--backend', 'lean', '--leaf-backend', 'stream', '--json']) == 0
    document = json.loads(capsys.readouterr().out)
    assert document['status'] == 'counted' and document['count'] == '196'
    assert 'pointwise correction-pair leaves' in document['algorithm']
    assert document['stream_pairs'] > 0 and document['leaf_backend'] == 'stream'
    refused = main(['count', str(source), '--backend', 'lean', '--leaf-backend', 'stream',
                    '--max-stream-pairs', '10', '--json'])
    assert refused == 2
    document = json.loads(capsys.readouterr().out)
    assert document['status'] == 'refused' and document['count'] is None
