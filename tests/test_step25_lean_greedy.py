"""Step 25 independent checks for sparse codes and the greedy selector."""

import json
from itertools import product

import pytest

from counting_sat_algebra.__main__ import main
from counting_sat_algebra.cnf import CNF, clause_poly_bool
from counting_sat_algebra.lean_code import (
    CodeBudgets,
    SparseTannerCode,
    canonical_hash_candidates,
    certify_sparse_code,
    construct_sparse_code,
    radix_words,
    search_local_ternary_matrix,
    search_walsh_shifts,
    sparse_code_parameters,
)
from counting_sat_algebra.lean_greedy import (
    GreedyBudgets,
    count_cnf_greedy,
    occurrence_bound,
    select_history_rows,
)
from counting_sat_algebra.lean_ring import LeanAlgebraBudgets
from counting_sat_algebra.lean_sparsify import count_cnf_streaming


def _models(cnf):
    total = 0
    for bits in product((0, 1), repeat=cnf.n):
        if all(clause_poly_bool(bits, clause) for clause in cnf.clauses):
            total += 1
    return total


def _accepted(cnf, rows, bits):
    bad = [1 - clause_poly_bool(bits, clause) for clause in cnf.clauses]
    return all(sum(weight * bad[index] for index, weight in row.terms) % 3 == 0 for row in rows)


def _brute_mass(cnf, rows):
    total = 0
    for bits in product((0, 1), repeat=cnf.n):
        if _accepted(cnf, rows, bits):
            total += sum(1 - clause_poly_bool(bits, clause) for clause in cnf.clauses)
    return total


def _brute_trajectory(cnf, initial, depth, max_rounds):
    rows = list(initial)
    trail = [_brute_mass(cnf, rows)]
    while trail[-1] > 0 and len(trail) <= max_rounds:
        best = best_score = None
        for _rank, candidate in canonical_hash_candidates(cnf.m, depth):
            score = _brute_mass(cnf, [candidate, *rows])
            if best is None or score < best_score:
                best, best_score = candidate, score
        rows = [best, *rows]
        trail.append(best_score)
    return rows, trail


CASES = ((4, 1), (16, 1), (16384, 1), (256, 2), (2, 3))


@pytest.mark.parametrize('m,B', CASES)
def test_scheduled_parameters_match_the_closed_form_integers(m, B):
    parameters = sparse_code_parameters(m, B)
    logarithm = m.bit_length()
    gap = 200 * B * logarithm
    degree = min(m, 8 * gap * gap * logarithm)
    dimension = (max(1, (m + degree - 1) // degree) - 1).bit_length()
    assert parameters.log_size == logarithm
    assert parameters.gap_denominator == gap
    assert parameters.degree_bound == 8 * gap * gap * logarithm
    assert parameters.degree == degree
    assert parameters.dimension == dimension
    assert parameters.local_rows == degree // (8 * B) + 1
    assert parameters.local_distance == degree // (100 * B * logarithm) + 1
    assert parameters.distance_denominator == 20000 * B * B * logarithm * logarithm
    assert parameters.total_rows == 2 * (1 << dimension) * parameters.local_rows


def test_local_search_returns_the_first_radix_order_witness():
    digits = [0, 0]
    tried = 0
    first = None
    while first is None:
        tried += 1
        matrix = ((digits[0], digits[1]),)
        codeword_weights = []
        for coefficients in product((0, 1, 2), repeat=2):
            if not any(coefficients):
                continue
            dot = matrix[0][0] * coefficients[0] + matrix[0][1] * coefficients[1]
            if dot % 3 == 0:
                codeword_weights.append(sum(1 for c in coefficients if c))
        if all(weight >= 2 for weight in codeword_weights):
            first = matrix
        else:
            for position in range(2):
                digits[position] += 1
                if digits[position] < 3:
                    break
                digits[position] = 0
    witness = search_local_ternary_matrix(2, 1, 2)
    assert witness.status == 'checked'
    assert witness.witness == first == ((1, 1),)
    assert witness.candidates == tried
    assert search_local_ternary_matrix(2, 1, 3).status == 'not_found'


def test_walsh_search_first_witness_matches_brute_characters():
    def good(dimension, shifts, denominator):
        return all(
            denominator
            * abs(sum(-1 if (character & shift).bit_count() % 2 else 1 for shift in shifts))
            <= len(shifts)
            for character in range(1, 1 << dimension)
        )

    candidate = search_walsh_shifts(1, 2, 2)
    assert candidate.status == 'checked' and good(1, candidate.witness, 2)
    assert candidate.witness == (1, 0) and candidate.candidates == 2
    wide = search_walsh_shifts(2, 4, 2)
    assert wide.status == 'checked' and good(2, wide.witness, 2)
    refused = search_walsh_shifts(2, 4, 2, budgets=CodeBudgets(max_candidates=1))
    assert refused.status == 'refused' and refused.witness is None


def test_rows_regenerate_from_slots_within_bounds():
    shifts = (0, 1, 2, 3)
    local = ((1, 2, 1, 0), (0, 1, 2, 1))
    code = SparseTannerCode(10, 2, shifts, local, 600, 1)
    assert code.row_count == 2 * (1 << 2) * len(local) == 16
    for index in range(code.row_count):
        chronical = code.row_count - 1 - index
        vertex, within = divmod(chronical, 2 * len(local))
        side, local_index = divmod(within, len(local))
        expected = tuple(
            (column, coefficient)
            for port, coefficient in enumerate(local[local_index])
            if (column := (vertex ^ shifts[port] if side else vertex) * len(shifts) + port) < 10
        )
        assert code.row(index).terms == expected
        assert len(expected) <= code.degree
        assert all(coefficient in (0, 1, 2) for _, coefficient in expected)
        assert code.rows[index].terms == expected


def test_certify_rejects_bad_distance_and_bad_walsh():
    with pytest.raises(ValueError, match='distance'):
        certify_sparse_code(SparseTannerCode(2, 1, (1, 0), ((1, 1),), 600, 3))
    with pytest.raises(ValueError, match='Walsh'):
        certify_sparse_code(SparseTannerCode(2, 1, (1, 1), ((1, 1),), 600, 2))
    certificate = certify_sparse_code(SparseTannerCode(2, 1, (1, 0), ((1, 1),), 600, 2))
    assert certificate.minimum_weight == 2 and certificate.depth == 1


def test_construct_refuses_scheduled_and_large_instances_without_a_code():
    scheduled = construct_sparse_code(16384, 1)
    assert scheduled.status == 'refused' and scheduled.certificate is None
    assert scheduled.candidates == 0 and 'no partial sparse code' in scheduled.reason
    large = construct_sparse_code(1 << 20, 1)
    assert large.status == 'refused' and large.certificate is None
    assert large.candidates == 0
    capped = construct_sparse_code(64, 1, budgets=CodeBudgets(max_candidates=1))
    assert capped.status == 'refused' and capped.certificate is None
    with pytest.raises(ValueError):
        construct_sparse_code(0, 1)


def test_candidate_stream_is_one_increasing_rank_per_ternary_tuple():
    stream = list(canonical_hash_candidates(2, 1))
    assert [rank for rank, _ in stream] == list(range(6))
    assert [row.terms for _, row in stream] == [
        ((0, 0),), ((1, 0),), ((0, 1),), ((1, 1),), ((0, 2),), ((1, 2),)
    ]
    deeper = list(canonical_hash_candidates(3, 2))
    assert len(deeper) == 81
    ranks = [rank for rank, _ in deeper]
    assert ranks == sorted(ranks) and len(set(ranks)) == 81
    assert {row.terms for _, row in deeper} == {
        tuple((digit % 3, digit // 3) for digit in digits)
        for digits in product(range(9), repeat=2)
    }
    deep_width = (3 * 3).bit_length()
    assert all(
        rank == sum(
            (index + 3 * weight) << (deep_width * position)
            for position, (index, weight) in enumerate(row.terms)
        )
        for rank, row in deeper
    )
    assert list(radix_words(3, 2))[:5] == [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1)]


def test_greedy_star_reproduces_brute_trajectory_and_counts_models():
    star = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
    result = count_cnf_greedy(star, depth=1)
    assert result.counted and result.count == _models(star) == 9
    assert occurrence_bound(star) == 3
    initial = result.rows[result.rounds:]
    rows, trail = _brute_trajectory(star, initial, 1, result.max_rounds)
    assert tuple(rows) == result.rows
    assert trail[0] == result.initial_mass == 12
    assert len(trail) - 1 == result.rounds
    assert result.final_mass == 0 and trail[-1] == 0
    assert result.accepted_queries == 1
    assert _brute_mass(star, result.rows) == 0


def test_greedy_zero_mass_accepts_exactly_the_models():
    star = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
    selected = select_history_rows(star, depth=1)
    assert selected.selected and selected.final_mass == 0
    accepted = {bits for bits in product((0, 1), repeat=4) if _accepted(star, selected.rows, bits)}
    models = {
        bits for bits in product((0, 1), repeat=4)
        if all(clause_poly_bool(bits, clause) for clause in star.clauses)
    }
    assert accepted == models and len(models) == 9


def test_greedy_refusals_record_work_without_partial_counts():
    star = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
    forest = CNF(
        6,
        [[(0, False), (1, False)], [(0, False), (2, False)],
         [(3, False), (4, False)], [(3, False), (5, False)]],
    )
    capped = count_cnf_greedy(star, depth=1, budgets=GreedyBudgets(max_candidates=1))
    assert capped.status == 'refused' and capped.count is None
    assert capped.marked_queries == 0 and capped.candidates_scanned == 0
    assert capped.code_status == '' and capped.count_word_bits == 0
    rounds = count_cnf_greedy(forest, depth=1, budgets=GreedyBudgets(max_rounds=1))
    assert rounds.status == 'refused' and rounds.count is None
    assert rounds.rounds == 1 and rounds.final_mass > 0 and rounds.marked_queries == 52
    algebra = count_cnf_greedy(
        star, depth=1, native_budgets=LeanAlgebraBudgets(max_scalar_operations=2000)
    )
    assert algebra.status == 'refused' and algebra.count is None
    assert algebra.marked_queries > 0 and algebra.scalar_operations <= 2000


def test_greedy_round_cap_defaults_to_the_lean_loop_cap():
    star = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
    result = select_history_rows(star, depth=1)
    assert result.max_rounds == 2 * (star.n + star.m.bit_length() + 1) == 14


def test_greedy_terminal_cases_without_scans():
    free = count_cnf_greedy(CNF(8, []))
    assert free.count == 256 and free.rounds == 0 and free.marked_queries == 0
    empty = count_cnf_greedy(CNF(8, [[]]))
    assert empty.counted and empty.count == 0 and empty.accepted_queries == 1
    single_cnf = CNF(3, [[(0, False)]])
    single = count_cnf_greedy(single_cnf)
    assert single.count == 4 == _models(single_cnf)
    with pytest.raises(ValueError):
        GreedyBudgets(max_rounds=-1)
    with pytest.raises(ValueError):
        count_cnf_greedy(single_cnf, budgets=GreedyBudgets(depth=0))
    with pytest.raises(ValueError):
        count_cnf_greedy(single_cnf, depth=0)


def test_streamed_hash_leaves_count_and_share_budgets():
    star = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
    demo = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])
    forest = CNF(
        6,
        [[(0, False), (1, False)], [(0, False), (2, False)],
         [(3, False), (4, False)], [(3, False), (5, False)]],
    )
    for cnf, expected in ((star, 9), (demo, 196), (forest, 25)):
        result = count_cnf_streaming(cnf, M=2, leaf_backend='hash')
        assert result.counted and result.count == expected == _models(cnf)
        assert result.hash_final_mass == 0
        assert result.hash_accepted_queries == result.leaves_counted
        assert result.boolean_vectors == 0 and result.same_declared_n
    shared = count_cnf_streaming(
        star, M=2, leaf_backend='hash', greedy_budgets=GreedyBudgets(max_marked_queries=50)
    )
    assert shared.status == 'refused' and shared.count is None
    assert shared.hash_marked_queries == 50 and shared.leaves_counted == 1
    assert shared.count_word_bits == 0


def test_count_command_hash_leaf_reports_and_refuses(tmp_path, capsys):
    source = tmp_path / 'demo.cnf'
    source.write_text('p cnf 8 2\n1 2 3 0\n4 5 6 0\n', encoding='utf-8')
    assert main(['count', str(source), '--backend', 'lean', '--leaf-backend', 'hash', '--json']) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed['status'] == 'counted' and printed['count'] == '196'
    assert 'greedy hash leaves' in printed['algorithm']
    assert printed['leaf_backend'] == 'hash' and printed['hash_final_mass'] == 0
    assert printed['hash_accepted_queries'] == printed['leaves_counted']
    refused = main([
        'count', str(source), '--backend', 'lean', '--leaf-backend', 'hash',
        '--max-hash-candidates', '1', '--json',
    ])
    assert refused == 2
    document = json.loads(capsys.readouterr().out)
    assert document['status'] == 'refused' and document['count'] is None
