"""Step 26 independent checks for the high-variable frontend and rank loop."""

import json
from itertools import product
from math import isqrt

import pytest

from counting_sat_algebra.__main__ import main
from counting_sat_algebra.cnf import CNF, clause_poly_bool
from counting_sat_algebra.lean_greedy import count_cnf_greedy
from counting_sat_algebra.lean_high import (
    HighBudgets,
    count_cnf_high_loop,
    early_log,
    early_threshold,
    high_frontend,
)
from counting_sat_algebra.lean_query import TernaryHashRow, count_cnf_ternary, count_ternary_query
from counting_sat_algebra.lean_ring import LeanAlgebraBudgets
from counting_sat_algebra.lean_sparsify import count_cnf_streaming

HUB = CNF(8, [[(0, False), (i, False)] for i in range(1, 8)])
STAR = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
FOREST = CNF(
    6,
    [[(0, False), (1, False)], [(0, False), (2, False)],
     [(3, False), (4, False)], [(3, False), (5, False)]],
)
TWO_HUB = CNF(
    16,
    [[(0, False), (i, False)] for i in range(1, 16)]
    + [[(1, False), (j, False)] for j in range(2, 16)],
)


def _brute_models(cnf):
    total = 0
    for bits in product((0, 1), repeat=cnf.n):
        if all(clause_poly_bool(bits, clause) for clause in cnf.clauses):
            total += 1
    return total


@pytest.mark.parametrize('ambient', [0, 1, 8, 15, 16, 255, 256, 1040, 4096])
def test_early_log_and_threshold_match_the_lean_integers(ambient):
    assert early_log(ambient) == (ambient + 16).bit_length()
    assert early_threshold(ambient) == isqrt((ambient + 16).bit_length()) + 1
    assert early_threshold(ambient) >= 2


def test_frontend_cutoff_is_an_exact_ceiling_with_both_bounds():
    for cnf, expected in ((HUB, (6, (0,), 2)), (STAR, (5, (), 1)), (FOREST, (4, (), 1))):
        frontend = high_frontend(cnf)
        scopes = [sorted({variable for variable, _ in clause}) for clause in cnf.clauses]
        total_scope = sum(len(scope) for scope in scopes)
        expected_cutoff = (total_scope * early_threshold(cnf.n) + cnf.n - 1) // cnf.n
        degrees = [0] * cnf.n
        for scope in scopes:
            for variable in scope:
                degrees[variable] += 1
        assert (frontend.cutoff, frontend.high_variables, frontend.ranks) == expected
        assert frontend.cutoff == expected_cutoff
        assert frontend.degrees == tuple(degrees)
        assert frontend.cardinality * (frontend.cutoff + 1) <= frontend.total_scope
        assert frontend.cardinality * frontend.threshold <= frontend.ambient
        assert frontend.mask() == tuple(int(v in set(frontend.high_variables)) for v in range(cnf.n))


def test_frontend_boundary_excludes_a_variable_of_exactly_cutoff_degree():
    boundary = high_frontend(STAR, threshold=1)
    assert boundary.cutoff == (6 + 3) // 4 == 2
    assert boundary.degrees == (3, 1, 1, 1)
    assert boundary.high_variables == (0,)


def test_loop_accumulates_conditioned_queries_to_the_models():
    loop = count_cnf_high_loop(HUB)
    assert loop.counted and loop.count == _brute_models(HUB) == 129
    assert loop.count == count_cnf_ternary(HUB).count
    assert loop.ranks == 2 and loop.queries == 2 and loop.high_cardinality == 1
    assert loop.accumulator_bits <= HUB.n + 1
    assert loop.boolean_vectors == 0 and loop.same_declared_n


def test_loop_empty_history_retains_every_ambient_bit():
    assert count_cnf_high_loop(HUB, rows=()).count == 256
    assert count_cnf_high_loop(STAR, rows=()).count == 16
    assert count_cnf_high_loop(HUB, rows=()).count == count_ternary_query(HUB).count


def test_loop_matches_explicit_history_and_marked_counts():
    history = (TernaryHashRow(((0, 1),)),)
    assert count_cnf_high_loop(STAR, history).count == count_ternary_query(STAR, history).count
    for clause in range(HUB.m):
        assert count_cnf_high_loop(HUB, marked_clause=clause).count == count_ternary_query(
            HUB, None, marked_clause=clause
        ).count


def test_loop_with_two_high_coordinates_sums_four_ranks():
    two = count_cnf_high_loop(TWO_HUB)
    assert two.counted and two.ranks == 4 and two.queries == 4 and two.high_cardinality == 2
    assert two.count == count_cnf_ternary(TWO_HUB).count
    conditioned_total = sum(
        count_ternary_query(TWO_HUB, None, fixed={0: bit0, 1: bit1}).count
        for bit0, bit1 in product((0, 1), repeat=2)
    )
    assert two.count == conditioned_total


def test_loop_fixed_conditioning_and_high_overlap_rejection():
    conditioned = count_cnf_high_loop(HUB, fixed={7: 1})
    assert conditioned.count == count_ternary_query(HUB, None, fixed={7: 1}).count
    with pytest.raises(ValueError, match='high'):
        count_cnf_high_loop(HUB, fixed={0: 1})


def test_loop_refusals_record_work_without_partial_counts():
    rank_refused = count_cnf_high_loop(HUB, budgets=HighBudgets(max_ranks=1))
    assert rank_refused.status == 'refused' and rank_refused.count is None
    assert rank_refused.queries == 0 and rank_refused.count_word_bits == 0
    query_refused = count_cnf_high_loop(HUB, budgets=HighBudgets(max_queries=1))
    assert query_refused.status == 'refused' and query_refused.queries == 0
    algebra = count_cnf_high_loop(HUB, native_budgets=LeanAlgebraBudgets(max_scalar_operations=100))
    assert algebra.status == 'refused' and algebra.count is None
    assert algebra.queries >= 1 and algebra.scalar_operations <= 100


def test_greedy_routing_counts_identically_to_direct_queries():
    big = LeanAlgebraBudgets(max_scalar_operations=200_000_000)
    plain = count_cnf_greedy(HUB)
    routed = count_cnf_greedy(HUB, high_loop=True, native_budgets=big)
    assert plain.count == routed.count == _brute_models(HUB) == 129
    assert routed.high_loop and routed.high_cardinality == 1 and routed.high_cutoff == 6
    assert routed.high_ranks == (
        routed.marked_queries + (1 << routed.high_cardinality) * routed.accepted_queries
    )
    assert routed.marker_sums == ()
    star_plain = count_cnf_greedy(STAR)
    star_routed = count_cnf_greedy(STAR, high_loop=True)
    assert star_plain.count == star_routed.count == 9
    assert star_routed.high_cardinality == 0
    assert star_routed.high_ranks == star_routed.marked_queries + 1


def test_streamed_routed_leaves_keep_the_shared_frontier():
    big = LeanAlgebraBudgets(max_scalar_operations=200_000_000)
    streamed = count_cnf_streaming(HUB, M=2, leaf_backend='hash', high_loop=True, native_budgets=big)
    assert streamed.counted and streamed.count == 129 and streamed.high_loop
    assert streamed.leaves_counted == 2 and streamed.hash_final_mass == 0
    assert streamed.hash_high_ranks > 0
    star_streamed = count_cnf_streaming(STAR, M=2, leaf_backend='hash', high_loop=True)
    assert star_streamed.count == 9 and star_streamed.hash_final_mass == 0
    demo = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])
    demo_streamed = count_cnf_streaming(demo, M=2, leaf_backend='hash', high_loop=True)
    assert demo_streamed.count == 196


def test_loop_terminal_cases():
    free = count_cnf_high_loop(CNF(8, []))
    zero = count_cnf_high_loop(CNF(0, []))
    empty_clause = count_cnf_high_loop(CNF(8, [[]]))
    assert free.count == 256 and free.ranks == 1
    assert zero.count == 1 and zero.ranks == 1
    assert empty_clause.count == 0 and empty_clause.queries == 1


def test_count_command_high_loop_flag(tmp_path, capsys):
    demo = tmp_path / 'demo.cnf'
    demo.write_text('p cnf 8 2\n1 2 3 0\n4 5 6 0\n', encoding='utf-8')
    assert main(['count', str(demo), '--backend', 'lean', '--leaf-backend', 'hash',
                 '--high-loop', '--json']) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed['status'] == 'counted' and printed['count'] == '196'
    assert printed['high_loop'] is True and printed['hash_final_mass'] == 0
    hub = tmp_path / 'hub.cnf'
    hub.write_text(
        'p cnf 8 7\n' + '\n'.join(f'1 {i} 0' for i in range(2, 9)) + '\n', encoding='utf-8'
    )
    assert main(['count', str(hub), '--backend', 'lean', '--leaf-backend', 'hash',
                 '--high-loop', '--sparsify-M', '2', '--max-algebra-operations', '50000000',
                 '--json']) == 0
    document = json.loads(capsys.readouterr().out)
    assert document['count'] == '129' and document['hash_high_ranks'] > 0
    refused = main(['count', str(hub), '--backend', 'lean', '--leaf-backend', 'hash',
                    '--high-loop', '--json'])
    assert refused == 2
    document = json.loads(capsys.readouterr().out)
    assert document['status'] == 'refused' and document['count'] is None
