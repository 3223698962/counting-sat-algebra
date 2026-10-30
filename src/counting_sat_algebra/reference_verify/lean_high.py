"""Step 26 finite witnesses for the high-variable frontend and rank loop."""

from __future__ import annotations

from itertools import product
from math import isqrt

from .protocol import _row, stamp_completion
from ..cnf import CNF, clause_poly_bool
from ..lean_greedy import count_cnf_greedy
from ..lean_high import (
    STEP26_SCOPE,
    HighBudgets,
    count_cnf_high_loop,
    early_log,
    early_threshold,
    high_frontend,
)
from ..lean_ring import LeanAlgebraBudgets
from ..lean_query import TernaryHashRow, count_cnf_ternary, count_ternary_query
from ..lean_sparsify import count_cnf_streaming

HUB = CNF(8, [[(0, False), (i, False)] for i in range(1, 8)])
STAR = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
FOREST = CNF(
    6,
    [
        [(0, False), (1, False)],
        [(0, False), (2, False)],
        [(3, False), (4, False)],
        [(3, False), (5, False)],
    ],
)
TWO_HUB = CNF(
    16,
    [[(0, False), (i, False)] for i in range(1, 16)]
    + [[(1, False), (j, False)] for j in range(2, 16)],
)


def _brute_count(cnf: CNF) -> int:
    total = 0
    for bits in product((0, 1), repeat=cnf.n):
        if all(clause_poly_bool(bits, clause) for clause in cnf.clauses):
            total += 1
    return total


def _independent_frontend(cnf: CNF):
    ambient = cnf.n
    scopes = [sorted({variable for variable, _ in clause}) for clause in cnf.clauses]
    total_scope = sum(len(scope) for scope in scopes)
    threshold = isqrt((ambient + 16).bit_length()) + 1
    cutoff = (total_scope * threshold + ambient - 1) // ambient if ambient else 0
    degrees = [0] * ambient
    for scope in scopes:
        for variable in scope:
            degrees[variable] += 1
    high = tuple(v for v in range(ambient) if cutoff < degrees[v])
    return threshold, total_scope, cutoff, degrees, high


def run_lean_high_suite() -> list[dict]:
    print('=== (AB) High-variable frontend and streamed rank loop ===', flush=True)
    rows = []

    log_ok = all(
        early_log(n) == (n + 16).bit_length()
        and early_threshold(n) == isqrt((n + 16).bit_length()) + 1
        and early_threshold(n) >= 2
        for n in (0, 1, 8, 15, 16, 255, 256, 1040, 4096)
    )
    rows.append(
        _row(
            'lean-high-early-log-and-threshold',
            0,
            0,
            log_ok,
            'CHECKED',
            count=None,
            samples=[early_log(n) for n in (0, 8, 1040)],
            thresholds=[early_threshold(n) for n in (0, 8, 1040)],
        )
    )

    frontend_ok = True
    details = []
    for cnf in (HUB, STAR, FOREST, TWO_HUB):
        frontend = high_frontend(cnf)
        threshold, total_scope, cutoff, degrees, high = _independent_frontend(cnf)
        valid = (
            (frontend.threshold, frontend.total_scope, frontend.cutoff)
            == (threshold, total_scope, cutoff)
            and frontend.degrees == tuple(degrees)
            and frontend.high_variables == high
            and frontend.cardinality * (cutoff + 1) <= total_scope
            and frontend.cardinality * threshold <= frontend.ambient
            and frontend.mask() == tuple(int(v in set(high)) for v in range(frontend.ambient))
        )
        frontend_ok = frontend_ok and valid
        details.append(dict(n=cnf.n, m=cnf.m, cutoff=cutoff, high=list(high), ranks=frontend.ranks))
    rows.append(
        _row(
            'lean-high-frontend-cutoff-bounds-and-mask',
            0,
            0,
            frontend_ok,
            'CHECKED',
            count=None,
            cases=details,
        )
    )

    boundary = high_frontend(
        CNF(4, [[(0, False), (1, False)], [(0, False), (2, False)]]), threshold=1
    )
    rows.append(
        _row(
            'lean-high-frontend-boundary-excludes-equal-degree',
            4,
            2,
            boundary.cutoff == 1
            and boundary.high_variables == (0,)
            and boundary.degrees == (2, 1, 1, 0),
            'CHECKED',
            count=None,
            cutoff=boundary.cutoff,
            high=list(boundary.high_variables),
        )
    )

    hub = count_cnf_high_loop(HUB)
    rows.append(
        _row(
            'lean-high-loop-unit-rows-count-models',
            HUB.n,
            HUB.m,
            hub.counted
            and hub.count == _brute_count(HUB) == 129
            and hub.count == count_cnf_ternary(HUB).count
            and hub.ranks == 2
            and hub.queries == 2
            and hub.high_cardinality == 1
            and hub.accumulator_bits <= HUB.n + 1,
            'COUNTED',
            count=hub.count,
            ranks=hub.ranks,
            queries=hub.queries,
            high_cardinality=hub.high_cardinality,
            accumulator_bits=hub.accumulator_bits,
        )
    )

    empty_hub = count_cnf_high_loop(HUB, rows=())
    empty_star = count_cnf_high_loop(STAR, rows=())
    rows.append(
        _row(
            'lean-high-loop-empty-history-retains-every-bit',
            HUB.n,
            HUB.m,
            empty_hub.count == 1 << HUB.n == 256 and empty_star.count == 1 << STAR.n == 16,
            'CHECKED',
            count=None,
            hub=empty_hub.count,
            star=empty_star.count,
        )
    )

    history = (TernaryHashRow(((0, 1),)),)
    explicit = count_cnf_high_loop(STAR, history)
    explicit_marked = count_cnf_high_loop(STAR, history, marked_clause=0)
    rows.append(
        _row(
            'lean-high-loop-explicit-history-identities',
            STAR.n,
            STAR.m,
            explicit.count == count_ternary_query(STAR, history).count
            and explicit_marked.count == count_ternary_query(STAR, history, marked_clause=0).count,
            'CHECKED',
            count=None,
            accepted=explicit.count,
            marked=explicit_marked.count,
        )
    )

    two = count_cnf_high_loop(TWO_HUB)
    rows.append(
        _row(
            'lean-high-loop-two-high-coordinates-sum-the-ranks',
            TWO_HUB.n,
            TWO_HUB.m,
            two.counted
            and two.count == count_cnf_ternary(TWO_HUB).count
            and two.ranks == 4
            and two.queries == 4
            and two.high_cardinality == 2
            and two.accumulator_bits <= TWO_HUB.n + 1,
            'COUNTED',
            count=two.count,
            ranks=two.ranks,
            queries=two.queries,
            high_cardinality=two.high_cardinality,
        )
    )

    conditioned = count_cnf_high_loop(HUB, fixed={7: 1})
    overlap = False
    try:
        count_cnf_high_loop(HUB, fixed={0: 1})
    except ValueError:
        overlap = True
    rows.append(
        _row(
            'lean-high-loop-fixed-conditioned-totals',
            HUB.n,
            HUB.m,
            conditioned.count == count_ternary_query(HUB, None, fixed={7: 1}).count and overlap,
            'CHECKED',
            count=None,
            conditioned=conditioned.count,
        )
    )

    rank_refused = count_cnf_high_loop(HUB, budgets=HighBudgets(max_ranks=1))
    query_refused = count_cnf_high_loop(HUB, budgets=HighBudgets(max_queries=1))
    algebra_refused = count_cnf_high_loop(
        HUB, native_budgets=LeanAlgebraBudgets(max_scalar_operations=100)
    )
    rows.append(
        _row(
            'lean-high-loop-refusals-without-partial-counts',
            0,
            0,
            rank_refused.status == 'refused'
            and rank_refused.count is None
            and rank_refused.queries == 0
            and query_refused.status == 'refused'
            and query_refused.count is None
            and query_refused.queries == 0
            and algebra_refused.status == 'refused'
            and algebra_refused.count is None
            and algebra_refused.queries >= 1
            and algebra_refused.count_word_bits == 0,
            'REFUSED',
            count=None,
            rank_reason=rank_refused.reason,
            algebra_queries=algebra_refused.queries,
            algebra_reason=algebra_refused.reason,
        )
    )

    big = LeanAlgebraBudgets(max_scalar_operations=200_000_000)
    plain = count_cnf_greedy(HUB)
    routed = count_cnf_greedy(HUB, high_loop=True, native_budgets=big)
    star_plain = count_cnf_greedy(STAR)
    star_routed = count_cnf_greedy(STAR, high_loop=True)
    rows.append(
        _row(
            'lean-high-greedy-routing-identical-counts',
            HUB.n,
            HUB.m,
            plain.counted
            and routed.counted
            and plain.count == routed.count == _brute_count(HUB) == 129
            and routed.high_loop
            and routed.high_cardinality == 1
            and routed.high_cutoff == 6
            and routed.high_ranks
            == routed.marked_queries + (1 << routed.high_cardinality) * routed.accepted_queries
            and star_plain.count == star_routed.count == 9
            and star_routed.high_cardinality == 0
            and star_routed.high_ranks
            == star_routed.marked_queries + (1 << star_routed.high_cardinality),
            'COUNTED',
            count=routed.count,
            high_ranks=routed.high_ranks,
            marked_queries=routed.marked_queries,
            rounds=routed.rounds,
        )
    )

    streamed = count_cnf_streaming(
        HUB, M=2, leaf_backend='hash', high_loop=True, native_budgets=big
    )
    star_streamed = count_cnf_streaming(STAR, M=2, leaf_backend='hash', high_loop=True)
    rows.append(
        _row(
            'lean-high-streamed-routed-leaves',
            0,
            0,
            streamed.counted
            and streamed.count == 129
            and streamed.high_loop
            and streamed.leaves_counted == 2
            and star_streamed.counted
            and star_streamed.count == 9
            and star_streamed.hash_final_mass == 0,
            'COUNTED',
            count=None,
            hub=streamed.count,
            hub_leaves=streamed.leaves_counted,
            hub_high_ranks=streamed.hash_high_ranks,
            star=star_streamed.count,
        )
    )

    free = count_cnf_high_loop(CNF(8, []))
    zero = count_cnf_high_loop(CNF(0, []))
    empty_clause = count_cnf_high_loop(CNF(8, [[]]))
    rows.append(
        _row(
            'lean-high-loop-terminal-cases',
            0,
            0,
            free.counted
            and free.count == 256
            and free.ranks == 1
            and zero.counted
            and zero.count == 1
            and empty_clause.counted
            and empty_clause.count == 0
            and empty_clause.queries == 1,
            'COUNTED',
            count=None,
            free_bits=free.count,
            zero_variables=zero.count,
            empty_clause=empty_clause.count,
        )
    )

    for row in rows:
        stamp_completion(row)
        print(f'  [AB] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP26_SCOPE}', flush=True)
    return rows
