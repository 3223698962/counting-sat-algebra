"""Step 29 finite witnesses for the streamed pointwise correction-pair leaf."""

from __future__ import annotations

from .protocol import _row, stamp_completion
from ..cnf import CNF, brute_count
from ..lean_query import TernaryHashRow, count_cnf_ternary, count_ternary_query
from ..lean_ring import LeanAlgebraBudgets
from ..lean_sparsify import count_cnf_streaming
from ..lean_stream import (
    STEP29_SCOPE,
    count_cnf_streamed,
    streamed_query_count,
)
from ..lean_taylor import TaylorBudgets

FIXTURE = CNF(5, [[(v, False)] for v in range(3)])
HISTORY = (TernaryHashRow(((0, 1), (1, 1), (2, 1))),)
STAR = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
HUB = CNF(8, [[(0, False), (i, False)] for i in range(1, 8)])
FOREST = CNF(
    6,
    [
        [(0, False), (1, False)],
        [(0, False), (2, False)],
        [(3, False), (4, False)],
        [(3, False), (5, False)],
    ],
)
DEMO = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])


def run_lean_stream_suite() -> list[dict]:
    print('=== (AE) Streamed pointwise correction-pair leaf ===', flush=True)
    rows = []

    history_stream = streamed_query_count(FIXTURE, HISTORY)
    degree_ok = (
        history_stream.order == 1 * (3 * (history_stream.precision - 1) + 2) + 1
        and history_stream.degree_bound == 1
        and history_stream.pairs == 3**FIXTURE.n
        and history_stream.point_evaluations == history_stream.pairs
    )
    hub_fixed = streamed_query_count(HUB, None, fixed={7: 1})
    degree_ok = degree_ok and hub_fixed.degree_bound == 12  # clause (0 or 7) is forced violated
    rows.append(
        _row(
            'lean-stream-degree-bound-and-order',
            0,
            0,
            degree_ok,
            'CHECKED',
            count=None,
            fixture_order=history_stream.order,
            hub_fixed_degree=hub_fixed.degree_bound,
        )
    )

    sums_ok = True
    details = []
    for cnf, history, mark, fixed in (
        (FIXTURE, HISTORY, None, None),
        (FIXTURE, HISTORY, 0, None),
        (HUB, None, 0, None),
        (HUB, None, 0, {7: 1}),
        (HUB, None, None, {7: 1}),
    ):
        native = count_ternary_query(cnf, history, marked_clause=mark, fixed=fixed)
        stream = streamed_query_count(cnf, history, marked_clause=mark, fixed=fixed)
        valid = (
            stream.counted
            and stream.count == native.count
            and stream.marker_sums == native.marker_sums
            and stream.precision == native.precision
            and stream.modulus == native.modulus
            and stream.fixed_variables == native.fixed_variables
        )
        sums_ok = sums_ok and valid
        details.append(
            dict(n=cnf.n, marked=mark, fixed=tuple(sorted((fixed or {}))), count=stream.count)
        )
    rows.append(
        _row(
            'lean-stream-marker-sums-match-native',
            0,
            0,
            sums_ok,
            'CHECKED',
            count=None,
            cases=details,
        )
    )

    accepted = streamed_query_count(FIXTURE, HISTORY)
    marked = streamed_query_count(FIXTURE, HISTORY, marked_clause=0)
    rows.append(
        _row(
            'lean-stream-accepted-and-marked-counts',
            FIXTURE.n,
            FIXTURE.m,
            accepted.count == 8 and marked.count == 4 and marked.counted,
            'COUNTED',
            count=accepted.count,
            marked=marked.count,
            pairs=accepted.pairs,
            precision=accepted.precision,
        )
    )

    leaf_ok = True
    leaf_details = []
    for cnf, expected in ((STAR, 9), (FOREST, 25), (HUB, 129), (DEMO, 196)):
        stream = count_cnf_streamed(cnf)
        native = count_cnf_ternary(cnf)
        valid = (
            stream.counted
            and stream.count == native.count == brute_count(cnf) == expected
            and stream.polynomial_terms == 0
            and stream.boolean_vectors == 0
            and stream.same_declared_n
        )
        leaf_ok = leaf_ok and valid
        leaf_details.append(dict(n=cnf.n, m=cnf.m, count=stream.count, pairs=stream.pairs))
    rows.append(
        _row(
            'lean-stream-unit-row-leaf-counts',
            0,
            0,
            leaf_ok,
            'COUNTED',
            count=None,
            cases=leaf_details,
        )
    )

    tight = LeanAlgebraBudgets(max_polynomial_terms=1)
    native_refused = count_cnf_ternary(STAR, budgets=tight)
    stream_counted = count_cnf_streamed(STAR, algebra_budgets=tight)
    rows.append(
        _row(
            'lean-stream-replaces-a-refused-support-table',
            STAR.n,
            STAR.m,
            native_refused.status == 'refused'
            and native_refused.count is None
            and stream_counted.counted
            and stream_counted.count == 9
            and stream_counted.polynomial_terms == 0,
            'COUNTED',
            count=stream_counted.count,
            native_reason=native_refused.reason,
        )
    )

    pair_refused = streamed_query_count(STAR, None, budgets=TaylorBudgets(max_pairs=40))
    evaluation_refused = streamed_query_count(
        STAR, None, budgets=TaylorBudgets(max_pairs=1 << 20, max_assignments=40)
    )
    precision_refused = streamed_query_count(
        STAR,
        None,
        budgets=TaylorBudgets(max_pairs=1 << 30, max_assignments=1 << 30),
        algebra_budgets=LeanAlgebraBudgets(max_precision=1),
    )
    rows.append(
        _row(
            'lean-stream-budget-refusals',
            0,
            0,
            pair_refused.status == 'refused'
            and pair_refused.count is None
            and pair_refused.pairs == 0
            and evaluation_refused.status == 'refused'
            and evaluation_refused.count is None
            and precision_refused.status == 'refused'
            and precision_refused.count is None
            and precision_refused.count_word_bits == 0,
            'REFUSED',
            count=None,
            pair_reason=pair_refused.reason,
            precision_reason=precision_refused.reason,
        )
    )

    word_ok = (
        accepted.count_word_bits <= FIXTURE.n + 1
        and streamed_query_count(FIXTURE, None).count_word_bits <= FIXTURE.n + 1
    )
    rows.append(
        _row(
            'lean-stream-count-word-bounds',
            FIXTURE.n,
            FIXTURE.m,
            word_ok,
            'CHECKED',
            count=None,
            count_word_bits=accepted.count_word_bits,
        )
    )

    free = count_cnf_streamed(CNF(8, []))
    zero = count_cnf_streamed(CNF(8, [[]]))
    tautology = count_cnf_streamed(CNF(6, [[(0, False), (0, True)]]))
    rows.append(
        _row(
            'lean-stream-terminal-cases',
            0,
            0,
            free.counted
            and free.count == 256
            and free.pairs == 0
            and zero.counted
            and zero.count == 0
            and tautology.counted
            and tautology.count == 64
            and tautology.n == 6,
            'COUNTED',
            count=None,
            free_bits=free.count,
            tautology=tautology.count,
        )
    )

    backend_ok = True
    backend_details = []
    for cnf, expected in ((STAR, 9), (FOREST, 25), (DEMO, 196)):
        streamed = count_cnf_streaming(cnf, M=2, leaf_backend='stream')
        native = count_cnf_streaming(cnf, M=2, leaf_backend='native')
        valid = (
            streamed.counted
            and streamed.count == native.count == expected
            and streamed.stream_pairs > 0
        )
        backend_ok = backend_ok and valid
        backend_details.append(
            dict(
                n=cnf.n,
                count=streamed.count,
                leaves=streamed.leaves_counted,
                pairs=streamed.stream_pairs,
            )
        )
    shared = count_cnf_streaming(
        STAR, M=2, leaf_backend='stream', taylor_budgets=TaylorBudgets(max_pairs=40)
    )
    rows.append(
        _row(
            'lean-stream-streamed-leaf-backend',
            0,
            0,
            backend_ok
            and shared.status == 'refused'
            and shared.count is None
            and shared.stream_pairs <= 40,
            'COUNTED',
            count=None,
            cases=backend_details,
            shared_refusal=shared.stream_pairs,
        )
    )

    for row in rows:
        stamp_completion(row)
        print(f'  [AE] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP29_SCOPE}', flush=True)
    return rows
