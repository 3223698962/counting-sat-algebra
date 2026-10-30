"""Step 24 finite witnesses for native ternary accepted and marked queries."""

from __future__ import annotations

from itertools import product
from math import comb

from .protocol import _row, stamp_completion
from ..cnf import CNF, clause_poly_bool
from ..lean_query import STEP24_SCOPE, TernaryHashRow, count_cnf_ternary, count_ternary_query
from ..lean_ring import AlgebraMeter, LeanAlgebraBudgets, TruncatedRing, construct_ternary_indicator
from ..lean_sparsify import count_cnf_streaming
from ..lean_squarefree import SquarefreeAlgebra


def _oracle(cnf, rows, marked=None, fixed=None):
    total = 0
    for bits in product((0, 1), repeat=cnf.n):
        if any(bits[v] != bit for v, bit in (fixed or {}).items()):
            continue
        bad = [1 - clause_poly_bool(bits, clause) for clause in cnf.clauses]
        accepted = all(
            sum(weight * bad[index] for index, weight in row.terms) % 3 == 0 for row in rows
        )
        total += accepted and (marked is None or bad[marked])
    return total


def _count_row(name, cnf, result, expected):
    return _row(
        name,
        cnf.n,
        cnf.m,
        result.counted
        and result.count == expected
        and result.boolean_vectors == 0
        and result.same_declared_n,
        'COUNTED',
        count=result.count,
        expected=expected,
        boolean_vectors=result.boolean_vectors,
        query_kind=result.query_kind,
        marked_clause=result.marked_clause,
        fixed_variables=result.fixed_variables,
        precision=result.precision,
        modulus=result.modulus,
        marker_sums=result.marker_sums,
        source_polynomials=result.source_polynomials,
        monomial_products=result.monomial_products,
        scalar_operations=result.scalar_operations,
        peak_polynomial_terms=result.peak_polynomial_terms,
        peak_live_coefficients=result.peak_live_coefficients,
        peak_algebra_bits=result.peak_algebra_bits,
        count_word_bits=result.count_word_bits,
    )


def run_lean_native_suite() -> list[dict]:
    print('=== (Z) Lean native ternary algebra and counting queries ===', flush=True)
    rows = []
    indicator = construct_ternary_indicator(6)
    hermite = all(
        sum(
            comb(j, derivative) * c * root ** (j - derivative)
            for j, c in enumerate(indicator.polynomial)
            if j >= derivative
        )
        % indicator.modulus
        == int(root == 0 and derivative == 0)
        for root in range(3)
        for derivative in range(6)
    )
    rows.append(
        _row(
            'lean-ternary-newton-hermite-coordinates',
            0,
            0,
            hermite and all(indicator.evaluate(x) == int(x % 3 == 0) for x in range(-81, 82)),
            'CHECKED',
            count=None,
            precision=6,
            scalar_words=18,
            newton_steps=indicator.newton_steps,
            scalar_operations=indicator.scalar_operations,
        )
    )
    small = construct_ternary_indicator(2)
    rows.append(
        _row(
            'lean-ternary-fixed-width-coefficient-words',
            0,
            0,
            small.digits == ((1, 3, 5), (3, 7, 3)) and len(small.encode()) == 24,
            'CHECKED',
            count=None,
            coefficient_word_bits=4,
            encoded_bits=len(small.encode()),
        )
    )
    ring = TruncatedRing(2, 3, 27)
    U, Z = ring.monomial(1, 0), ring.monomial(0, 1)
    rows.append(
        _row(
            'lean-shared-marker-readout-after-product',
            0,
            0,
            (U * U).readout() == (0, 0, 0)
            and U.readout() == (1, 0, 0)
            and (U * Z).readout() == (0, 1, 0)
            and (Z**3).readout() == (0, 0, 0),
            'CHECKED',
            count=None,
            dense_slots=ring.slots,
            readout_is_multiplicative=False,
        )
    )
    algebra = SquarefreeAlgebra(3, 27, AlgebraMeter(LeanAlgebraBudgets()))
    X, Y = algebra.variable(0), algebra.variable(1)
    poly = algebra.one - X * Y
    rows.append(
        _row(
            'lean-support-union-boolean-sum',
            3,
            0,
            (X * X).terms == X.terms
            and not (X * (algebra.one - X)).terms
            and algebra.sum_boolean(poly) == 6,
            'CHECKED',
            count=None,
            boolean_polynomial_sum=6,
            boolean_vectors=0,
        )
    )

    formula = CNF(5, [[(v, False)] for v in range(3)])
    history = (TernaryHashRow(((0, 1), (1, 1), (2, 1))),)
    for mark, name, expected in ((None, 'accepted', 8), (0, 'marked-violation', 4)):
        result = count_ternary_query(formula, history, marked_clause=mark)
        row = _count_row(f'lean-ternary-{name}-count', formula, result, expected)
        row['ok'] = row['ok'] and _oracle(formula, history, mark) == expected
        if not row['ok']:
            row['label'] = 'FAIL'
        rows.append(row)
    conditioned = count_ternary_query(formula, history, fixed={3: 0})
    rows.append(
        _count_row(
            'lean-ternary-high-conditioning-keeps-n',
            formula,
            conditioned,
            _oracle(formula, history, fixed={3: 0}),
        )
    )
    repeated = CNF(4, [[(0, False), (0, False)], [(1, False), (1, True)], []])
    repeated_history = (TernaryHashRow(((0, 1), (0, 2))),)
    rows.append(
        _count_row(
            'lean-ternary-raw-source-repetitions',
            repeated,
            count_ternary_query(repeated, repeated_history, marked_clause=0),
            8,
        )
    )
    demo = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])
    rows.append(_count_row('lean-native-cnf-count-196', demo, count_cnf_ternary(demo), 196))

    star = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
    forest = CNF(9, [[(3 * i, False), (3 * i + j, False)] for i in range(3) for j in (1, 2)])
    for name, cnf, expected, leaves in (('star', star, 9, 2), ('forest', forest, 125, 8)):
        result = count_cnf_streaming(cnf, M=2, leaf_backend='native', total_vector_budget=0)
        rows.append(
            _row(
                f'lean-native-streamed-{name}',
                cnf.n,
                cnf.m,
                result.count == expected
                and result.leaves_counted == leaves
                and result.native_queries == leaves
                and result.completed_leaves_retained == 0
                and result.max_pending_nodes <= result.depth + 1
                and result.boolean_vectors == 0,
                'COUNTED',
                count=result.count,
                leaves=result.leaves_counted,
                native_queries=result.native_queries,
                depth=result.depth,
                max_pending_nodes=result.max_pending_nodes,
                completed_leaves_retained=0,
                boolean_vectors=result.boolean_vectors,
                scalar_operations=result.scalar_operations,
                monomial_products=result.monomial_products,
                peak_polynomial_terms=result.peak_polynomial_terms,
            )
        )
    for cnf, expected, name in ((CNF(8, []), 256, 'free-bits'), (CNF(8, [[]]), 0, 'empty-clause')):
        result = count_cnf_ternary(cnf)
        rows.append(_count_row(f'lean-native-terminal-{name}', cnf, result, expected))
    for name, limits in (
        ('precision', LeanAlgebraBudgets(max_precision=1)),
        ('storage', LeanAlgebraBudgets(max_live_coefficients=1)),
        ('product', LeanAlgebraBudgets(max_monomial_products=0)),
    ):
        result = count_cnf_ternary(star, budgets=limits)
        rows.append(
            _row(
                f'lean-native-{name}-budget-refusal',
                star.n,
                star.m,
                result.status == 'refused'
                and result.count is None
                and result.count_word_bits == 0
                and result.boolean_vectors == 0,
                'REFUSED',
                count=None,
                indicator_constructed=result.indicator_constructed,
                scalar_operations=result.scalar_operations,
                monomial_products=result.monomial_products,
                reason=result.reason,
            )
        )
    result = count_cnf_streaming(
        star,
        M=2,
        leaf_backend='native',
        native_budgets=LeanAlgebraBudgets(max_scalar_operations=1000),
    )
    rows.append(
        _row(
            'lean-native-global-budget-no-partial-total',
            star.n,
            star.m,
            result.status == 'refused'
            and result.count is None
            and result.leaves_counted == 1
            and result.native_queries == 2
            and 0 < result.scalar_operations <= 1000
            and result.count_word_bits == 0,
            'REFUSED',
            count=None,
            leaves_counted=result.leaves_counted,
            native_queries=result.native_queries,
            scalar_operations=result.scalar_operations,
            reason=result.reason,
        )
    )
    for row in rows:
        stamp_completion(row)
        print(f'  [Z] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP24_SCOPE}', flush=True)
    return rows
