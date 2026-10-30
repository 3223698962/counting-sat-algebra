"""Streamed pointwise correction-pair evaluation of the native ternary queries.

The native leaf builds full squarefree support tables before summing Boolean
values. ``BooleanDifferences`` offers the alternative the Lean program runs:
a stream over correction pairs ``(E, U)`` with ``|E| < J``, where every pair
contributes exactly one pointwise evaluation. For a polynomial of X-degree
below ``J`` the identity is exact:

    sum_x p(x) = sum_{(E,U)} (-1)^(|E|-|U|) * 2^(free-|E|) * p(U)

Here ``p`` is the marker-sum polynomial of the query. Its degree is bounded
by ``(sum of row degrees) * (3*(t-1) + 2)``, where a row degree is the
largest post-substitution distinct-variable count of its clauses, so
``J = that bound + 1`` suffices. Evaluating ``p(U)`` needs only the scalar
row values at one point: no support table is ever materialized. Pair and
point-evaluation budgets refuse before any stream work; a refusal publishes
no count. The full machine's degree reduction, the native placement and the
Python bit-time/space transfer remain open.
"""

from __future__ import annotations

from dataclasses import dataclass

from .cnf import CNF, Literal
from .lean_encoding import formula_bit_length
from .lean_query import TernaryHashRow
from .lean_ring import (
    LeanAlgebraBudgets,
    LeanAlgebraRefusal,
    construct_ternary_indicator,
    natural,
)
from .lean_taylor import TaylorBudgets, TaylorRefusal, decode_pair, pair_count, pair_walk

STEP29_SCOPE = (
    'Streamed pointwise correction-pair evaluation of the native queries without '
    'squarefree support tables, plus the streamed CNF leaf backend. Costs scale with '
    'pairs times point evaluations; the full machine degree reduction, the native '
    'placement and the Python bit-time/space transfer remain open. '
    'No Theorem 1.1, no 2^{o(n)} claim.'
)


@dataclass(frozen=True)
class StreamedQueryResult:
    status: str
    count: int | None
    reason: str
    n: int = 0
    m: int = 0
    query_kind: str = ''
    rows: int = 0
    row_terms: int = 0
    input_bits: int = 0
    marked_clause: int | None = None
    fixed_variables: tuple[int, ...] = ()
    precision: int = 0
    modulus: int = 0
    degree_bound: int = 0
    order: int = 0
    pairs: int = 0
    point_evaluations: int = 0
    marker_sums: tuple[int, ...] = ()
    count_word_bits: int = 0
    boolean_vectors: int = 0
    polynomial_terms: int = 0
    same_declared_n: bool = True
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False
    recursive_contraction: bool = False

    @property
    def counted(self) -> bool:
        return self.status == 'counted'


def _violates(clause: tuple[Literal, ...], point: dict[int, int]) -> int:
    """One iff every literal of the clause is false at the point."""
    for variable, negated in clause:
        if point[variable] != int(negated):
            return 0
    return 1


def _clause_degree(clause: tuple[Literal, ...], fixed: dict[int, int]) -> int:
    """The X-degree of one clause violation polynomial after fixed substitution."""
    polarity: dict[int, int] = {}
    for variable, negated in clause:
        bit = int(negated)
        if variable in polarity and polarity[variable] != bit:
            return 0
        polarity[variable] = bit
    required: set[int] = set()
    for variable, bit in polarity.items():
        if variable in fixed:
            if fixed[variable] != bit:
                return 0
        else:
            required.add(variable)
    return len(required)


def _row_degree(row: TernaryHashRow, clauses, fixed: dict[int, int]) -> int:
    degree = 0
    for index, weight in row.terms:
        if weight:
            degree = max(degree, _clause_degree(clauses[index], fixed))
    return degree


def _query_precision(free_bits: int, limits: LeanAlgebraBudgets, requested):
    bound = 1 << free_bits
    if requested is not None:
        requested = natural(requested, 'precision')
        if not requested:
            raise ValueError('query precision must be positive')
        if requested > limits.max_precision:
            raise LeanAlgebraRefusal('ternary precision budget exceeded')
        modulus = 3**requested
        if modulus <= bound:
            raise ValueError('3^precision must exceed 2^free_bits for an exact count')
        return requested, modulus
    steps, modulus = 0, 1
    while modulus <= bound:
        if steps >= limits.max_precision:
            raise LeanAlgebraRefusal(
                'exact count reconstruction exceeds the ternary precision budget'
            )
        steps, modulus = steps + 1, 3 * modulus
    return steps, modulus


def streamed_query_count(
    cnf: CNF,
    rows=None,
    *,
    marked_clause=None,
    fixed=None,
    precision=None,
    budgets=None,
    algebra_budgets=None,
) -> StreamedQueryResult:
    """Count a ternary query by streaming correction pairs with point evaluations.

    ``rows=None`` streams one unit history row per original clause. Every pair
    evaluates the marker-sum polynomial at exactly one Boolean point using the
    scalar row values; no squarefree support table is built. Pair and
    point-evaluation budgets are preflighted, and a refusal returns no count.
    """
    limits = budgets or TaylorBudgets()
    algebra = algebra_budgets or LeanAlgebraBudgets()
    n = natural(cnf.n, 'n')
    m = natural(cnf.m, 'clause count')
    clauses = tuple(tuple(clause) for clause in cnf.clauses)
    kind = 'cnf' if rows is None else ('marked' if marked_clause is not None else 'accepted')
    receipt: dict = dict(n=n, m=m, query_kind=kind)
    try:
        if n > algebra.max_variables or m > algebra.max_clauses:
            raise LeanAlgebraRefusal('query variable or clause budget exceeded')
        if sum(map(len, clauses)) > algebra.max_literals:
            raise LeanAlgebraRefusal('query literal occurrence budget exceeded')
        input_bits = formula_bit_length(cnf, canonical=False)
        receipt['input_bits'] = input_bits
        if input_bits > algebra.max_formula_bits:
            raise LeanAlgebraRefusal('query formula encoding bit budget exceeded')
        fixed_bits: dict[int, int] = {}
        for variable, bit in (fixed or {}).items():
            variable, bit = natural(variable, 'fixed variable'), int(bit)
            if variable >= n or bit not in (0, 1):
                raise ValueError(
                    'fixed coordinates must be declared variables with Boolean values'
                )
            fixed_bits[variable] = bit
        receipt['fixed_variables'] = tuple(sorted(fixed_bits))
        if marked_clause is not None:
            marked_clause = natural(marked_clause, 'marked clause')
            if marked_clause >= len(clauses):
                raise ValueError('marked clause index exceeds the original clause list')
        receipt['marked_clause'] = marked_clause
        row_list = (
            tuple(TernaryHashRow(((index, 1),)) for index in range(m))
            if rows is None
            else tuple(rows)
        )
        for row in row_list:
            if not isinstance(row, TernaryHashRow):
                raise ValueError('history rows must be TernaryHashRow values')
            if any(index >= len(clauses) for index, _ in row.terms):
                raise ValueError('history row clause index exceeds the original clause list')
        row_terms = sum(len(row.terms) for row in row_list)
        receipt.update(rows=len(row_list), row_terms=row_terms)
        if len(row_list) > algebra.max_query_rows:
            raise LeanAlgebraRefusal('query row budget exceeded')
        if row_terms > algebra.max_row_terms:
            raise LeanAlgebraRefusal('query weighted source term budget exceeded')
        if rows is None and any(not clause for clause in clauses):
            return _finish('counted', 0, 'an original empty clause has no satisfying assignments', receipt)
        if not row_list and marked_clause is None:
            return _finish(
                'counted', 1 << (n - len(fixed_bits)),
                'empty history retains every free ambient bit', receipt,
            )
        t, modulus = _query_precision(n, algebra, precision)
        receipt.update(precision=t, modulus=modulus)
        indicator = construct_ternary_indicator(t, budgets=algebra)
        free = tuple(variable for variable in range(n) if variable not in fixed_bits)
        degree_total = sum(_row_degree(row, clauses, fixed_bits) for row in row_list)
        mark_degree = (
            _clause_degree(clauses[marked_clause], fixed_bits)
            if marked_clause is not None
            else 0
        )
        order = degree_total * (3 * (t - 1) + 2) + mark_degree + 1
        receipt.update(degree_bound=degree_total, order=order)
        total_pairs = pair_count(len(free), order)
        if total_pairs > limits.max_pairs:
            raise TaylorRefusal('correction pair budget exceeded')
        if total_pairs > limits.max_assignments:
            raise TaylorRefusal('point evaluation budget exceeded')
        marker_sums = [0] * t
        pairs = 0
        digits = indicator.digits
        for word in pair_walk(len(free), order):
            outer, inner = decode_pair(word)
            point = dict(fixed_bits)
            for position, variable in enumerate(free):
                point[variable] = 1 if position in inner else 0
            product = [0] * t
            product[0] = 1
            for row in row_list:
                value = 0
                for index, weight in row.terms:
                    if weight and _violates(clauses[index], point):
                        value += weight
                defect = value * (value - 1) * (value - 2) % modulus
                row_vector = [0] * t
                power = 1
                for s, (a, b, c) in enumerate(digits):
                    row_vector[s] = (a + b * value + c * value * value) % modulus * power % modulus
                    power = power * defect % modulus
                convoluted = [0] * t
                for u, cu in enumerate(product):
                    if not cu:
                        continue
                    for v in range(t - u):
                        if row_vector[v]:
                            convoluted[u + v] = (
                                convoluted[u + v] + cu * row_vector[v]
                            ) % modulus
                product = convoluted
            if marked_clause is not None and not _violates(clauses[marked_clause], point):
                product = [0] * t
            sign = -1 if (len(outer) - len(inner)) % 2 else 1
            weight = (1 << (len(free) - len(outer))) % modulus
            scale = sign * weight
            for s in range(t):
                marker_sums[s] = (marker_sums[s] + scale * product[s]) % modulus
            pairs += 1
        receipt.update(
            pairs=pairs, point_evaluations=pairs, marker_sums=tuple(marker_sums)
        )
        count = 0
        for value in marker_sums:
            count = (count + value) % modulus
        if count > (1 << (n - len(fixed_bits))):
            raise ArithmeticError(
                'the reconstructed query count exceeds its ambient assignment bound'
            )
        return _finish(
            'counted', count,
            'pointwise correction-pair readout without any squarefree support table',
            receipt,
        )
    except (LeanAlgebraRefusal, TaylorRefusal) as exc:
        return _finish('refused', None, f'{exc}; no partial count returned', receipt)


def _finish(status: str, count: int | None, reason: str, receipt: dict) -> StreamedQueryResult:
    return StreamedQueryResult(
        status,
        count,
        reason,
        **receipt,
        count_word_bits=0 if count is None else max(1, count.bit_length()),
    )


def count_cnf_streamed(
    cnf: CNF,
    *,
    precision=None,
    budgets=None,
    algebra_budgets=None,
) -> StreamedQueryResult:
    """Identity unit history rows give an exact CNF leaf without any support table."""
    return streamed_query_count(
        cnf, None, precision=precision, budgets=budgets, algebra_budgets=algebra_budgets
    )
