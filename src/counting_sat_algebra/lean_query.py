"""Finite native ternary counting queries for the Lean-guided frontier.

History rows are zero tests of weighted clause violations over F_3. A marked
query adds the violation of one original clause; an accepted query adds no
mark. Their counts include every bit of the original declared n. We use the
constructed ternary digits and one shared marker cutoff, then sum the sparse
Boolean polynomial. No satisfying-assignment list or truth-table oracle is
used. The complete Lean code search, packing and streamed Taylor evaluator
are still separate implementation obligations.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass
from typing import Mapping, Sequence

from .cnf import CNF, Literal
from .lean_encoding import formula_bit_length
from .lean_ring import (
    AlgebraMeter,
    LeanAlgebraBudgets,
    LeanAlgebraRefusal,
    TernaryIndicator,
    construct_ternary_indicator,
    natural,
)
from .lean_squarefree import SquarefreeAlgebra, SquarefreePolynomial

STEP24_SCOPE = (
    'Lean ternary Newton digits, shared marker truncation and exact accepted/marked '
    'queries through finite squarefree polynomial sums. The native sparse-code search, '
    'greedy selector, packing, streamed Taylor reconstruction and Python bit-time/space '
    'transfer remain open. No Theorem 1.1, no 2^{o(n)} claim.'
)


@dataclass(frozen=True)
class TernaryHashRow:
    """Original clause IDs and canonical prime-subfield weights, in source order."""

    terms: tuple[tuple[int, int], ...]

    def __post_init__(self) -> None:
        result = []
        for index, coefficient in self.terms:
            index = natural(index, 'clause index')
            coefficient = operator.index(coefficient)
            if coefficient not in (0, 1, 2):
                raise ValueError('a ternary row weight must be 0, 1 or 2')
            result.append((index, coefficient))
        object.__setattr__(self, 'terms', tuple(result))


@dataclass(frozen=True)
class TernaryQueryResult:
    status: str
    count: int | None
    reason: str
    n: int
    m: int
    query_kind: str
    input_bits: int = 0
    rows: int = 0
    marked_clause: int | None = None
    fixed_variables: tuple[int, ...] = ()
    precision: int = 0
    modulus: int = 0
    coefficient_word_bits: int = 0
    indicator_constructed: bool = False
    newton_steps: int = 0
    source_polynomials: int = 0
    marker_sums: tuple[int, ...] = ()
    scalar_operations: int = 0
    monomial_products: int = 0
    peak_polynomial_terms: int = 0
    peak_live_coefficients: int = 0
    peak_algebra_bits: int = 0
    count_word_bits: int = 0
    boolean_vectors: int = 0
    same_declared_n: bool = True
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False
    recursive_contraction: bool = False

    @property
    def counted(self) -> bool:
        return self.status == 'counted'


def _violation(
    algebra: SquarefreeAlgebra,
    clause: Sequence[Literal],
    fixed: Mapping[int, int],
) -> SquarefreePolynomial:
    required: dict[int, int] = {}
    for variable, negated in clause:
        bit = int(negated)
        if variable in fixed:
            if fixed[variable] != bit:
                return algebra.zero
        elif variable in required and required[variable] != bit:
            return algebra.zero
        else:
            required[variable] = bit
    result = algebra.one
    for variable, bit in sorted(required.items()):
        literal = algebra.variable(variable)
        result = result * (literal if bit else algebra.one - literal)
    return result


def _indicator_row(
    algebra: SquarefreeAlgebra,
    indicator: TernaryIndicator,
    value: SquarefreePolynomial,
) -> tuple[SquarefreePolynomial, ...]:
    square = value * value
    defect = value * (value - algebra.one) * (value - algebra.constant(2))
    power = algebra.one
    output = [algebra.zero] * indicator.t
    for s, (a, b, c) in enumerate(indicator.digits):
        if not power.terms:
            break
        digit = algebra.constant(a) + value.scale(b) + square.scale(c)
        output[s] = digit * power
        if s + 1 < indicator.t:
            power = power * defect
    return tuple(output)


def _precision(n: int, limits: LeanAlgebraBudgets, requested: int | None) -> tuple[int, int]:
    bound = 1 << n
    if requested is not None:
        t = natural(requested, 'precision')
        if not t:
            raise ValueError('query precision must be positive')
        if t > limits.max_precision:
            raise LeanAlgebraRefusal('ternary precision budget exceeded')
        modulus = 3**t
        if modulus <= bound:
            raise ValueError('3^precision must exceed 2^n to reconstruct an exact count')
        return t, modulus
    t, modulus = 0, 1
    while modulus <= bound:
        if t >= limits.max_precision:
            raise LeanAlgebraRefusal(
                'exact count reconstruction exceeds the ternary precision budget'
            )
        t, modulus = t + 1, 3 * modulus
    return t, modulus


def count_ternary_query(
    cnf: CNF,
    rows: Sequence[TernaryHashRow] = (),
    *,
    marked_clause: int | None = None,
    fixed: Mapping[int, int] | None = None,
    precision: int | None = None,
    budgets: LeanAlgebraBudgets | None = None,
) -> TernaryQueryResult:
    """Count history acceptance, optionally intersected with one clause violation.

    ``fixed`` conditions ambient coordinates without renumbering the source.
    A marked query is a violation score, so its count is not the CNF's #SAT.
    Budget refusals publish no count, including after partial algebra work.
    """
    return _count_query(cnf, rows, marked_clause, fixed, precision, budgets)


def count_cnf_ternary(
    cnf: CNF,
    *,
    precision: int | None = None,
    budgets: LeanAlgebraBudgets | None = None,
) -> TernaryQueryResult:
    """Identity history rows give an exact CNF leaf; rows are streamed one at a time."""
    return _count_query(cnf, None, None, None, precision, budgets)


def _count_query(
    cnf: CNF,
    rows: Sequence[TernaryHashRow] | None,
    mark: int | None,
    fixed: Mapping[int, int] | None,
    precision: int | None,
    budgets: LeanAlgebraBudgets | None,
) -> TernaryQueryResult:
    limits = budgets or LeanAlgebraBudgets()
    meter = AlgebraMeter(limits)
    n = natural(cnf.n, 'n')
    kind = 'cnf' if rows is None else ('marked' if mark is not None else 'accepted')
    receipt: dict = dict(n=n, m=cnf.m, query_kind=kind)
    source_polynomials = 0

    def finish(status: str, count: int | None, reason: str) -> TernaryQueryResult:
        return TernaryQueryResult(
            status,
            count,
            reason,
            **receipt,
            source_polynomials=source_polynomials,
            scalar_operations=meter.scalar_operations,
            monomial_products=meter.monomial_products,
            peak_polynomial_terms=meter.peak_polynomial_terms,
            peak_live_coefficients=meter.peak_live_coefficients,
            peak_algebra_bits=meter.peak_algebra_bits,
            count_word_bits=0 if count is None else max(1, count.bit_length()),
        )

    try:
        if n > limits.max_variables or cnf.m > limits.max_clauses:
            raise LeanAlgebraRefusal('query variable or clause budget exceeded')
        if sum(map(len, cnf.clauses)) > limits.max_literals:
            raise LeanAlgebraRefusal('query literal occurrence budget exceeded')
        input_bits = formula_bit_length(cnf, canonical=False)
        receipt['input_bits'] = input_bits
        if input_bits > limits.max_formula_bits:
            raise LeanAlgebraRefusal('query formula encoding bit budget exceeded')
        clauses = tuple(tuple(clause) for clause in cnf.clauses)
        fixed_bits: dict[int, int] = {}
        for variable, bit in (fixed or {}).items():
            variable, bit = natural(variable, 'fixed variable'), operator.index(bit)
            if variable >= n or bit not in (0, 1):
                raise ValueError('fixed coordinates must be declared variables with Boolean values')
            fixed_bits[variable] = bit
        receipt['fixed_variables'] = tuple(sorted(fixed_bits))
        if mark is not None:
            mark = natural(mark, 'marked clause')
            if mark >= len(clauses):
                raise ValueError('marked clause index exceeds the original clause list')
        receipt['marked_clause'] = mark
        row_count = len(clauses) if rows is None else len(rows)
        receipt['rows'] = row_count
        if row_count > limits.max_query_rows:
            raise LeanAlgebraRefusal('query row budget exceeded')
        if rows is not None:
            for row in rows:
                if not isinstance(row, TernaryHashRow):
                    raise ValueError('history rows must be TernaryHashRow values')
                if any(index >= len(clauses) for index, _ in row.terms):
                    raise ValueError('history row clause index exceeds the original clause list')
        row_terms = len(clauses) if rows is None else sum(len(row.terms) for row in rows)
        if row_terms > limits.max_row_terms:
            raise LeanAlgebraRefusal('query weighted source term budget exceeded')
        if rows is None and any(not clause for clause in clauses):
            return finish('counted', 0, 'an original empty clause has no satisfying assignments')
        if row_count == 0 and mark is None:
            return finish(
                'counted',
                1 << (n - len(fixed_bits)),
                'empty history retains every free ambient bit',
            )
        t, modulus = _precision(n, limits, precision)
        receipt.update(precision=t, modulus=modulus, coefficient_word_bits=2 * t)
        # An identity row is generated and discarded; only one such row is parked.
        row_bits = (
            (8 * max(len(clauses) - 1, 0) + 24)
            if rows is None
            else sum(8 * index + 8 * weight + 16 for row in rows for index, weight in row.terms)
        )
        row_bits += sum(8 * v + 16 for v in fixed_bits)
        workspace = 60 * t + 32
        meter.storage(workspace, input_bits + row_bits + workspace * (2 * t))
        indicator = construct_ternary_indicator(t, meter=meter)
        receipt.update(indicator_constructed=True, newton_steps=indicator.newton_steps)
        parked = 6 * t
        algebra = SquarefreeAlgebra(
            n,
            modulus,
            meter,
            coefficient_bits=2 * t,
            parked_coefficients=parked,
            parked_bits=input_bits + row_bits + parked * (2 * t),
        )
        accumulated = (algebra.one,) + (algebra.zero,) * (t - 1)
        source_rows = (
            (TernaryHashRow(((i, 1),)) for i in range(len(clauses))) if rows is None else iter(rows)
        )
        for source in source_rows:
            value = algebra.zero
            for index, weight in source.terms:
                if weight:
                    clause_value = _violation(algebra, clauses[index], fixed_bits)
                    source_polynomials += 1
                    value = value + clause_value.scale(weight)
                    del clause_value
            row = _indicator_row(algebra, indicator, value)
            accumulated = algebra.marker_product(accumulated, row)
            del row, value
        if mark is not None:
            violation = _violation(algebra, clauses[mark], fixed_bits)
            source_polynomials += 1
            accumulated = tuple(coefficient * violation for coefficient in accumulated)
            del violation
        algebra.reserve_scalars(t + 6)
        marker_sums = tuple(
            algebra.sum_boolean(coefficient, fixed=frozenset(fixed_bits))
            for coefficient in accumulated
        )
        receipt['marker_sums'] = marker_sums
        meter.operations(t)
        count = 0
        for value in marker_sums:
            count = (count + value) % modulus
        if count > (1 << (n - len(fixed_bits))):
            raise ArithmeticError(
                'the reconstructed query count exceeds its ambient assignment bound'
            )
        return finish('counted', count, 'shared-marker readout of the exact sparse Boolean sum')
    except LeanAlgebraRefusal as exc:
        return finish('refused', None, f'{exc}; no partial count returned')
