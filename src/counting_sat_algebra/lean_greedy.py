"""Streamed greedy hash-row selection from the supplied Lean program.

``HashGreedyCNF`` starts from the assembled sparse-code rows and, at every
round, prepends the candidate row of least violation mass over the canonical
ternary family of ``HashCandidates``. The mass of a history is the
accumulated marked-query score over every original clause: the total number
of (accepted assignment, violated clause) pairs. When the mass reaches zero
the accepted set is exactly the CNF model set, and the final accepted query
is executed separately from every trial query. Candidates are regenerated
rank by rank, so no candidate library is stored. The machine runs the fixed
``loopCap`` rounds; once the mass is zero, further rounds only prepend
zero-mass rows and cannot shrink the accepted set, so this port stops at the
first zero-mass history. The Lean schedule's candidate depth is
``distanceDenom`` and lies beyond the Python budgets, so the depth is an
explicit finite parameter that refuses when it exceeds its budget. The
high-variable loop, counting extensions, placement and the streamed Taylor
leaf remain open obligations.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace

from .cnf import CNF
from .lean_code import CodeBudgets, canonical_hash_candidates, construct_sparse_code
from .lean_query import TernaryHashRow, count_ternary_query
from .lean_ring import LeanAlgebraBudgets, natural

STEP25_SCOPE = (
    'Scheduled sparse ternary code construction with row regeneration and '
    'violation-mass greedy selection scored by real marked queries; the final '
    'accepted query runs separately. The high-variable loop, counting extensions, '
    'placement, streamed Taylor leaves and the Python bit-time/space transfer '
    'remain open. No Theorem 1.1, no 2^{o(n)} claim.'
)


@dataclass(frozen=True)
class GreedyBudgets:
    depth: int = 1
    max_candidates: int = 300_000
    max_rounds: int | None = None
    max_marked_queries: int = 1_000_000
    max_history_rows: int = 65_536
    max_history_terms: int = 2_000_000

    def __post_init__(self) -> None:
        for item in fields(self):
            value = getattr(self, item.name)
            if item.name == 'max_rounds' and value is None:
                continue
            value = natural(value, item.name)
            object.__setattr__(self, item.name, value)


@dataclass(frozen=True)
class GreedyCountResult:
    status: str
    count: int | None
    reason: str
    n: int = 0
    m: int = 0
    depth: int = 0
    occurrence_bound: int = 0
    max_rounds: int = 0
    code_status: str = ''
    code_candidates: int = 0
    code_word_checks: int = 0
    code_character_checks: int = 0
    code_operations: int = 0
    code_peak_bits: int = 0
    relative_distance_numerator: int = 0
    relative_distance_denominator: int = 1
    minimum_weight: int = 0
    kernel_depth: int = 0
    initial_mass: int = 0
    final_mass: int = 0
    rounds: int = 0
    candidates_scanned: int = 0
    marked_queries: int = 0
    accepted_queries: int = 0
    history_rows: int = 0
    precision: int = 0
    modulus: int = 0
    marker_sums: tuple[int, ...] = ()
    source_polynomials: int = 0
    scalar_operations: int = 0
    monomial_products: int = 0
    peak_polynomial_terms: int = 0
    peak_live_coefficients: int = 0
    peak_algebra_bits: int = 0
    max_precision: int = 0
    count_word_bits: int = 0
    boolean_vectors: int = 0
    same_declared_n: bool = True
    rows: tuple[TernaryHashRow, ...] = ()
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False
    recursive_contraction: bool = False

    @property
    def counted(self) -> bool:
        return self.status == 'counted'

    @property
    def selected(self) -> bool:
        return self.status in ('selected', 'counted')


class _GreedyRefused(Exception):
    """An explicit limit stops selection before a completed count."""


class _Ledger:
    """Accumulate per-query algebra receipts against one shared budget."""

    def __init__(self, budgets: LeanAlgebraBudgets):
        self.budgets = budgets
        self.products = 0
        self.operations = 0
        self.peak_terms = 0
        self.peak_live = 0
        self.peak_bits = 0
        self.max_precision = 0

    def remaining(self) -> LeanAlgebraBudgets:
        return replace(
            self.budgets,
            max_monomial_products=self.budgets.max_monomial_products - self.products,
            max_scalar_operations=self.budgets.max_scalar_operations - self.operations,
        )

    def record(self, result) -> None:
        self.products += result.monomial_products
        self.operations += result.scalar_operations
        self.peak_terms = max(self.peak_terms, result.peak_polynomial_terms)
        self.peak_live = max(self.peak_live, result.peak_live_coefficients)
        self.peak_bits = max(self.peak_bits, result.peak_algebra_bits)
        self.max_precision = max(self.max_precision, result.precision)


def occurrence_bound(cnf: CNF) -> int:
    """Max number of clauses containing one variable, as the Lean occurrence input."""
    counts: dict[int, int] = {}
    for clause in cnf.clauses:
        for variable in {literal[0] for literal in clause}:
            counts[variable] = counts.get(variable, 0) + 1
    return max(counts.values(), default=1)


def _greedy_core(
    cnf: CNF,
    depth: int,
    occurrence: int | None,
    budgets: GreedyBudgets | None,
    native_budgets: LeanAlgebraBudgets | None,
    code_budgets: CodeBudgets | None,
    *,
    final: bool,
) -> GreedyCountResult:
    limits = budgets or GreedyBudgets()
    n = natural(cnf.n, 'n')
    m = natural(cnf.m, 'clause count')
    depth = natural(limits.depth if depth is None else depth, 'depth')
    if depth < 1:
        raise ValueError('candidate depth must be positive')
    ledger = _Ledger(native_budgets or LeanAlgebraBudgets())
    stats: dict = dict(n=n, m=m, depth=depth)
    rows: tuple[TernaryHashRow, ...] = ()

    def finish(status: str, count: int | None, reason: str) -> GreedyCountResult:
        stats.update(
            final_mass=stats.pop('current_mass', 0),
            history_rows=len(rows),
            precision=stats.get('precision', 0),
            modulus=stats.get('modulus', 0),
            marker_sums=tuple(stats.get('marker_sums', ())),
            source_polynomials=stats.get('source_polynomials', 0),
            scalar_operations=ledger.operations,
            monomial_products=ledger.products,
            peak_polynomial_terms=ledger.peak_terms,
            peak_live_coefficients=ledger.peak_live,
            peak_algebra_bits=ledger.peak_bits,
            max_precision=ledger.max_precision,
            count_word_bits=0 if count is None else max(1, count.bit_length()),
            rows=rows,
        )
        return GreedyCountResult(status, count, reason, **stats)

    stats.update(rounds=0, candidates_scanned=0, marked_queries=0, accepted_queries=0,
                 initial_mass=0, current_mass=0, source_polynomials=0)
    try:
        if m:
            if occurrence is None:
                occurrence = occurrence_bound(cnf)
            occurrence = natural(occurrence, 'occurrence bound')
            if occurrence < 1:
                raise ValueError('the occurrence bound must be positive')
            stats['occurrence_bound'] = occurrence
            if (3 * m) ** depth > limits.max_candidates:
                raise _GreedyRefused('candidate rank budget exceeded')
            code = construct_sparse_code(m, occurrence, budgets=code_budgets)
            stats.update(code_status=code.status, code_candidates=code.candidates,
                         code_word_checks=code.word_checks,
                         code_character_checks=code.character_checks,
                         code_operations=code.operations, code_peak_bits=code.peak_bits)
            if not code.checked:
                raise _GreedyRefused(f'the scheduled sparse code refused: {code.reason}')
            certificate = code.certificate
            stats.update(
                relative_distance_numerator=certificate.relative_distance.numerator,
                relative_distance_denominator=certificate.relative_distance.denominator,
                minimum_weight=certificate.minimum_weight,
                kernel_depth=certificate.depth,
            )
            rows = tuple(certificate.code.rows)
            if len(rows) > limits.max_history_rows:
                raise _GreedyRefused('initial history row budget exceeded')
            if sum(len(row.terms) for row in rows) > limits.max_history_terms:
                raise _GreedyRefused('initial history term budget exceeded')
        max_rounds = limits.max_rounds
        if max_rounds is None:
            max_rounds = 2 * (n + m.bit_length() + 1)
        stats['max_rounds'] = max_rounds

        def mass(history: tuple[TernaryHashRow, ...]) -> int:
            total = 0
            for clause in range(m):
                if stats['marked_queries'] >= limits.max_marked_queries:
                    raise _GreedyRefused('marked query budget exceeded')
                query = count_ternary_query(
                    cnf, history, marked_clause=clause, budgets=ledger.remaining()
                )
                stats['marked_queries'] += 1
                ledger.record(query)
                if not query.counted:
                    raise _GreedyRefused(f'a marked query refused: {query.reason}')
                stats['source_polynomials'] += query.source_polynomials
                total += query.count
            return total

        current = mass(rows)
        stats['initial_mass'] = current
        stats['current_mass'] = current
        while current > 0:
            if stats['rounds'] >= max_rounds:
                raise _GreedyRefused(
                    'the greedy round cap was reached with positive violation mass'
                )
            best: TernaryHashRow | None = None
            best_score: int | None = None
            for _rank, candidate in canonical_hash_candidates(m, depth):
                score = mass((candidate,) + rows)
                stats['candidates_scanned'] += 1
                if best is None or score < best_score:
                    best, best_score = candidate, score
            rows = (best,) + rows
            if len(rows) > limits.max_history_rows:
                raise _GreedyRefused('history row budget exceeded')
            if sum(len(row.terms) for row in rows) > limits.max_history_terms:
                raise _GreedyRefused('history term budget exceeded')
            stats['rounds'] += 1
            current = stats['current_mass'] = best_score
        if final:
            query = count_ternary_query(cnf, rows, budgets=ledger.remaining())
            stats['accepted_queries'] = 1
            ledger.record(query)
            if not query.counted:
                raise _GreedyRefused(f'the final accepted query refused: {query.reason}')
            stats['source_polynomials'] += query.source_polynomials
            stats.update(precision=query.precision, modulus=query.modulus,
                         marker_sums=query.marker_sums)
            if query.count > (1 << n):
                raise ArithmeticError('the accepted query count exceeds its assignment bound')
            return finish(
                'counted', query.count,
                f'zero violation mass after {stats["rounds"]} greedy rounds; '
                'the separate accepted query counts the exact model set',
            )
        return finish(
            'selected', None,
            f'zero violation mass after {stats["rounds"]} greedy rounds',
        )
    except _GreedyRefused as exc:
        return finish('refused', None, f'{exc}; no partial count returned')


def select_history_rows(
    cnf: CNF,
    *,
    depth: int | None = None,
    occurrence: int | None = None,
    budgets: GreedyBudgets | None = None,
    native_budgets: LeanAlgebraBudgets | None = None,
    code_budgets: CodeBudgets | None = None,
) -> GreedyCountResult:
    """Run the greedy rounds only; the final accepted query is not executed."""
    return _greedy_core(cnf, depth, occurrence, budgets, native_budgets, code_budgets, final=False)


def count_cnf_greedy(
    cnf: CNF,
    *,
    depth: int | None = None,
    occurrence: int | None = None,
    budgets: GreedyBudgets | None = None,
    native_budgets: LeanAlgebraBudgets | None = None,
    code_budgets: CodeBudgets | None = None,
) -> GreedyCountResult:
    """Select history rows by real marked scores, then count the accepted query.

    The quality gate is the accumulated violation mass: it reaches zero
    exactly when the accepted set has no clause-violating assignment, so the
    separate final query counts exactly the CNF models. A refusal returns no
    count, including after partial marked-query work.
    """
    return _greedy_core(cnf, depth, occurrence, budgets, native_budgets, code_budgets, final=True)
