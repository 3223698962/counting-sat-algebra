"""Streamed high-variable loops from the supplied Lean program.

``HighQueryProgram`` conditions every query on the early high-incidence set:
the coordinates whose clause degree exceeds the paid cutoff
``ceil(totalScope * threshold / ambient)`` with the square-root threshold of
``EarlyThreshold``. The loop walks all ``2^|H|`` assignments of the high
coordinates through a binary odometer, runs the conditioned query at each
rank and accumulates the exact totals. ``rankSum_fullCount`` makes the
accumulated total equal the unconditioned query count on the original
declared n. The finite port bounds the rank count per query and shares the
algebra work; a refusal publishes no count.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
from math import isqrt

from .cnf import CNF
from .lean_query import TernaryHashRow, count_ternary_query
from .lean_ring import LeanAlgebraBudgets, natural

STEP26_SCOPE = (
    'Early high-incidence threshold, cutoff arithmetic, the high-variable frontend '
    'and the streamed conditioned rank loop through real queries. The counting '
    'extension ring, placement, streamed Taylor leaves and the Python bit-time/space '
    'transfer remain open. No Theorem 1.1, no 2^{o(n)} claim.'
)


def early_log(ambient: int) -> int:
    """``Asymptotic.L``: floor(log2 (ambient + 16)) + 1 without floating point."""
    ambient = natural(ambient, 'ambient')
    return (ambient + 16).bit_length()


def early_threshold(ambient: int) -> int:
    """``EarlyThreshold.threshold``: the paid integer square-root threshold."""
    return isqrt(early_log(ambient)) + 1


@dataclass(frozen=True)
class HighFrontend:
    ambient: int
    threshold: int
    total_scope: int
    cutoff: int
    high_variables: tuple[int, ...]
    degrees: tuple[int, ...]

    @property
    def cardinality(self) -> int:
        return len(self.high_variables)

    @property
    def ranks(self) -> int:
        return 1 << len(self.high_variables)

    def mask(self) -> tuple[int, ...]:
        """Ambient-order bits: coordinate v is high exactly when its bit is one."""
        high = set(self.high_variables)
        return tuple(int(v in high) for v in range(self.ambient))


def high_frontend(cnf: CNF, *, threshold: int | None = None) -> HighFrontend:
    """``HighFrontendMask``: the exact cutoff and the coordinates above it."""
    ambient = natural(cnf.n, 'n')
    scopes = [sorted({variable for variable, _ in clause}) for clause in cnf.clauses]
    total_scope = sum(len(scope) for scope in scopes)
    if threshold is None:
        threshold = early_threshold(ambient)
    else:
        threshold = natural(threshold, 'threshold')
        if threshold < 1:
            raise ValueError('the high-incidence threshold must be positive')
    degrees = [0] * ambient
    for scope in scopes:
        for variable in scope:
            degrees[variable] += 1
    if ambient:
        cutoff = (total_scope * threshold + ambient - 1) // ambient
        high = tuple(v for v in range(ambient) if cutoff < degrees[v])
    else:
        cutoff = 0
        high = ()
    return HighFrontend(ambient, threshold, total_scope, cutoff, high, tuple(degrees))


@dataclass(frozen=True)
class HighBudgets:
    max_ranks: int = 1 << 20
    max_queries: int = 1_000_000

    def __post_init__(self) -> None:
        for item in fields(self):
            object.__setattr__(self, item.name, natural(getattr(self, item.name), item.name))


@dataclass(frozen=True)
class HighLoopResult:
    status: str
    count: int | None
    reason: str
    n: int = 0
    threshold: int = 0
    total_scope: int = 0
    cutoff: int = 0
    high_cardinality: int = 0
    high_variables: tuple[int, ...] = ()
    ranks: int = 0
    queries: int = 0
    accumulator_bits: int = 0
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
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False
    recursive_contraction: bool = False

    @property
    def counted(self) -> bool:
        return self.status == 'counted'


class _Ledger:
    def __init__(self, budgets: LeanAlgebraBudgets):
        self.budgets = budgets
        self.products = 0
        self.operations = 0
        self.peak_terms = 0
        self.peak_live = 0
        self.peak_bits = 0
        self.max_precision = 0
        self.source_polynomials = 0

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
        self.source_polynomials += result.source_polynomials


def count_cnf_high_loop(
    cnf: CNF,
    rows: tuple[TernaryHashRow, ...] | None = None,
    *,
    marked_clause: int | None = None,
    fixed: dict[int, int] | None = None,
    threshold: int | None = None,
    budgets: HighBudgets | None = None,
    native_budgets: LeanAlgebraBudgets | None = None,
) -> HighLoopResult:
    """Accumulate the conditioned query over every high assignment.

    ``rows=None`` streams one unit history row per original clause, so the
    accumulation counts the CNF models exactly like ``count_cnf_ternary``;
    an explicit tuple evaluates those history rows instead. Ranks increase
    with the binary odometer over the sorted high coordinates and every rank
    runs one real conditioned query. The accumulated total equals the
    unconditioned query count on the original declared n. Budget refusals
    return no count, including after completed ranks.
    """
    limits = budgets or HighBudgets()
    frontend = high_frontend(cnf, threshold=threshold)
    ledger = _Ledger(native_budgets or LeanAlgebraBudgets())
    base_fixed: dict[int, int] = {}
    for variable, bit in (fixed or {}).items():
        variable, bit = natural(variable, 'fixed variable'), int(bit)
        if variable >= frontend.ambient or bit not in (0, 1):
            raise ValueError('fixed coordinates must be declared variables with Boolean values')
        base_fixed[variable] = bit
    if set(base_fixed) & set(frontend.high_variables):
        raise ValueError('a pre-fixed coordinate cannot also be a high coordinate')
    total = 0
    queries = 0
    accumulator_bits = max(1, (1 << frontend.ambient).bit_length())

    def finish(status: str, count: int | None, reason: str) -> HighLoopResult:
        return HighLoopResult(
            status,
            count,
            reason,
            n=frontend.ambient,
            threshold=frontend.threshold,
            total_scope=frontend.total_scope,
            cutoff=frontend.cutoff,
            high_cardinality=frontend.cardinality,
            high_variables=frontend.high_variables,
            ranks=frontend.ranks,
            queries=queries,
            accumulator_bits=max(1, accumulator_bits if count is None else count.bit_length()),
            precision=ledger.max_precision,
            modulus=3**ledger.max_precision if ledger.max_precision else 0,
            source_polynomials=ledger.source_polynomials,
            scalar_operations=ledger.operations,
            monomial_products=ledger.products,
            peak_polynomial_terms=ledger.peak_terms,
            peak_live_coefficients=ledger.peak_live,
            peak_algebra_bits=ledger.peak_bits,
            max_precision=ledger.max_precision,
            count_word_bits=0 if count is None else max(1, count.bit_length()),
        )

    if frontend.ranks > limits.max_ranks:
        return finish('refused', None, 'high rank budget exceeded; no partial count returned')
    if frontend.ranks > limits.max_queries:
        return finish('refused', None, 'high query budget exceeded; no partial count returned')
    bound = 1 << (frontend.ambient - len(base_fixed))
    for rank in range(frontend.ranks):
        conditioned = dict(base_fixed)
        for position, variable in enumerate(frontend.high_variables):
            conditioned[variable] = (rank >> position) & 1
        query = count_ternary_query(
            cnf,
            rows,
            marked_clause=marked_clause,
            fixed=conditioned,
            budgets=ledger.remaining(),
        )
        queries += 1
        ledger.record(query)
        if not query.counted:
            return finish(
                'refused', None,
                f'rank {rank} refused: {query.reason}; no partial count returned',
            )
        total += query.count
        if total > bound:
            raise ArithmeticError('the accumulation exceeds its ambient assignment bound')
        accumulator_bits = max(accumulator_bits, total.bit_length())
    return finish(
        'counted', total,
        f'accumulated {frontend.ranks} high ranks through conditioned queries on the original n',
    )
