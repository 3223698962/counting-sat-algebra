"""Finite Section 6 equation selection and disjoint sparsification.

A clause violation bit is a square-free polynomial over F2: the factor of a
positive literal is ``1+x``, the factor of a negative literal is ``x``, a
tautology is 0, and an empty clause is 1. A binary row is a list of
``(clause index, coefficient)`` pairs. Repeated indices cancel modulo 2, and
the sum expands to a sparse equation compatible with ``sparse_systems.py``.

Lemma 6.1's integers ``μ``, ``P`` and ``D0`` are exact. The matrix ``A(0)`` is
built only when the shift search and the local-code search fit the budgets.
On those small instances the local distance bound is at most 1, so the
deterministic even-parity checks meet it and the shortened Tanner rows are
explicit. Larger parameters are refused before that search is allocated.

Lemma 6.2 appends rows of a given list length by minimizing ``N(A)``. The
candidate set ``([m]×{0,1})^P`` is enumerated only when ``(2m)^P`` fits a
budget. The counter is injected: the default is an F2 cube scan, and a caller
may pass ``count_common_zeros``. A budget failure returns no count.

Lemma 6.3 / Appendix A builds a disjoint sparsification tree on the same
declared ``n``. Normalization deletes repeated literals, tautologies,
duplicate clauses and clauses that contain another clause. It does not
propagate units and it does not copy variables. Step 5's implication-cycle
copying increases ``n``; this tree does not. Every satisfying assignment of a
parent meets exactly one child. Leaf counts are summed by an injected
counter.

Proposition 6.4's integers ``λ_n``, ``ν_n`` and ``s`` are computed by integer
comparisons. ``σ=1/s`` only when ``s≥2``; otherwise the schedule is the
exhaustive-search case. This module does not claim the ``O_k(ρ(n))`` bound,
does not run Appendix D, and does not implement Theorem 1.1.
``theorem_1_1_complete`` and ``lean_proof_complete`` stay false.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass, fields
from fractions import Fraction
from itertools import combinations, product
from typing import Callable, Sequence

from .cnf import CNF, brute_count
from .sparse_systems import (
    Section5Budgets,
    SparsePolynomial,
    SparseSystem,
    count_common_zeros,
    sparse_polynomial,
)

STEP19_SCOPE = (
    'Finite Section 6 checks: violation polynomials, a budgeted Lemma 6.1 '
    'matrix A(0), greedy row selection, disjoint sparsification on the same '
    'n, and Proposition 6.4 integers. Step 5 copying is not this tree. '
    'Appendix D is not run. No O_k(rho(n)) bound, no Theorem 1.1, and no '
    '2^{o(n)} claim. theorem_1_1_complete and lean_proof_complete stay false.'
)

Literal = tuple[int, bool]
ClauseSet = frozenset[Literal]
LeafCounter = Callable[[CNF], int | None]
SystemCounter = Callable[[tuple[SparsePolynomial, ...]], int | None]


class Section6Refusal(ValueError):
    """A Section 6 budget stopped the work. No partial count is available."""


@dataclass(frozen=True)
class Section6Budgets:
    """Explicit caps. Oversized searches, candidate lists and cubes stop first."""

    max_candidate_lists: int = 10_000
    max_rounds: int = 64
    max_counts: int = 100_000
    max_boolean_vectors: int = 4_096
    max_local_search_bits: int = 16
    max_shift_search_bits: int = 16
    max_codeword_tests: int = 4_096
    max_tree_nodes: int = 4_096
    max_branching_bits: int = 4_096

    def __post_init__(self) -> None:
        for field in fields(self):
            value = operator.index(getattr(self, field.name))
            if value < 0:
                raise ValueError('Section 6 budgets must be nonnegative integers')
            object.__setattr__(self, field.name, value)


@dataclass(frozen=True)
class Lemma61Parameters:
    """Exact integers from Lemma 6.1. They are not a running-time bound."""

    n: int
    m: int
    B: int
    mu: int
    P: int
    D0: int
    d_G: int
    q_G: int
    ell_G: int
    c_G: int
    h: int
    epsilon: str
    delta: str
    shift_search_bits: int
    local_search_bits: int

    @property
    def support_numerator(self) -> int:
        return self.m

    @property
    def support_denominator(self) -> int:
        return self.P


@dataclass(frozen=True)
class BinaryMatrix:
    """Rows are strictly increasing lists of nonzero column indices."""

    columns: int
    rows: tuple[tuple[int, ...], ...]

    @property
    def nrows(self) -> int:
        return len(self.rows)

    @property
    def max_row_weight(self) -> int:
        return max((len(row) for row in self.rows), default=0)

    def accepts(self, mask: int) -> bool:
        for row in self.rows:
            parity = 0
            for column in row:
                parity ^= (mask >> column) & 1
            if parity:
                return False
        return True


@dataclass(frozen=True)
class InitialMatrixResult:
    status: str
    parameters: Lemma61Parameters | None
    matrix: BinaryMatrix | None
    reason: str
    rows_allocated: bool = False
    search_started: bool = False
    note: str = STEP19_SCOPE

    @property
    def constructed(self) -> bool:
        return self.status == 'constructed'

    @property
    def refused(self) -> bool:
        return self.status == 'refused'


@dataclass(frozen=True)
class KernelCheckResult:
    status: str
    ok: bool
    count: None
    reason: str
    vectors_scanned: int = 0
    nonzero_kernel_vectors: int = 0
    min_nonzero_weight: int | None = None
    threshold: str = ''
    note: str = STEP19_SCOPE

    @property
    def refused(self) -> bool:
        return self.status == 'refused'


@dataclass(frozen=True)
class GreedySelectionResult:
    status: str
    count: int | None
    reason: str
    rounds_completed: int = 0
    list_length: int = 0
    candidates_per_round: int = 0
    counts_performed: int = 0
    candidates_enumerated: int = 0
    boolean_vectors_scanned: int = 0
    equations: tuple[SparsePolynomial, ...] | None = None
    matrix_allocated: bool = False
    note: str = STEP19_SCOPE

    @property
    def counted(self) -> bool:
        return self.status == 'counted'

    @property
    def refused(self) -> bool:
        return self.status == 'refused'


@dataclass(frozen=True)
class BranchChoice:
    width: int
    shared: frozenset[Literal]
    petals: tuple[frozenset[Literal], ...]

    @property
    def gap(self) -> int:
        return self.width - len(self.shared)


@dataclass(frozen=True)
class SparsificationNode:
    clauses: tuple[ClauseSet, ...]
    branch: BranchChoice | None
    children: tuple['SparsificationNode', 'SparsificationNode'] | None


@dataclass(frozen=True)
class SparsificationTree:
    status: str
    n: int
    k: int
    M: int
    occurrence_bound: int
    nodes: int
    depth: int
    leaves: tuple[CNF, ...]
    root: SparsificationNode | None
    reason: str
    nodes_allocated: int = 0
    note: str = STEP19_SCOPE

    @property
    def built(self) -> bool:
        return self.status == 'built'

    @property
    def refused(self) -> bool:
        return self.status == 'refused'


@dataclass(frozen=True)
class DisjointCountResult:
    status: str
    count: int | None
    reason: str
    n: int
    leaves: int = 0
    nodes: int = 0
    depth: int = 0
    occurrence_bound: int = 0
    max_occurrence: int | None = None
    same_declared_n: bool = False
    boolean_vectors_scanned: int = 0
    leaf_counts: tuple[int, ...] = ()
    tree: SparsificationTree | None = None
    note: str = STEP19_SCOPE

    @property
    def counted(self) -> bool:
        return self.status == 'counted'

    @property
    def refused(self) -> bool:
        return self.status == 'refused'


@dataclass(frozen=True)
class Proposition64Schedule:
    """Integer schedule from Proposition 6.4. The rho bound is not claimed."""

    n: int
    lambda_n: int
    nu_n: int
    s: int
    sigma: str | None
    exhaustive_search: bool
    rho_bound_claimed: bool = False
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False
    note: str = STEP19_SCOPE


@dataclass(frozen=True)
class BranchingChoice:
    status: str
    k: int
    sigma: str
    M: int | None
    C: int | None
    gamma: int
    h: int | None
    reason: str
    powers_tested: int = 0
    note: str = STEP19_SCOPE

    @property
    def chosen(self) -> bool:
        return self.status == 'chosen'

    @property
    def refused(self) -> bool:
        return self.status == 'refused'


def _index(value: int, name: str) -> int:
    number = operator.index(value)
    if number < 0:
        raise ValueError(f'{name} must be a nonnegative integer')
    return number


def _positive(value: int, name: str) -> int:
    number = operator.index(value)
    if number < 1:
        raise ValueError(f'{name} must be a positive integer')
    return number


def ceil_log2(value: int) -> int:
    """``ceil(log2(value))`` for a positive integer. ``ceil(log2(1))`` is 0."""
    number = _positive(value, 'log argument')
    if number == 1:
        return 0
    return (number - 1).bit_length()


def cube_exceeds(variables: int, limit: int) -> bool:
    """``2^variables > limit`` without building the cube."""
    width = _index(variables, 'n')
    cap = _index(limit, 'limit')
    if cap == 0:
        return True
    return width > cap.bit_length() - 1


def _pow_within(base: int, exponent: int, limit: int) -> int | None:
    """``base**exponent`` when it is ``<= limit``; otherwise None, before the value is stored."""
    if exponent < 0 or base < 0 or limit < 0:
        raise ValueError('power bounds must be nonnegative')
    if exponent == 0:
        return 1 if 1 <= limit else None
    if base == 0:
        return 0 if limit >= 0 else None
    result = 1
    for _ in range(exponent):
        if result > limit // base:
            return None
        result *= base
    return result


def _toggle(present: set[frozenset[int]], monomial: frozenset[int]) -> None:
    if monomial in present:
        present.remove(monomial)
    else:
        present.add(monomial)


def _polynomial_from(n: int, present: set[frozenset[int]]) -> SparsePolynomial:
    return sparse_polynomial(n, [tuple(sorted(monomial)) for monomial in present])


def literal_key(literal: Literal) -> tuple[int, int]:
    return literal[0], int(literal[1])


def clause_key(clause: ClauseSet) -> tuple[int, tuple[tuple[int, int], ...]]:
    ordered = tuple(sorted(clause, key=literal_key))
    return len(ordered), ordered


def clause_violation_polynomial(n: int, clause: Sequence[Literal]) -> SparsePolynomial:
    """Violation bit of one clause, as square-free monomials over F2.

    Repeated literals are deleted first. Both polarities of one variable make
    a tautology, whose violation polynomial is 0. The empty product is 1.
    """
    variables = _index(n, 'n')
    seen: dict[int, bool] = {}
    for raw_var, raw_neg in clause:
        var = operator.index(raw_var)
        neg = bool(raw_neg)
        if var < 0 or var >= variables:
            raise ValueError(f'variable {var} is outside 0..{variables - 1}')
        if var in seen and seen[var] != neg:
            return sparse_polynomial(variables, [])
        seen[var] = neg
    present: set[frozenset[int]] = {frozenset()}
    for var in sorted(seen):
        updated: set[frozenset[int]] = set()
        if seen[var]:
            for monomial in present:
                updated.add(monomial | {var})
        else:
            for monomial in present:
                _toggle(updated, monomial)
                _toggle(updated, monomial | {var})
        present = updated
    return _polynomial_from(variables, present)


def expand_binary_row(
    n: int,
    violations: Sequence[SparsePolynomial],
    pairs: Sequence[tuple[int, int]],
) -> SparsePolynomial:
    """Expand ``sum c_e v_e`` over F2. Repeated clause indices cancel modulo 2.

    A coefficient of 0 contributes nothing. The list length counts every pair,
    including those cancelled pairs; this function only builds the polynomial.
    """
    variables = _index(n, 'n')
    if any(item.n != variables for item in violations):
        raise ValueError('every violation polynomial must use the same n')
    coefficients = [0] * len(violations)
    for raw_index, raw_coeff in pairs:
        index = operator.index(raw_index)
        if index < 0 or index >= len(violations):
            raise ValueError(f'clause index {index} is outside 0..{len(violations) - 1}')
        coefficients[index] = (coefficients[index] + operator.index(raw_coeff)) & 1
    present: set[frozenset[int]] = set()
    for index, bit in enumerate(coefficients):
        if not bit:
            continue
        for monomial in violations[index].monomials:
            _toggle(present, monomial)
    return _polynomial_from(variables, present)


def equation_is_zero(equation: SparsePolynomial, mask: int) -> bool:
    """Evaluate one sparse F2 equation on the bits of ``mask``."""
    parity = 0
    for monomial in equation.monomials:
        alive = 1
        for variable in monomial:
            if (mask >> variable) & 1 == 0:
                alive = 0
                break
        parity ^= alive
    return parity == 0


def count_equations_f2(
    equations: Sequence[SparsePolynomial],
    *,
    max_boolean_vectors: int,
) -> int | None:
    """Independent cube scan. Returns None before the scan when the cube is over budget."""
    if not equations:
        raise ValueError('the cube counter needs n from at least one polynomial')
    variables = equations[0].n
    if any(equation.n != variables for equation in equations):
        raise ValueError('equations must share n')
    if cube_exceeds(variables, max_boolean_vectors):
        return None
    total = 0
    for mask in range(1 << variables):
        if all(equation_is_zero(equation, mask) for equation in equations):
            total += 1
    return total


def violations_system(cnf: CNF, *, k: int | None = None, c: int = 1) -> SparseSystem:
    """One violation equation per clause. An empty clause is the constant 1."""
    polynomials = tuple(clause_violation_polynomial(cnf.n, clause) for clause in cnf.clauses)
    degree = max((polynomial.degree for polynomial in polynomials), default=0)
    width = degree if k is None else _index(k, 'k')
    if width < degree:
        raise ValueError(f'k={width} is below the violation degree {degree}')
    power = _positive(c, 'c')
    return SparseSystem(cnf.n, width, power, polynomials)


def count_violations_common_zeros(
    cnf: CNF,
    budgets: Section5Budgets | None = None,
) -> int | None:
    """Leaf or system count through Step 17. None means that call refused; no partial count."""
    result = count_common_zeros(violations_system(cnf), budgets)
    if result.refused or result.count is None:
        return None
    return result.count


def lemma_6_1_parameters(n: int, m: int, B: int) -> Lemma61Parameters:
    """``μ = floor(log2 m)+1``, ``P = 20000 B^2 μ^2``, ``D0 = 320000 B^2 μ^3``."""
    variables = _positive(n, 'n')
    clauses = _positive(m, 'm')
    occurrence = _positive(B, 'B')
    if clauses > occurrence * variables:
        raise ValueError(f'm={clauses} exceeds B n = {occurrence * variables}')
    mu = (clauses.bit_length() - 1) + 1
    P = 20000 * occurrence * occurrence * mu * mu
    D0 = 320000 * occurrence * occurrence * mu * mu * mu
    d_G = min(clauses, D0)
    need = (clauses + d_G - 1) // d_G
    ell_G = 0 if need <= 1 else (need - 1).bit_length()
    q_G = 1 << ell_G
    if not (clauses <= q_G * d_G < 2 * clauses):
        raise RuntimeError('q_G failed the Lemma 6.1 edge-count window')
    c_G = (d_G + 8 * occurrence - 1) // (8 * occurrence)
    denominator = 100 * occurrence * mu
    minimum_weight = (d_G + denominator - 1) // denominator
    h = minimum_weight - 1
    shift_bits = 0 if q_G == 1 else d_G * ell_G
    local_bits = 0 if h <= 0 else c_G * d_G
    return Lemma61Parameters(
        n=variables,
        m=clauses,
        B=occurrence,
        mu=mu,
        P=P,
        D0=D0,
        d_G=d_G,
        q_G=q_G,
        ell_G=ell_G,
        c_G=c_G,
        h=h,
        epsilon=str(Fraction(1, 200 * occurrence * mu)),
        delta=str(Fraction(1, denominator)),
        shift_search_bits=shift_bits,
        local_search_bits=local_bits,
    )


def _explicit_local_rows(c_G: int, d_G: int) -> tuple[tuple[int, ...], ...]:
    """Even-parity local code. Valid whenever the required distance is at most 1.

    ``h <= 0`` means ``ceil(δ d_G) <= 1``. Every nonzero word then has weight
    at least ``δ d_G``. The first row is the all-ones parity check when
    ``d_G >= 1``; the remaining rows are zero. The kernel is the even-weight
    words, whose nonzero weight is at least 2 when ``d_G >= 2``, and is empty
    of nonzero words when ``d_G = 1``.
    """
    if c_G < 1:
        raise RuntimeError('c_G is positive for every positive d_G in Lemma 6.1')
    rows = [tuple(range(d_G))]
    rows.extend(() for _ in range(c_G - 1))
    return tuple(rows)


def _assemble_regular_checks(
    parameters: Lemma61Parameters,
    local_rows: Sequence[Sequence[int]],
) -> BinaryMatrix:
    """Tanner checks for ``q_G = 1``. Both sides see the ``m`` parallel edges.

    The single left vertex and the single right vertex each receive every
    local row on the port order ``0..d_G-1``. There are no extra edges to
    shorten. Rows may repeat; the construction does not delete them.
    """
    if parameters.q_G != 1 or parameters.d_G != parameters.m:
        raise RuntimeError('the explicit assembly is only the q_G = 1 case')
    rows: list[tuple[int, ...]] = []
    for _side in range(2):
        for local in local_rows:
            support = tuple(column for column in local if 0 <= column < parameters.m)
            rows.append(support)
    if len(rows) != 2 * parameters.q_G * parameters.c_G:
        raise RuntimeError('Tanner row count disagrees with 2 q_G c_G')
    matrix = BinaryMatrix(parameters.m, tuple(rows))
    if matrix.max_row_weight > parameters.D0:
        raise RuntimeError('a constructed row is heavier than D0')
    if matrix.nrows > 3 * parameters.n:
        raise RuntimeError('constructed row count exceeds 3n')
    return matrix


def construct_initial_matrix(
    n: int,
    m: int,
    B: int,
    budgets: Section6Budgets | None = None,
) -> InitialMatrixResult:
    """Deterministic ``A(0)`` when the exhaustive searches fit; otherwise refuse.

    The parameter integers are computed either way. No shift tuple, local
    matrix, or kernel vector is allocated on the refusal path.
    """
    limits = budgets or Section6Budgets()
    parameters = lemma_6_1_parameters(n, m, B)
    if parameters.shift_search_bits > limits.max_shift_search_bits:
        return InitialMatrixResult(
            'refused', parameters, None,
            reason=(
                f'shift search uses {parameters.shift_search_bits} bits, above budget '
                f'{limits.max_shift_search_bits}; no shifts allocated'
            ),
        )
    if parameters.local_search_bits > limits.max_local_search_bits:
        return InitialMatrixResult(
            'refused', parameters, None,
            reason=(
                f'local-code search uses {parameters.local_search_bits} bits, above budget '
                f'{limits.max_local_search_bits}; no local matrix allocated'
            ),
        )
    if parameters.q_G != 1:
        return InitialMatrixResult(
            'refused', parameters, None,
            reason=(
                'q_G > 1; the shift tuples are not allocated'
            ),
        )
    if parameters.h > 0:
        found = search_local_code(
            parameters.d_G,
            parameters.c_G,
            parameters.h + 1,
            max_bits=limits.max_local_search_bits,
            max_words=limits.max_codeword_tests,
        )
        if found is None:
            return InitialMatrixResult(
                'refused', parameters, None,
                reason=(
                    'local-code search did not return a matrix; no Tanner rows allocated'
                ),
                search_started=False,
            )
        matrix = _assemble_regular_checks(parameters, found)
        return InitialMatrixResult(
            'constructed', parameters, matrix,
            reason='exhaustive local-code search produced a matrix of distance at least ceil(δ d_G)',
            rows_allocated=True,
            search_started=True,
        )
    local_rows = _explicit_local_rows(parameters.c_G, parameters.d_G)
    matrix = _assemble_regular_checks(parameters, local_rows)
    return InitialMatrixResult(
        'constructed', parameters, matrix,
        reason=(
            'q_G = 1 and h <= 0, so every local matrix meets the distance bound; '
            'A(0) is the two copies of the even-parity checks'
        ),
        rows_allocated=True,
    )


def search_local_code(
    d_G: int,
    c_G: int,
    min_weight: int,
    *,
    max_bits: int,
    max_words: int,
) -> tuple[tuple[int, ...], ...] | None:
    """First binary ``c_G × d_G`` matrix whose kernel has no nonzero word lighter than ``min_weight``.

    Returns None before allocating the matrices or the test words when either
    budget is too small, and None after a finished search that finds nothing.
    The integer order of matrices is the deterministic enumeration order.
    """
    columns = _index(d_G, 'd_G')
    rows = _index(c_G, 'c_G')
    distance = _index(min_weight, 'min_weight')
    bits = rows * columns
    if bits > max_bits or rows > max_bits:
        return None
    if distance <= 1:
        return _explicit_local_rows(max(rows, 1), columns)[:rows] if rows else ()
    cap = _index(max_words, 'max_words')
    if _codeword_count(columns, distance - 1, cap) is None:
        return None
    words: list[int] = []
    for weight in range(1, distance):
        for indexes in combinations(range(columns), weight):
            mask = 0
            for column in indexes:
                mask |= 1 << column
            words.append(mask)
            if len(words) > cap:
                return None
    total = 1 << bits
    for code in range(total):
        if _local_code_accepts(code, columns, rows, words):
            decoded = []
            for row in range(rows):
                decoded.append(tuple(
                    column for column in range(columns)
                    if (code >> (row * columns + column)) & 1
                ))
            return tuple(decoded)
    return None


def _local_code_accepts(code: int, columns: int, rows: int, words: Sequence[int]) -> bool:
    for word in words:
        in_kernel = True
        for row in range(rows):
            parity = 0
            bit = word
            column = 0
            while bit:
                if bit & 1 and (code >> (row * columns + column)) & 1:
                    parity ^= 1
                bit >>= 1
                column += 1
            if parity:
                in_kernel = False
                break
        if in_kernel:
            return False
    return True


def _codeword_count(d_G: int, h: int, limit: int) -> int | None:
    if h <= 0:
        return 0
    total = 0
    term = 1
    for j in range(1, h + 1):
        term = term * (d_G - j + 1) // j
        total += term
        if total > limit:
            return None
    return total


def kernel_support_check(
    matrix: BinaryMatrix,
    *,
    numerator: int,
    denominator: int,
    max_boolean_vectors: int,
) -> KernelCheckResult:
    """Enumerate ``A z = 0``. Refuse before the loop when ``2^m`` exceeds the budget.

    A nonzero accepted vector must have weight at least ``numerator/denominator``.
    The zero vector is not a violation. No partial list of weights is returned
    on a refusal.
    """
    clauses = _positive(matrix.columns, 'm') if matrix.columns else 0
    num = _positive(numerator, 'numerator')
    den = _positive(denominator, 'denominator')
    threshold = f'{num}/{den}'
    if matrix.columns == 0 or cube_exceeds(matrix.columns, max_boolean_vectors):
        return KernelCheckResult(
            'refused', False, None,
            reason=(
                f'kernel cube 2^{matrix.columns} exceeds budget {max_boolean_vectors}; '
                'no vectors scanned'
            ),
            threshold=threshold,
        )
    scanned = 0
    accepted = 0
    minimum: int | None = None
    failed = False
    for mask in range(1 << clauses):
        scanned += 1
        if mask == 0 or not matrix.accepts(mask):
            continue
        weight = mask.bit_count()
        accepted += 1
        if minimum is None or weight < minimum:
            minimum = weight
        if weight * den < num:
            failed = True
    if failed or minimum is None and accepted:
        return KernelCheckResult(
            'checked', False, None,
            reason='a nonzero kernel vector is lighter than m/P',
            vectors_scanned=scanned,
            nonzero_kernel_vectors=accepted,
            min_nonzero_weight=minimum,
            threshold=threshold,
        )
    reason = (
        'every nonzero vector accepted by A(0) has weight at least m/P'
        if accepted
        else 'the kernel is {0}, so the nonzero-support claim holds vacuously'
    )
    return KernelCheckResult(
        'checked', True, None, reason,
        vectors_scanned=scanned,
        nonzero_kernel_vectors=accepted,
        min_nonzero_weight=minimum,
        threshold=threshold,
    )


def _matrix_equations(cnf: CNF, matrix: BinaryMatrix) -> tuple[SparsePolynomial, ...]:
    violations = [clause_violation_polynomial(cnf.n, clause) for clause in cnf.clauses]
    if matrix.columns != len(violations):
        raise ValueError('A(0) columns must equal the number of clauses')
    return tuple(
        expand_binary_row(cnf.n, violations, [(column, 1) for column in row])
        for row in matrix.rows
    )


def _candidate_pairs(m: int, list_length: int):
    choices = [(index, bit) for index in range(m) for bit in (0, 1)]
    return product(choices, repeat=list_length)


def greedy_select_equations(
    cnf: CNF,
    *,
    list_length: int,
    B: int | None = None,
    rounds: int | None = None,
    budgets: Section6Budgets | None = None,
    counter: SystemCounter | None = None,
    initial_matrix: BinaryMatrix | None = None,
) -> GreedySelectionResult:
    """Append ``rounds`` rows of the given list length, each minimizing ``N``.

    Ties keep the first candidate in product order, clause index slow and
    coefficient ``0`` before ``1``. The paper's round count is ``3n``. The
    candidate budget is applied before ``A(0)`` is built and before any
    counter call. A refusal returns ``count is None``.
    """
    limits = budgets or Section6Budgets()
    variables = _index(cnf.n, 'n')
    clauses = _index(cnf.m, 'm')
    length = _positive(list_length, 'list length')
    round_count = 3 * variables if rounds is None else _index(rounds, 'rounds')
    occurrence = max(1, (clauses + variables - 1) // variables) if B is None and clauses else B
    candidate_count = _pow_within(2 * clauses, length, limits.max_candidate_lists)
    if clauses == 0:
        return GreedySelectionResult(
            'refused', None,
            reason='Lemma 6.2 requires at least one clause; no candidates enumerated',
        )
    if occurrence is None:
        return GreedySelectionResult('refused', None, reason='B is required when m > 0')
    if candidate_count is None:
        return GreedySelectionResult(
            'refused', None,
            reason=(
                f'(2m)^P = ({2 * clauses})^{length} exceeds candidate budget '
                f'{limits.max_candidate_lists}; no candidates enumerated'
            ),
            list_length=length,
        )
    if round_count > limits.max_rounds:
        return GreedySelectionResult(
            'refused', None,
            reason=(
                f'{round_count} rounds exceed budget {limits.max_rounds}; '
                'no candidates enumerated'
            ),
            list_length=length,
            candidates_per_round=candidate_count,
        )
    if round_count > 0 and candidate_count > limits.max_counts // round_count:
        return GreedySelectionResult(
            'refused', None,
            reason=(
                f'{round_count * candidate_count} counts exceed budget {limits.max_counts}; '
                'no candidates enumerated'
            ),
            list_length=length,
            candidates_per_round=candidate_count,
        )
    if counter is None and cube_exceeds(variables, limits.max_boolean_vectors):
        return GreedySelectionResult(
            'refused', None,
            reason=(
                f'Boolean cube 2^{variables} exceeds budget {limits.max_boolean_vectors}; '
                'no candidates enumerated and no assignments scanned'
            ),
            list_length=length,
            candidates_per_round=candidate_count,
        )
    if initial_matrix is None:
        built = construct_initial_matrix(variables, clauses, _positive(occurrence, 'B'), limits)
        if built.refused or built.matrix is None:
            return GreedySelectionResult(
                'refused', None,
                reason=built.reason,
                list_length=length,
                candidates_per_round=candidate_count,
                matrix_allocated=False,
            )
        matrix = built.matrix
    else:
        matrix = initial_matrix
    violations = [clause_violation_polynomial(variables, clause) for clause in cnf.clauses]
    equations = list(_matrix_equations(cnf, matrix))
    default_counter = counter is None
    if counter is None:
        cap = limits.max_boolean_vectors

        def counter(system: tuple[SparsePolynomial, ...]) -> int | None:
            return count_equations_f2(system, max_boolean_vectors=cap)

    counts_performed = 0
    vectors = 0
    cube = 0 if cube_exceeds(variables, limits.max_boolean_vectors) else 1 << variables
    current_count: int | None = None
    for _round in range(round_count):
        best_count: int | None = None
        best_row: SparsePolynomial | None = None
        for gamma in _candidate_pairs(clauses, length):
            row = expand_binary_row(variables, violations, gamma)
            try:
                value = counter(tuple(equations) + (row,))
            except Section6Refusal:
                value = None
            counts_performed += 1
            if value is None:
                return GreedySelectionResult(
                    'refused', None,
                    reason='the injected counter refused; no partial count returned',
                    rounds_completed=0,
                    list_length=length,
                    candidates_per_round=candidate_count,
                    counts_performed=0,
                    candidates_enumerated=0,
                    boolean_vectors_scanned=0,
                    matrix_allocated=True,
                )
            if best_count is None or value < best_count:
                best_count = value
                best_row = row
        if best_row is None or best_count is None:
            return GreedySelectionResult(
                'refused', None,
                reason='the candidate list was empty; no count returned',
                matrix_allocated=True,
            )
        equations.append(best_row)
        current_count = best_count
        if default_counter:
            vectors += cube * candidate_count
    if current_count is None:
        return GreedySelectionResult(
            'refused', None,
            reason='no round was run; no count returned',
            list_length=length,
            candidates_per_round=candidate_count,
            matrix_allocated=True,
        )
    return GreedySelectionResult(
        'counted', current_count,
        reason='minimum-N rows appended in enumeration order',
        rounds_completed=round_count,
        list_length=length,
        candidates_per_round=candidate_count,
        counts_performed=counts_performed,
        candidates_enumerated=counts_performed,
        boolean_vectors_scanned=vectors,
        equations=tuple(equations),
        matrix_allocated=True,
    )


def normalize_clauses(n: int, clauses: Sequence[Sequence[Literal]]) -> tuple[ClauseSet, ...]:
    """Delete repeated literals, tautologies, duplicates and containing clauses.

    All ``n`` declared variables stay. Unit clauses are not propagated.
    """
    variables = _index(n, 'n')
    cleaned: list[ClauseSet] = []
    for clause in clauses:
        seen: dict[int, bool] = {}
        tautology = False
        for raw_var, raw_neg in clause:
            var = operator.index(raw_var)
            neg = bool(raw_neg)
            if var < 0 or var >= variables:
                raise ValueError(f'variable {var} is outside 0..{variables - 1}')
            if var in seen and seen[var] != neg:
                tautology = True
                break
            seen[var] = neg
        if tautology:
            continue
        cleaned.append(frozenset((var, neg) for var, neg in seen.items()))
    unique = sorted(set(cleaned), key=clause_key)
    kept = [clause for clause in unique if not any(other < clause for other in unique)]
    return tuple(kept)


def normalize_cnf(cnf: CNF) -> CNF:
    return _clauses_to_cnf(cnf.n, normalize_clauses(cnf.n, cnf.clauses))


def _clauses_to_cnf(n: int, clauses: Sequence[ClauseSet]) -> CNF:
    rendered = [
        sorted(clause, key=literal_key)
        for clause in sorted(clauses, key=clause_key)
    ]
    return CNF(n, rendered)


def lemma_6_3_occurrence_bound(k: int, M: int) -> int:
    """``B = 2 sum_{j=0}^{k-1} M^j`` from Appendix A."""
    width = _positive(k, 'k')
    if width < 2:
        raise ValueError('k must be at least 2')
    branching = _positive(M, 'M')
    if branching == 1:
        total = width
    else:
        total = (branching ** width - 1) // (branching - 1)
    return 2 * total


def variable_occurrences(cnf: CNF) -> tuple[int, ...]:
    counts = [0] * cnf.n
    for clause in cnf.clauses:
        seen: set[int] = set()
        for var, _neg in clause:
            if var not in seen:
                counts[var] += 1
                seen.add(var)
    return tuple(counts)


def choose_branching_M(
    k: int,
    sigma: Fraction,
    budgets: Section6Budgets | None = None,
) -> BranchingChoice:
    """Smallest power of two ``M >= Γ`` with ``Γ h <= σ M``.

    ``h = ceil(log2(4 C M))`` and ``C = 4 sum_{j<k-1} (3M)^j``. The search
    stops before the next power when that integer would exceed the bit budget.
    This does not build the tree and does not claim a node-count runtime.
    """
    limits = budgets or Section6Budgets()
    width = _positive(k, 'k')
    if width < 2:
        raise ValueError('k must be at least 2')
    ratio = Fraction(sigma)
    if ratio <= 0 or ratio > 1:
        raise ValueError('sigma must satisfy 0 < sigma <= 1')
    gamma = 2 * (3 ** (width - 1) - 1)
    M = 1
    while M < gamma:
        if (M * 2).bit_length() > limits.max_branching_bits:
            return BranchingChoice(
                'refused', width, str(ratio), None, None, gamma, None,
                reason=(
                    f'the next power of two exceeds {limits.max_branching_bits} bits; '
                    'M was not stored'
                ),
            )
        M *= 2
    tested = 0
    while True:
        tested += 1
        if M.bit_length() > limits.max_branching_bits:
            return BranchingChoice(
                'refused', width, str(ratio), None, None, gamma, None,
                reason=(
                    f'M exceeds {limits.max_branching_bits} bits before Γ h <= σ M; '
                    'the tree was not built'
                ),
                powers_tested=tested,
            )
        C = 4 if width == 2 else 4 * ((3 * M) ** (width - 1) - 1) // (3 * M - 1)
        product = 4 * C * M
        h = 0 if product <= 1 else (product - 1).bit_length()
        if gamma * h * ratio.denominator <= ratio.numerator * M:
            return BranchingChoice(
                'chosen', width, str(ratio), M, C, gamma, h,
                reason='first power of two with Γ h <= σ M',
                powers_tested=tested,
            )
        if (M * 2).bit_length() > limits.max_branching_bits:
            return BranchingChoice(
                'refused', width, str(ratio), None, None, gamma, None,
                reason=(
                    f'the next power of two exceeds {limits.max_branching_bits} bits; '
                    'the tree was not built'
                ),
                powers_tested=tested,
            )
        M *= 2


def _choose_branch(clauses: Sequence[ClauseSet], k: int, M: int) -> BranchChoice | None:
    grouped: dict[int, list[ClauseSet]] = {}
    for clause in clauses:
        width = len(clause)
        if width > k:
            raise ValueError(f'clause width {width} exceeds k={k}')
        grouped.setdefault(width, []).append(clause)
    best_key: tuple[int, int, tuple[tuple[int, int], ...]] | None = None
    best_H: ClauseSet | None = None
    best_c = 0
    for c in range(2, k + 1):
        group = grouped.get(c, [])
        if len(group) < 2 and M > 1:
            # A single clause cannot meet θ_j = M^j >= M >= 2.
            continue
        counts: dict[ClauseSet, int] = {}
        for clause in group:
            ordered = tuple(sorted(clause, key=literal_key))
            for size in range(1, c):
                for combo in combinations(ordered, size):
                    shared = frozenset(combo)
                    counts[shared] = counts.get(shared, 0) + 1
        for shared, degree in counts.items():
            gap = c - len(shared)
            if degree < M ** gap:
                continue
            key = (c, -len(shared), tuple(sorted(shared, key=literal_key)))
            if best_key is None or key < best_key:
                best_key = key
                best_H = shared
                best_c = c
    if best_H is None:
        return None
    selected = tuple(
        sorted(
            (clause for clause in grouped[best_c] if best_H <= clause),
            key=clause_key,
        )
    )
    petals = tuple(clause - best_H for clause in selected)
    return BranchChoice(best_c, best_H, petals)


def _split_children(
    clauses: Sequence[ClauseSet],
    branch: BranchChoice,
) -> tuple[tuple[ClauseSet, ...], tuple[ClauseSet, ...]]:
    selected = {branch.shared | petal for petal in branch.petals}
    others = [clause for clause in clauses if clause not in selected]
    left = _renormalize(others + [branch.shared])
    units = [frozenset(((var, not neg),)) for var, neg in sorted(branch.shared, key=literal_key)]
    right = _renormalize([*others, *branch.petals, *units])
    return left, right


def _renormalize(clauses: Sequence[ClauseSet]) -> tuple[ClauseSet, ...]:
    unique = sorted({clause for clause in clauses}, key=clause_key)
    return tuple(clause for clause in unique if not any(other < clause for other in unique))


class _TreeBudget(Exception):
    """Internal marker. The partial tree is discarded by the caller."""


def build_disjoint_sparsification(
    cnf: CNF,
    *,
    k: int,
    M: int,
    budgets: Section6Budgets | None = None,
) -> SparsificationTree:
    """Appendix A tree. Children keep ``n``. A refusal returns no leaves to sum."""
    limits = budgets or Section6Budgets()
    width = _positive(k, 'k')
    branching = _positive(M, 'M')
    if width < 2 or branching < 2:
        raise ValueError('Lemma 6.3 requires k >= 2 and M >= 2')
    bound = lemma_6_3_occurrence_bound(width, branching)
    root_clauses = normalize_clauses(cnf.n, cnf.clauses)
    state = {'nodes': 0, 'depth': 0, 'leaves': []}

    def rec(clauses: tuple[ClauseSet, ...], depth: int) -> SparsificationNode:
        state['nodes'] += 1
        if state['nodes'] > limits.max_tree_nodes:
            raise _TreeBudget
        branch = _choose_branch(clauses, width, branching)
        if branch is None:
            state['depth'] = max(state['depth'], depth)
            state['leaves'].append(_clauses_to_cnf(cnf.n, clauses))
            return SparsificationNode(clauses, None, None)
        if len(branch.petals) < branching ** branch.gap:
            raise RuntimeError('the selected shared set does not meet θ_j')
        if any(len(petal) != branch.gap for petal in branch.petals):
            raise RuntimeError('a petal does not have width c-|H|')
        left_clauses, right_clauses = _split_children(clauses, branch)
        if branch.shared not in left_clauses:
            raise RuntimeError('the new A-clause did not survive normalization')
        for petal in branch.petals:
            if petal not in right_clauses:
                raise RuntimeError('a new B-clause did not survive normalization')
        for var, neg in branch.shared:
            if frozenset(((var, not neg),)) not in right_clauses:
                raise RuntimeError('a new unit clause did not survive normalization')
        left = rec(left_clauses, depth + 1)
        right = rec(right_clauses, depth + 1)
        return SparsificationNode(clauses, branch, (left, right))

    try:
        root = rec(root_clauses, 0)
    except _TreeBudget:
        return SparsificationTree(
            'refused', cnf.n, width, branching, bound, 0, 0, (), None,
            reason=(
                f'tree node budget {limits.max_tree_nodes} exceeded; '
                'no leaves returned and no assignments scanned'
            ),
            nodes_allocated=0,
        )
    except ValueError as exc:
        raise ValueError(str(exc)) from exc
    leaves = tuple(state['leaves'])
    if any(leaf.n != cnf.n for leaf in leaves):
        raise RuntimeError('a leaf changed the declared n')
    return SparsificationTree(
        'built', cnf.n, width, branching, bound,
        state['nodes'], state['depth'], leaves, root,
        reason='children partition the parent satisfying assignments and keep n',
        nodes_allocated=state['nodes'],
    )


def count_via_disjoint_sparsification(
    cnf: CNF,
    *,
    k: int,
    M: int,
    leaf_counter: LeafCounter | None = None,
    budgets: Section6Budgets | None = None,
    verify_disjoint: bool = False,
) -> DisjointCountResult:
    """Sum leaf counts. Empty clause is 0 and the empty formula is ``2^n``, without a scan.

    The cube budget is checked before the tree is built when the default
    counter is used. A refusal returns ``count is None`` and an empty leaf list.
    """
    limits = budgets or Section6Budgets()
    normalized = normalize_cnf(cnf)
    if any(len(clause) == 0 for clause in normalized.clauses):
        return DisjointCountResult(
            'counted', 0,
            reason='an empty clause is unsatisfiable; no tree and no scan',
            n=cnf.n,
            leaves=1,
            nodes=1,
            same_declared_n=True,
            max_occurrence=max(variable_occurrences(normalized), default=0),
        )
    if normalized.m == 0:
        return DisjointCountResult(
            'counted', 1 << cnf.n,
            reason='a formula with no clauses has 2^n assignments; no tree and no scan',
            n=cnf.n,
            leaves=1,
            nodes=1,
            same_declared_n=True,
            max_occurrence=0,
        )
    if leaf_counter is None and cube_exceeds(cnf.n, limits.max_boolean_vectors):
        return DisjointCountResult(
            'refused', None,
            reason=(
                f'Boolean cube 2^{cnf.n} exceeds budget {limits.max_boolean_vectors}; '
                'the tree was not built and no assignments were scanned'
            ),
            n=cnf.n,
        )
    tree = build_disjoint_sparsification(cnf, k=k, M=M, budgets=limits)
    if tree.refused:
        return DisjointCountResult(
            'refused', None, tree.reason, cnf.n,
            occurrence_bound=tree.occurrence_bound,
        )
    counter = leaf_counter
    if counter is None:
        cap = limits.max_boolean_vectors

        def counter(leaf: CNF) -> int | None:
            if cube_exceeds(leaf.n, cap):
                return None
            if leaf.n > 20:
                return None
            return brute_count(leaf)

    counts: list[int] = []
    scanned = 0
    cube = 0 if cube_exceeds(cnf.n, limits.max_boolean_vectors) else 1 << cnf.n
    for leaf in tree.leaves:
        if leaf.n != cnf.n:
            return DisjointCountResult(
                'refused', None,
                reason='a leaf changed n; no partial count returned',
                n=cnf.n,
            )
        try:
            value = counter(leaf)
        except Section6Refusal:
            value = None
        if value is None:
            return DisjointCountResult(
                'refused', None,
                reason='a leaf counter refused; no partial count returned',
                n=cnf.n,
                occurrence_bound=tree.occurrence_bound,
            )
        counts.append(value)
        scanned += cube
    occurrences = [max(variable_occurrences(leaf), default=0) for leaf in tree.leaves]
    maximum = max(occurrences, default=0)
    if maximum > tree.occurrence_bound:
        return DisjointCountResult(
            'refused', None,
            reason=(
                f'leaf occurrence {maximum} exceeds the computed bound {tree.occurrence_bound}; '
                'no count returned'
            ),
            n=cnf.n,
            occurrence_bound=tree.occurrence_bound,
            max_occurrence=maximum,
        )
    disjoint_note = ''
    if verify_disjoint:
        if cube_exceeds(cnf.n, limits.max_boolean_vectors):
            return DisjointCountResult(
                'refused', None,
                reason='disjointness scan exceeds the cube budget; no count returned',
                n=cnf.n,
            )
        parent = set(_satisfying_masks(normalized))
        seen: set[int] = set()
        for leaf in tree.leaves:
            models = set(_satisfying_masks(leaf))
            if seen & models:
                return DisjointCountResult(
                    'refused', None,
                    reason='two leaves share a satisfying assignment; no count returned',
                    n=cnf.n,
                )
            seen |= models
        if seen != parent:
            return DisjointCountResult(
                'refused', None,
                reason='the leaves do not partition the normalized satisfying assignments',
                n=cnf.n,
            )
        scanned += cube * (len(tree.leaves) + 1)
        disjoint_note = ' The satisfying assignments partition across the leaves.'
    return DisjointCountResult(
        'counted', sum(counts),
        reason='sum of injected leaf counts on one declared n.' + disjoint_note,
        n=cnf.n,
        leaves=len(tree.leaves),
        nodes=tree.nodes,
        depth=tree.depth,
        occurrence_bound=tree.occurrence_bound,
        max_occurrence=maximum,
        same_declared_n=True,
        boolean_vectors_scanned=scanned,
        leaf_counts=tuple(counts),
        tree=tree,
    )


def _satisfying_masks(cnf: CNF) -> list[int]:
    masks = []
    for mask in range(1 << cnf.n):
        good = True
        for clause in cnf.clauses:
            if not clause:
                good = False
                break
            satisfied = False
            for var, neg in clause:
                bit = (mask >> var) & 1
                if (bit == 0) if neg else (bit == 1):
                    satisfied = True
                    break
            if not satisfied:
                good = False
                break
        if good:
            masks.append(mask)
    return masks


def largest_s(lambda_n: int, nu_n: int) -> int:
    """Largest nonnegative integer ``s`` with ``s^2 ν^4 <= λ``.

    ``ν = 0`` makes the displayed floor ``sqrt(λ)/ν^2`` undefined. The
    exhaustive-search case applies, so the returned ``s`` is 0.
    """
    lam = _index(lambda_n, 'lambda')
    nu = _index(nu_n, 'nu')
    if nu == 0:
        return 0
    nu4 = nu ** 4
    if lam < nu4:
        return 0
    lo = 0
    hi = lam
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid * mid * nu4 <= lam:
            lo = mid
        else:
            hi = mid - 1
    return lo


def proposition_6_4_schedule(n: int) -> Proposition64Schedule:
    """``λ_n``, ``ν_n``, largest ``s`` with ``s^2 ν_n^4 <= λ_n``, and ``σ=1/s`` when ``s≥2``.

    ``s < 2`` is the exhaustive-search case in the proof. No rho bound is claimed.
    """
    variables = _positive(n, 'n')
    lam = ceil_log2(variables)
    nu = 0 if lam <= 1 else ceil_log2(lam)
    s = largest_s(lam, nu)
    exhaustive = s < 2
    sigma = None if exhaustive else str(Fraction(1, s))
    return Proposition64Schedule(
        n=variables,
        lambda_n=lam,
        nu_n=nu,
        s=s,
        sigma=sigma,
        exhaustive_search=exhaustive,
    )
