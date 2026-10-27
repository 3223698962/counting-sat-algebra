"""Finite Section 5.1–5.2 objects for sparse systems over F2.

Each equation is a short list of square-free monomials. Duplicate monomials
cancel modulo two, and the empty monomial is the constant 1. The integer lift
of an equation is even on a Boolean point exactly when the equation is zero
over F2. Lemma 5.1 builds the complemented Beigel–Tarui polynomial and divides
it by the monic quadratic h(Y) = Y(Y-1). The truncated product keeps only the
terms whose total number of h factors is less than ell.

The Boolean sum below is an exhaustive scan of a budgeted cube. It checks the
residue identity on tiny systems. It does not build replacement forms, select
Section 6 equations, or run Appendix D. It is not Theorem 2.2's algorithm and
it is not Theorem 1.1. ``theorem_1_1_complete`` stays false.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass, fields
from typing import Sequence

STEP17_SCOPE = (
    'Finite Section 5.1–5.2 checks: sparse F2 monomials, integer lifts, '
    'Lemma 5.1 coefficient division, and the truncated product in Z[V]/(V^ell). '
    'Common-zero checks on tiny cubes are exhaustive. Section 5.3 replacement '
    'forms, Section 6 equation selection, and Appendix D are not run. '
    'This does not implement Theorem 1.1 and it is not a 2^{o(n)} algorithm.'
)


class BudgetRefusal(ValueError):
    """A Section 5 budget stopped the work. No partial count is available."""


@dataclass(frozen=True)
class Section5Budgets:
    """Explicit caps. Oversized ell, coefficient bits, term tables, and cubes stop first."""

    max_ell: int = 24
    max_coefficient_bits: int = 200_000
    max_terms: int = 4096
    max_equations: int = 32
    max_boolean_vectors: int = 4096
    max_accumulator_bits: int = 200_000

    def __post_init__(self):
        for field in fields(self):
            value = operator.index(getattr(self, field.name))
            if value < 0:
                raise ValueError('Section 5 budgets must be nonnegative integers')
            object.__setattr__(self, field.name, value)


@dataclass(frozen=True)
class SparsePolynomial:
    """Square-free monomials over F2, already combined modulo two."""

    n: int
    monomials: tuple[frozenset[int], ...]

    @property
    def degree(self) -> int:
        return max((len(monomial) for monomial in self.monomials), default=0)

    @property
    def support(self) -> frozenset[int]:
        if not self.monomials:
            return frozenset()
        return frozenset().union(*self.monomials)

    @property
    def is_zero(self) -> bool:
        return not self.monomials

    @property
    def is_one(self) -> bool:
        return len(self.monomials) == 1 and not self.monomials[0]


@dataclass(frozen=True)
class SparseSystem:
    """A system P1=...=Pm_test=0 with the Theorem 2.2 sparsity bounds.

    ``k >= 2``, ``c >= 1`` and ``n >= 20`` are the theorem's hypotheses. Smaller
    finite systems are accepted when their monomials still meet degree ``k`` and
    the integer allowance ``floor((log2 n)^c)``. Meeting those bounds on a tiny
    cube is not the asymptotic statement.
    """

    n: int
    k: int
    c: int
    equations: tuple[SparsePolynomial, ...]

    @property
    def m_test(self) -> int:
        return len(self.equations)

    @property
    def meets_theorem_size(self) -> bool:
        return self.k >= 2 and self.c >= 1 and self.n >= 20


@dataclass(frozen=True)
class ParityPolynomials:
    """``e_j(Y) = a_j + b_j Y`` for ``0 <= j < ell``, from Lemma 5.1."""

    ell: int
    coefficients: tuple[tuple[int, int], ...]
    coefficient_bits: int

    def __post_init__(self):
        if len(self.coefficients) != self.ell:
            raise ValueError('Lemma 5.1 returns one linear polynomial per power below ell')

    def evaluate(self, y: int) -> int:
        """Exact integer value of ``p_ell(y) = sum_j e_j(y) h(y)^j``."""
        value = operator.index(y)
        height = value * (value - 1)
        total = 0
        height_power = 1
        for slope_constant, slope in self.coefficients:
            total += (slope_constant + slope * value) * height_power
            height_power *= height
        return total

    def indicator_agrees(self, y: int) -> bool:
        """``p_ell(y) ≡ 1{2 divides y} (mod 2^ell)``."""
        value = operator.index(y)
        expected = 0 if value % 2 else 1
        return self.evaluate(value) % (1 << self.ell) == expected


@dataclass(frozen=True)
class RetainedTerm:
    """One product term. ``powers[i]`` is the number of ``h(L_i)`` factors."""

    powers: tuple[int, ...]

    @property
    def quadratic_factors(self) -> int:
        return sum(self.powers)

    @property
    def linear_factors(self) -> int:
        return len(self.powers)

    @property
    def lift_degree(self) -> int:
        """Degree in the lifts: one linear ``e_j`` and ``2 j_i`` from each ``h^{j_i}``."""
        return self.linear_factors + 2 * self.quadratic_factors


@dataclass(frozen=True)
class TruncatedProduct:
    """Terms of ``prod_i p_ell(L_i)`` whose total h-power is strictly less than ell."""

    ell: int
    m_test: int
    terms: tuple[RetainedTerm, ...]

    @property
    def linear_factor_bound(self) -> int:
        return self.m_test

    @property
    def quadratic_factor_bound(self) -> int:
        return self.ell - 1

    def shape_ok(self) -> bool:
        quadratic_bound = self.quadratic_factor_bound
        linear_bound = self.linear_factor_bound
        degree_bound = linear_bound + 2 * quadratic_bound
        seen: set[tuple[int, ...]] = set()
        for term in self.terms:
            powers = term.powers
            if len(powers) != self.m_test or powers in seen:
                return False
            if any(power < 0 for power in powers):
                return False
            if term.quadratic_factors > quadratic_bound or term.linear_factors > linear_bound:
                return False
            if term.lift_degree > degree_bound:
                return False
            seen.add(powers)
        return True


@dataclass(frozen=True)
class CommonZeroResult:
    status: str
    count: int | None
    reason: str
    ell: int
    boolean_vectors: int = 0
    equations_used: int = 0
    equations_discarded: int = 0
    terms: int = 0
    term_table_allocated: bool = False
    coefficients_constructed: bool = False
    coefficient_bits: int = 0
    linear_factor_bound: int = 0
    quadratic_factor_bound: int = 0
    input_degree_bound: int = 0
    enumerated: bool = False
    note: str = STEP17_SCOPE

    @property
    def counted(self) -> bool:
        return self.status == 'counted'

    @property
    def refused(self) -> bool:
        return self.status == 'refused'


def monomial_allowance(n: int, c: int) -> int:
    """Largest integer ``m`` with ``m <= (log2 n)^c``.

    The comparison is exact integer arithmetic. For ``n < 2`` the logarithm is
    not positive, so the allowance is 0. This is the integer form of the
    monomial bound in Theorem 2.2; it does not assert the running-time claim.
    """
    variables = operator.index(n)
    power = operator.index(c)
    if power < 1:
        raise ValueError('c must be a positive integer')
    if variables < 2:
        return 0
    floor_log = variables.bit_length() - 1
    upper = (floor_log + 1) ** power
    lo = 0
    hi = upper
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if _log_power_fits(variables, power, mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


def sparse_polynomial(
    n: int,
    monomials: Sequence[Sequence[int]],
    *,
    degree_bound: int | None = None,
    monomial_bound: int | None = None,
) -> SparsePolynomial:
    """Combine duplicate monomials modulo two.

    A repeated index inside one monomial is the same Boolean variable, so it is
    stored once. The empty monomial is the constant 1. The empty list is the
    zero polynomial.
    """
    variables = operator.index(n)
    if variables < 0:
        raise ValueError('n must be nonnegative')
    if degree_bound is not None:
        degree_bound = operator.index(degree_bound)
        if degree_bound < 0:
            raise ValueError('degree bound must be nonnegative')
    if monomial_bound is not None:
        monomial_bound = operator.index(monomial_bound)
        if monomial_bound < 0:
            raise ValueError('monomial bound must be nonnegative')
        if len(monomials) > monomial_bound:
            raise ValueError(
                f'{len(monomials)} monomials were supplied; the bound is {monomial_bound}'
            )
    present: set[frozenset[int]] = set()
    for raw in monomials:
        indices = []
        for item in raw:
            index = operator.index(item)
            if index < 0 or index >= variables:
                raise ValueError(f'variable index {index} is outside 0..{variables - 1}')
            indices.append(index)
        monomial = frozenset(indices)
        if degree_bound is not None and len(monomial) > degree_bound:
            raise ValueError(f'monomial degree {len(monomial)} exceeds k={degree_bound}')
        if monomial in present:
            present.remove(monomial)
        else:
            present.add(monomial)
    if monomial_bound is not None and len(present) > monomial_bound:
        raise ValueError(
            f'{len(present)} monomials remain after cancellation; the bound is {monomial_bound}'
        )
    ordered = tuple(sorted(present, key=lambda item: (len(item), tuple(sorted(item)))))
    return SparsePolynomial(variables, ordered)


def sparse_system(
    n: int,
    equations: Sequence[Sequence[Sequence[int]]],
    *,
    k: int,
    c: int,
) -> SparseSystem:
    """Build a system whose equations meet degree ``k`` and the ``(log2 n)^c`` bound."""
    variables = operator.index(n)
    degree = operator.index(k)
    power = operator.index(c)
    if variables < 0:
        raise ValueError('n must be nonnegative')
    if degree < 0:
        raise ValueError('k must be nonnegative')
    if power < 1:
        raise ValueError('c must be a positive integer')
    allowance = monomial_allowance(variables, power)
    built = []
    for index, monomials in enumerate(equations):
        try:
            polynomial = sparse_polynomial(
                variables, monomials, degree_bound=degree, monomial_bound=allowance,
            )
        except ValueError as exc:
            raise ValueError(f'equation {index}: {exc}') from exc
        built.append(polynomial)
    return SparseSystem(variables, degree, power, tuple(built))


def integer_lift(polynomial: SparsePolynomial, assignment: Sequence[int]) -> int:
    """Evaluate ``L(X) = sum_M prod_{v in M} X_v``.

    On a Boolean assignment this counts the monomials equal to one, so the F2
    equation is zero exactly when the lift is even.
    """
    if len(assignment) != polynomial.n:
        raise ValueError('assignment length must equal n')
    values = [operator.index(bit) for bit in assignment]
    total = 0
    for monomial in polynomial.monomials:
        factor = 1
        for variable in monomial:
            factor *= values[variable]
            if factor == 0:
                break
        total += factor
    return total


def beigel_tarui_complement(ell: int, max_bits: int | None = None) -> tuple[int, ...]:
    """Coefficients of ``p_ell(Y) = (1-Y)^ell * sum_{j<ell} C(ell+j-1, j) Y^j``.

    Index 0 is the constant term. The coefficients are exact integers. A bit
    budget stops the construction before the next coefficient is stored.
    """
    width = operator.index(ell)
    if width < 1:
        raise ValueError('ell must be a positive integer')
    series = [1]
    for j in range(1, width):
        numerator = series[-1] * (width + j - 1)
        if numerator % j:
            raise RuntimeError('Beigel–Tarui binomial coefficient was not an integer')
        coefficient = numerator // j
        _require_bits(coefficient, max_bits, 'Beigel–Tarui coefficient')
        series.append(coefficient)
    factor = [1]
    for _ in range(width):
        factor = _multiply(factor, (1, -1))
        for coefficient in factor:
            _require_bits(coefficient, max_bits, 'Beigel–Tarui factor')
    product = _multiply(factor, series)
    for coefficient in product:
        _require_bits(coefficient, max_bits, 'p_ell coefficient')
    product = _trim(product)
    if len(product) != 2 * width or product[0] != 1 or any(product[j] for j in range(1, width)):
        raise RuntimeError('p_ell is not congruent to 1 modulo Y^ell')
    return tuple(product)


def construct_parity_polynomials(
    ell: int, budgets: Section5Budgets | None = None,
) -> ParityPolynomials:
    """Divide ``p_ell`` by ``h(Y) = Y^2 - Y`` until every remainder has degree at most 1.

    The returned list is complete or the call raises ``BudgetRefusal``. A refusal
    does not return the coefficients constructed so far.
    """
    limits = budgets or Section5Budgets()
    width = operator.index(ell)
    if width < 1:
        raise ValueError('ell must be a positive integer')
    if width > limits.max_ell:
        raise BudgetRefusal(
            f'ell={width} exceeds budget {limits.max_ell}; parity coefficients not constructed'
        )
    raw = beigel_tarui_complement(width, limits.max_coefficient_bits)
    coefficients = _successive_division(raw, width, limits.max_coefficient_bits)
    expanded = _expand_divided(coefficients)
    if tuple(_trim(expanded)) != raw:
        raise RuntimeError('successive division did not reproduce p_ell')
    bits = max(max(abs(constant).bit_length(), abs(slope).bit_length()) for constant, slope in coefficients)
    return ParityPolynomials(width, tuple(coefficients), bits)


def count_retained_terms(m_test: int, ell: int, limit: int) -> int | None:
    """Number of nonnegative solutions of ``j_1 + ... + j_m < ell``, or None if above ``limit``.

    The binomial coefficient is computed without allocating the tuples. The
    running value stays at most ``limit``.
    """
    equations = operator.index(m_test)
    width = operator.index(ell)
    cap = operator.index(limit)
    if equations < 0 or width < 1 or cap < 0:
        raise ValueError('term count expects m_test >= 0, ell >= 1 and a nonnegative limit')
    # Solutions of j_1 + ... + j_m + slack = ell - 1.
    return _bounded_binomial(equations + width - 1, equations, cap)


def build_truncated_product(
    m_test: int, ell: int, budgets: Section5Budgets | None = None,
) -> TruncatedProduct:
    """Multiply equation factors, keeping a running h-count and dropping totals ``>= ell``.

    The term count is checked before each new list is allocated.
    """
    limits = budgets or Section5Budgets()
    equations = operator.index(m_test)
    width = operator.index(ell)
    if equations < 0:
        raise ValueError('m_test must be nonnegative')
    if width < 1:
        raise ValueError('ell must be a positive integer')
    if width > limits.max_ell:
        raise BudgetRefusal(
            f'ell={width} exceeds budget {limits.max_ell}; term table not allocated'
        )
    if equations > limits.max_equations:
        raise BudgetRefusal(
            f'm_test={equations} exceeds budget {limits.max_equations}; term table not allocated'
        )
    predicted = count_retained_terms(equations, width, limits.max_terms)
    if predicted is None:
        raise BudgetRefusal(
            f'retained term count exceeds budget {limits.max_terms}; term table not allocated'
        )
    terms = [RetainedTerm(())]
    for _ in range(equations):
        next_count = 0
        for term in terms:
            next_count += width - term.quadratic_factors
        if next_count > limits.max_terms:
            raise BudgetRefusal(
                f'retained term count exceeds budget {limits.max_terms}; term table not allocated'
            )
        grown: list[RetainedTerm] = []
        for term in terms:
            room = width - term.quadratic_factors
            powers = term.powers
            for power in range(room):
                grown.append(RetainedTerm(powers + (power,)))
        terms = grown
    if len(terms) != predicted or not TruncatedProduct(width, equations, tuple(terms)).shape_ok():
        raise RuntimeError('truncated term table disagrees with the h-count bound')
    return TruncatedProduct(width, equations, tuple(terms))


def evaluate_truncated_product(
    product: TruncatedProduct,
    parity: ParityPolynomials,
    lifts: Sequence[int],
) -> int:
    """Sum the retained terms at integer lifts. Like terms in ``V`` are not required."""
    if parity.ell != product.ell:
        raise ValueError('the parity polynomials and the product must use the same ell')
    if len(lifts) != product.m_test:
        raise ValueError('one integer lift is required per equation')
    columns = [_factor_column(parity, operator.index(value)) for value in lifts]
    total = 0
    for term in product.terms:
        value = 1
        for column, power in zip(columns, term.powers):
            value *= column[power]
        total += value
    return total


def truncated_indicator_residue(
    parity: ParityPolynomials, lifts: Sequence[int],
) -> int:
    """Least nonnegative residue of ``S_ell(prod B_i)`` modulo ``2^ell``.

    Multiplication is in ``Z[V]/(V^ell)``: the exponent of ``V`` is the running
    h-count, and a product that would reach ``ell`` is dropped. This is the
    combined form of the same truncation. It allocates ``ell`` integers, not the
    multi-index table.
    """
    width = parity.ell
    accumulator = [1] + [0] * (width - 1)
    for lift in lifts:
        column = _factor_column(parity, operator.index(lift))
        combined = [0] * width
        for left, left_value in enumerate(accumulator):
            if left_value == 0:
                continue
            room = width - left
            for right in range(room):
                combined[left + right] += left_value * column[right]
        accumulator = combined
    return sum(accumulator) % (1 << width)


def count_common_zeros(
    system: SparseSystem,
    budgets: Section5Budgets | None = None,
    *,
    ell: int | None = None,
) -> CommonZeroResult:
    """Return ``N(P)`` as the least nonnegative residue of the truncated indicator.

    ``ell`` defaults to ``n + 1``. Constant equations are answered without a term
    table: ``1 = 0`` contributes no points, ``0 = 0`` is discarded, and an empty
    remaining list contributes ``2^n``. A budget failure returns ``count is None``
    and ``boolean_vectors == 0`` before the term table or the cube is built.
    A counted non-constant system has scanned the whole cube. This scan is not
    the subexponential algorithm.
    """
    limits = budgets or Section5Budgets()
    width = system.n + 1 if ell is None else operator.index(ell)
    if width < 1:
        raise ValueError('ell must be a positive integer')
    if width <= system.n:
        raise ValueError('ell must be at least n+1 so the residue equals N(P)')
    active = []
    discarded = 0
    unsatisfiable = False
    for equation in system.equations:
        if equation.is_zero:
            discarded += 1
        elif equation.is_one:
            unsatisfiable = True
        else:
            active.append(equation)
    if unsatisfiable:
        return _counted(
            0, width, used=0, discarded=discarded,
            reason='constant equation 1 = 0 has no Boolean solutions',
        )
    if not active:
        return _counted(
            1 << system.n, width, used=0, discarded=discarded,
            reason='empty equation list has 2^n common zeros',
        )
    cube = 1 << system.n
    if cube > limits.max_boolean_vectors:
        return _refused(
            width, used=len(active), discarded=discarded,
            reason=(
                f'Boolean cube 2^{system.n} exceeds budget {limits.max_boolean_vectors}; '
                'no assignments scanned'
            ),
        )
    if width > limits.max_ell:
        return _refused(
            width, used=len(active), discarded=discarded,
            reason=(
                f'ell={width} exceeds budget {limits.max_ell}; '
                'parity coefficients not constructed'
            ),
        )
    if len(active) > limits.max_equations:
        return _refused(
            width, used=len(active), discarded=discarded,
            reason=(
                f'm_test={len(active)} exceeds budget {limits.max_equations}; '
                'term table not allocated'
            ),
        )
    predicted = count_retained_terms(len(active), width, limits.max_terms)
    if predicted is None:
        return _refused(
            width, used=len(active), discarded=discarded,
            reason=(
                f'retained term count exceeds budget {limits.max_terms}; '
                'term table not allocated'
            ),
        )
    try:
        parity = construct_parity_polynomials(width, limits)
    except BudgetRefusal as exc:
        return _refused(
            width, used=len(active), discarded=discarded, reason=str(exc),
        )
    bound = _accumulator_bit_bound(parity, active, system.n, predicted)
    if bound > limits.max_accumulator_bits:
        return _refused(
            width, used=len(active), discarded=discarded,
            reason=(
                f'accumulator bit bound {bound} exceeds budget {limits.max_accumulator_bits}; '
                'term table not allocated and no assignments scanned'
            ),
            coefficients_constructed=True, coefficient_bits=parity.coefficient_bits,
        )
    try:
        product = build_truncated_product(len(active), width, limits)
    except BudgetRefusal as exc:
        return _refused(
            width, used=len(active), discarded=discarded, reason=str(exc),
            coefficients_constructed=True, coefficient_bits=parity.coefficient_bits,
        )
    total = 0
    for mask in range(cube):
        lifts = [_lift_mask(equation, mask) for equation in active]
        total += evaluate_truncated_product(product, parity, lifts)
    residue = total % (1 << width)
    degree_bound = system.k * (len(active) + 2 * (width - 1))
    return CommonZeroResult(
        status='counted',
        count=residue,
        reason='least nonnegative residue of the truncated indicator on the full cube',
        ell=width,
        boolean_vectors=cube,
        equations_used=len(active),
        equations_discarded=discarded,
        terms=len(product.terms),
        term_table_allocated=True,
        coefficients_constructed=True,
        coefficient_bits=parity.coefficient_bits,
        linear_factor_bound=product.linear_factor_bound,
        quadratic_factor_bound=product.quadratic_factor_bound,
        input_degree_bound=degree_bound,
        enumerated=True,
    )


def brute_common_zeros(system: SparseSystem, *, max_variables: int = 16) -> int:
    """Independent F2 scan used only as a finite oracle.

    This evaluates monomials modulo two. It does not use ``p_ell``. It refuses
    ``n`` above ``max_variables`` instead of scanning a large cube.
    """
    if system.n > max_variables:
        raise ValueError(f'brute-force oracle refuses n={system.n} > {max_variables}')
    total = 0
    equations = system.equations
    for mask in range(1 << system.n):
        if all(_f2_zero(equation, mask) for equation in equations):
            total += 1
    return total


def _counted(count, ell, *, used, discarded, reason) -> CommonZeroResult:
    return CommonZeroResult(
        status='counted', count=count, reason=reason, ell=ell,
        equations_used=used, equations_discarded=discarded,
        enumerated=False, term_table_allocated=False, coefficients_constructed=False,
    )


def _refused(
    ell, *, used, discarded, reason,
    coefficients_constructed=False, coefficient_bits=0,
) -> CommonZeroResult:
    return CommonZeroResult(
        status='refused', count=None, reason=reason, ell=ell,
        boolean_vectors=0, equations_used=used, equations_discarded=discarded,
        terms=0, term_table_allocated=False,
        coefficients_constructed=coefficients_constructed,
        coefficient_bits=coefficient_bits, enumerated=False,
    )


def _require_bits(value: int, max_bits: int | None, label: str) -> None:
    if max_bits is not None and value.bit_length() > max_bits:
        raise BudgetRefusal(
            f'{label} bit length {value.bit_length()} exceeds budget {max_bits}; '
            'parity coefficients not constructed'
        )


def _log_power_fits(n: int, c: int, m: int) -> bool:
    """Return True when ``m <= (log2 n)^c``.

    ``(n**scale).bit_length() - 1`` is ``floor(scale * log2 n)``. Squaring
    tightens that bracket until it sits entirely on one side of ``m``.
    """
    if m <= 0:
        return True
    if n < 2:
        return False
    floor_log = n.bit_length() - 1
    if m <= floor_log ** c:
        return True
    if n == 1 << floor_log:
        return False
    if m >= (floor_log + 1) ** c:
        return False
    powered = n
    scale = 1
    while scale <= 1 << 20:
        floor_scaled = powered.bit_length() - 1
        denominator = scale ** c
        if powered.bit_count() == 1:
            return m * denominator <= floor_scaled ** c
        if m * denominator <= floor_scaled ** c:
            return True
        if m * denominator >= (floor_scaled + 1) ** c:
            return False
        powered *= powered
        scale *= 2
    raise RuntimeError('exact log comparison did not separate')


def _trim(coefficients: Sequence[int]) -> list[int]:
    trimmed = list(coefficients)
    while len(trimmed) > 1 and trimmed[-1] == 0:
        trimmed.pop()
    return trimmed or [0]


def _multiply(left: Sequence[int], right: Sequence[int]) -> list[int]:
    if not left or not right:
        return [0]
    product = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        if left_value == 0:
            continue
        for j, right_value in enumerate(right):
            product[i + j] += left_value * right_value
    return _trim(product)


def _successive_division(
    coefficients: Sequence[int], ell: int, max_bits: int,
) -> tuple[tuple[int, int], ...]:
    current = list(coefficients)
    remainders: list[tuple[int, int]] = []
    for _ in range(ell - 1):
        current, remainder = _divide_once(current, max_bits)
        remainders.append(remainder)
    final = _linear_remainder(current, max_bits)
    remainders.append(final)
    return tuple(remainders)


def _divide_once(coefficients: Sequence[int], max_bits: int) -> tuple[list[int], tuple[int, int]]:
    """One monic division by ``h = Y^2 - Y``. The quotient and remainder stay in ``Z[Y]``."""
    work = list(coefficients)
    if len(work) <= 2:
        return [0], _linear_remainder(work, max_bits)
    quotient = [0] * (len(work) - 2)
    for degree in range(len(work) - 1, 1, -1):
        lead = work[degree]
        if lead == 0:
            continue
        _require_bits(lead, max_bits, 'division coefficient')
        updated = work[degree - 1] + lead
        _require_bits(updated, max_bits, 'division coefficient')
        work[degree] = 0
        work[degree - 1] = updated
        quotient[degree - 2] = lead
    return _trim(quotient), _linear_remainder(work, max_bits)


def _linear_remainder(coefficients: Sequence[int], max_bits: int) -> tuple[int, int]:
    trimmed = _trim(coefficients)
    if len(trimmed) > 2:
        raise RuntimeError('monic division left a remainder of degree greater than 1')
    constant = trimmed[0] if trimmed else 0
    slope = trimmed[1] if len(trimmed) > 1 else 0
    _require_bits(constant, max_bits, 'remainder coefficient')
    _require_bits(slope, max_bits, 'remainder coefficient')
    return constant, slope


def _expand_divided(coefficients: Sequence[tuple[int, int]]) -> list[int]:
    height = (0, -1, 1)
    total = [0]
    height_power = [1]
    for constant, slope in coefficients:
        linear = [constant] if slope == 0 else [constant, slope]
        term = _multiply(height_power, linear)
        width = max(len(total), len(term))
        total = _trim([
            (total[index] if index < len(total) else 0) + (term[index] if index < len(term) else 0)
            for index in range(width)
        ])
        height_power = _multiply(height_power, height)
    return total


def _factor_column(parity: ParityPolynomials, lift: int) -> tuple[int, ...]:
    height = lift * (lift - 1)
    column = []
    height_power = 1
    for constant, slope in parity.coefficients:
        column.append((constant + slope * lift) * height_power)
        height_power *= height
    return tuple(column)


def _lift_mask(polynomial: SparsePolynomial, mask: int) -> int:
    total = 0
    for monomial in polynomial.monomials:
        include = True
        for variable in monomial:
            if (mask >> variable) & 1 == 0:
                include = False
                break
        if include:
            total += 1
    return total


def _f2_zero(polynomial: SparsePolynomial, mask: int) -> bool:
    return _lift_mask(polynomial, mask) % 2 == 0


def _accumulator_bit_bound(
    parity: ParityPolynomials,
    equations: Sequence[SparsePolynomial],
    n: int,
    term_count: int,
) -> int:
    """Upper bound on the bit length of the unreduced cube sum, before any scan."""
    height_max = 1
    linear_bits = 0
    for equation in equations:
        monomials = len(equation.monomials)
        height = monomials * (monomials - 1) if monomials >= 2 else 0
        if height > height_max:
            height_max = height
        largest = 0
        for constant, slope in parity.coefficients:
            size = abs(constant) + abs(slope) * monomials
            if size > largest:
                largest = size
        linear_bits += max(largest.bit_length(), 1)
    height_bits = 0 if parity.ell <= 1 else (parity.ell - 1) * max(height_max.bit_length(), 1)
    return n + term_count.bit_length() + linear_bits + height_bits + 2


def _bounded_binomial(n: int, k: int, limit: int) -> int | None:
    if k < 0 or n < 0 or k > n:
        return 0
    k = min(k, n - k)
    if k == 0:
        return 1 if limit >= 1 else None
    result = 1
    for index in range(k):
        numerator = n - index
        denominator = index + 1
        cap = (limit * denominator) // numerator
        if result > cap:
            return None
        result = result * numerator // denominator
    if result > limit:
        return None
    return result
