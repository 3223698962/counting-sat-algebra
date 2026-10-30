"""Finite counting extension rings and Hermite placement from SqrtLog/RingCore.

``RingCore`` builds ``O = (ZMod 3^t)[X] / (F)`` where ``F`` is a monic lift of an
irreducible degree-``d`` polynomial over ``F_3``. ``IrreducibleSearch`` finds it
by complete enumeration: the first passing monic candidate of the ``3^d``
coefficient vectors, with an exhaustive factorization test. Nodes are the
fixed, distinct, nonzero residues with radix-3 index ``1..Q``; their pairwise
differences are units. The Hermite guard is ``width*R*(2*factors +
3*(t-1) + 1) < Q*J``: division plus one gives strictly more interpolation
conditions than the degree, and a natural count below ``3^t`` reads back
uniquely from coefficient zero. The scheduled ``Q`` is zero until ``L`` is in
the millions, so the finite witnesses use explicit ``(t, d, Q, J)`` that
satisfy the relation set, exactly like the explicit candidate depth of
Step 25. The streamed Taylor evaluation of the full leaf and the Python
bit-time/space transfer remain open.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import product
from math import isqrt

from .lean_code import radix_words
from .lean_high import early_log
from .lean_ring import natural

STEP27_SCOPE = (
    'Finite counting extension rings (searched irreducible lifts, residue units, '
    'unique natural readback), fixed radix nodes and the Q*J > degree Hermite '
    'condition with exact Taylor jets and reconstruction. The full streamed '
    'Taylor leaf, placement into the native ABI and the Python bit-time/space '
    'transfer remain open. No Theorem 1.1, no 2^{o(n)} claim.'
)


class ExtensionRefusal(ValueError):
    """An explicit limit stops the construction before completion."""


def asympt_w(ambient: int) -> int:
    """``Asymptotic.w``: floor(log2 L) + 1 for the schedule level ``L``."""
    return early_log(ambient).bit_length()


def asympt_s(ambient: int) -> int:
    """``Asymptotic.s``: isqrt(L) / w^2, the square-root-log schedule level."""
    return isqrt(early_log(ambient)) // asympt_w(ambient) ** 2


def layout_r(ambient: int) -> int:
    return max(2, asympt_s(ambient))


def layout_den(ambient: int) -> int:
    return 64 * (natural(layout_r(ambient), 'layout_r').bit_length())


def layout_q(ambient: int) -> int:
    """``Asymptotic.layoutQ``: the schedule's interpolation node count."""
    return early_log(ambient) // layout_den(ambient)


def layout_energy(ambient: int, occurrence: int) -> int:
    occurrence = natural(occurrence, 'occurrence')
    return (64 * occurrence * layout_r(ambient) * early_log(ambient) ** 2).bit_length()


def layout_cutoff(ambient: int, occurrence: int) -> int:
    """``Asymptotic.layoutCutoff``: the shared-layout radius R."""
    return isqrt(128 * layout_q(ambient) * layout_energy(ambient, occurrence)) + 4


@dataclass(frozen=True)
class ExtensionSchedule:
    ambient: int
    free_count: int
    clauses: int
    width: int
    factors: int
    occurrence: int
    arity: int
    incidences: int
    level: int
    w: int
    s: int
    r: int
    Q: int
    R: int
    precision: int
    degree: int
    nodes_squared: int
    jet_order: int

    @property
    def degree_bound(self) -> int:
        return self.Q * self.jet_order

    @property
    def guarded(self) -> bool:
        return (
            0 < self.free_count <= self.ambient
            and 0 < self.clauses
            and 0 < self.occurrence
            and 0 < self.arity
            and self.free_count <= self.incidences
            and 2 <= self.Q
            and 4 * self.R <= self.Q
            and 16 * self.Q * self.r**2 <= self.free_count
        )

    @property
    def hermite_conditions_exceed_degree(self) -> bool:
        return (
            self.width * self.R * (2 * self.factors + 3 * (self.precision - 1) + 1)
            < self.degree_bound
        )

    @property
    def field_cardinality_cap(self) -> bool:
        return 3**self.degree <= 16 * (self.Q + 1) ** 2


def extension_schedule(
    ambient: int,
    *,
    free_count: int | None = None,
    clauses: int,
    width: int,
    factors: int,
    occurrence: int,
    arity: int = 1,
    incidences: int | None = None,
) -> ExtensionSchedule:
    """The exact ``LeafParametersProgram.produce`` integers on a declared input."""
    ambient = natural(ambient, 'ambient')
    free_count = natural(ambient if free_count is None else free_count, 'free_count')
    clauses = natural(clauses, 'clauses')
    width = natural(width, 'width')
    factors = natural(factors, 'factors')
    occurrence = natural(occurrence, 'occurrence')
    arity = natural(arity, 'arity')
    incidences = natural(free_count if incidences is None else incidences, 'incidences')
    level = early_log(free_count)
    r = layout_r(free_count)
    q = layout_q(free_count)
    precision = ambient + natural(clauses, 'clauses').bit_length() + 1
    degree = natural(q, 'q').bit_length() + 1 if q else 1
    jet = width * layout_cutoff(free_count, occurrence) * (
        2 * factors + 3 * (precision - 1) + 1
    ) // q + 1 if q else 0
    return ExtensionSchedule(
        ambient=ambient,
        free_count=free_count,
        clauses=clauses,
        width=width,
        factors=factors,
        occurrence=occurrence,
        arity=arity,
        incidences=incidences,
        level=level,
        w=asympt_w(free_count),
        s=asympt_s(free_count),
        r=r,
        Q=q,
        R=layout_cutoff(free_count, occurrence),
        precision=precision,
        degree=degree,
        nodes_squared=r ** (2 * q) if q else 0,
        jet_order=jet,
    )


@dataclass(frozen=True)
class ExtensionBudgets:
    max_degree: int = 6
    max_candidates: int = 100_000
    max_division_checks: int = 1_000_000
    max_ring_operations: int = 2_000_000

    def __post_init__(self) -> None:
        for item in fields(self):
            object.__setattr__(self, item.name, natural(getattr(self, item.name), item.name))


class _Meter:
    def __init__(self, budgets: ExtensionBudgets):
        self.budgets = budgets
        self.candidates = 0
        self.division_checks = 0
        self.operations = 0

    def candidate(self) -> None:
        self.candidates += 1
        if self.candidates > self.budgets.max_candidates:
            raise ExtensionRefusal('candidate budget exceeded')

    def division(self) -> None:
        self.division_checks += 1
        if self.division_checks > self.budgets.max_division_checks:
            raise ExtensionRefusal('factorization test budget exceeded')

    def work(self, count: int) -> None:
        self.operations += count
        if self.operations > self.budgets.max_ring_operations:
            raise ExtensionRefusal('ring operation budget exceeded')


def _f3_divmod(numerator: tuple[int, ...], denominator: tuple[int, ...], meter: _Meter):
    """Monic division over F_3; the denominator must be monic and nonconstant."""
    remainder = list(numerator)
    quotient = [0] * max(0, len(numerator) - len(denominator) + 1)
    for k in range(len(remainder) - 1, len(denominator) - 2, -1):
        coefficient = remainder[k] % 3
        quotient[k - len(denominator) + 1] = coefficient
        meter.division()
        meter.work(len(denominator))
        for i, d in enumerate(denominator):
            remainder[k - len(denominator) + 1 + i] = (
                remainder[k - len(denominator) + 1 + i] - coefficient * d
            ) % 3
    while remainder and remainder[-1] == 0:
        remainder.pop()
    return tuple(quotient), tuple(remainder)


def is_irreducible_f3(coefficients: tuple[int, ...], meter: _Meter | None = None) -> bool:
    """Exhaustive factorization test for a monic polynomial over F_3."""
    used = meter or _Meter(ExtensionBudgets())
    degree = len(coefficients)
    poly = tuple(coefficients) + (1,)
    for factor_degree in range(1, degree):
        for candidate in product((0, 1, 2), repeat=factor_degree):
            divisor = tuple(candidate) + (1,)
            _quotient, remainder = _f3_divmod(poly, divisor, used)
            if not remainder:
                return False
    return True


def search_irreducible_polynomial(
    degree: int, *, budgets: ExtensionBudgets | None = None
) -> tuple[tuple[int, ...] | None, int, int]:
    """``IrreducibleSearch.search``: first passing monic candidate, in odometer order."""
    degree = natural(degree, 'degree')
    limits = budgets or ExtensionBudgets()
    if not degree:
        raise ValueError('the extension degree must be positive')
    if degree > limits.max_degree:
        raise ExtensionRefusal('extension degree budget exceeded')
    meter = _Meter(limits)
    for digits in radix_words(3, degree):
        meter.candidate()
        coefficients = tuple(digits)
        if is_irreducible_f3(coefficients, meter):
            return coefficients, meter.candidates, meter.division_checks
    return None, meter.candidates, meter.division_checks


@dataclass(frozen=True)
class CountingExtension:
    precision: int
    degree: int
    coefficients: tuple[int, ...]

    def __post_init__(self) -> None:
        precision, degree = natural(self.precision, 'precision'), natural(
            self.degree, 'degree'
        )
        if not precision:
            raise ValueError('the precision must be positive')
        if not degree:
            raise ValueError('the extension degree must be positive')
        coefficients = tuple(int(c) % 3 for c in self.coefficients)
        if len(coefficients) != degree:
            raise ValueError('the defining polynomial must have exactly degree coefficients')
        object.__setattr__(self, 'precision', precision)
        object.__setattr__(self, 'degree', degree)
        object.__setattr__(self, 'coefficients', coefficients)

    @property
    def modulus(self) -> int:
        return 3**self.precision

    def element(self, coefficients) -> tuple[int, ...]:
        values = [int(c) % self.modulus for c in coefficients]
        if len(values) > self.degree:
            raise ValueError('an element has at most d coefficients')
        return tuple(values + [0] * (self.degree - len(values)))

    zero = property(lambda self: (0,) * self.degree)
    one = property(lambda self: (1,) + (0,) * (self.degree - 1))

    def constant(self, value: int) -> tuple[int, ...]:
        return self.element([value])

    def add(self, left, right):
        return tuple((a + b) % self.modulus for a, b in zip(left, right))

    def sub(self, left, right):
        return tuple((a - b) % self.modulus for a, b in zip(left, right))

    def neg(self, value):
        return tuple((-a) % self.modulus for a in value)

    def scale(self, scalar: int, value):
        scalar %= self.modulus
        return tuple(scalar * a % self.modulus for a in value)

    def mul(self, left, right):
        d, modulus = self.degree, self.modulus
        raw = [0] * (2 * d - 1)
        for i, a in enumerate(left):
            if not a:
                continue
            for j, b in enumerate(right):
                raw[i + j] = (raw[i + j] + a * b) % modulus
        for k in range(len(raw) - 1, d - 1, -1):
            factor = raw[k]
            if factor:
                raw[k] = 0
                for i, c in enumerate(self.coefficients):
                    raw[k - d + i] = (raw[k - d + i] - factor * c) % modulus
        return tuple(raw[:d])

    def pow(self, value, exponent: int):
        exponent = natural(exponent, 'exponent')
        result, base = self.one, value
        while exponent:
            if exponent & 1:
                result = self.mul(result, base)
            exponent >>= 1
            if exponent:
                base = self.mul(base, base)
        return result

    def residue(self, value):
        return tuple(a % 3 for a in value)

    def residue_mul(self, left, right):
        d = self.degree
        raw = [0] * (2 * d - 1)
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                raw[i + j] = (raw[i + j] + a * b) % 3
        for k in range(len(raw) - 1, d - 1, -1):
            factor = raw[k]
            if factor:
                raw[k] = 0
                for i, c in enumerate(self.coefficients):
                    raw[k - d + i] = (raw[k - d + i] - factor * c) % 3
        return tuple(raw[:d])

    def residue_inverse(self, value):
        """Fermat inversion in the residue field, exponent ``3^d - 2``."""
        if not any(self.residue(value)):
            raise ValueError('the residue is zero and has no inverse')
        return self.residue_pow(self.residue(value), 3**self.degree - 2)

    def residue_pow(self, value, exponent: int):
        result, base = (1,) + (0,) * (self.degree - 1), value
        while exponent:
            if exponent & 1:
                result = self.residue_mul(result, base)
            exponent >>= 1
            if exponent:
                base = self.residue_mul(base, base)
        return result

    def is_unit(self, value) -> bool:
        """In the local ring exactly the nonzero residues are units."""
        return any(self.residue(value))

    def inverse(self, value):
        """Newton-lift the residue inverse to the full ``3^t`` precision."""
        if not self.is_unit(value):
            raise ValueError('a non-unit has no inverse')
        lifted = self.element(list(self.residue_inverse(value)))
        for _ in range(max(0, self.precision - 1)):
            two = self.constant(2)
            lifted = self.mul(lifted, self.sub(two, self.mul(value, lifted)))
        if self.mul(value, lifted) != self.one:
            raise ArithmeticError('the Newton inversion failed to reach precision')
        return lifted

    def nat_cast(self, value: int) -> tuple[int, ...]:
        """The constant element ``(value : O)``."""
        return self.constant(natural(value, 'value'))

    def decode_count(self, value) -> int:
        """``decodeCount``: the constant coefficient as a natural, below 3^t."""
        return int(value[0]) % self.modulus


def counting_extension(
    precision: int, degree: int, *, budgets: ExtensionBudgets | None = None
) -> CountingExtension:
    """Search the first irreducible monic lift; no supplied polynomial is trusted."""
    precision = natural(precision, 'precision')
    degree = natural(degree, 'degree')
    if not precision:
        raise ValueError('the precision must be positive')
    coefficients, _candidates, _checks = search_irreducible_polynomial(
        degree, budgets=budgets
    )
    if coefficients is None:
        raise ExtensionRefusal('the scheduled irreducible search found no witness')
    return CountingExtension(precision, degree, coefficients)


def node_element(extension: CountingExtension, index: int) -> tuple[int, ...]:
    """``nodeCode`` node ``index`` (1-based): digit j is ``(index / 3^j) % 3``."""
    index = natural(index, 'node index')
    if not index or index >= 3**extension.degree:
        raise ValueError('the node index must lie in 1..3^d - 1')
    digits = []
    remaining = index
    for _ in range(extension.degree):
        digits.append(remaining % 3)
        remaining //= 3
    return extension.element(digits)


def node_code_words(extension: CountingExtension, index: int) -> tuple[int, ...]:
    """The physical 2t-bit little-endian words of one node's coefficient code."""
    node = node_element(extension, index)
    return tuple(int(digit) % (1 << (2 * extension.precision)) for digit in node)


def node_stream(precision: int, degree: int, nodes: int) -> tuple[tuple[tuple[int, ...], ...], int]:
    """``RingNodeProgram.stream``: nodes 1..Q in ascending order, flat words."""
    coefficients, _candidates, _checks = search_irreducible_polynomial(degree)
    if coefficients is None:
        raise ExtensionRefusal('the scheduled irreducible search found no witness')
    extension = CountingExtension(precision, degree, coefficients)
    nodes = natural(nodes, 'nodes')
    if nodes >= 3**degree:
        raise ValueError('the node count must stay below 3^d')
    stream = tuple(node_code_words(extension, index) for index in range(1, nodes + 1))
    bits = nodes * degree * 2 * precision
    return stream, bits


def taylor_jets(extension: CountingExtension, node, target, jet_order: int):
    """The ``J`` coefficients of ``target(X + node)``, by integer binomials."""
    jet_order = natural(jet_order, 'jet order')
    modulus = extension.modulus
    jets = []
    for k in range(jet_order):
        coefficient = extension.zero
        for j in range(k, len(target)):
            binomial = 1
            for i in range(k):
                binomial = binomial * (j - i) // (i + 1)
            if binomial % modulus:
                term = extension.scale(binomial, extension.mul(target[j], extension.pow(node, j - k)))
                coefficient = extension.add(coefficient, term)
        jets.append(coefficient)
    return tuple(jets)


def hermite_reconstruct(extension: CountingExtension, nodes, jets, jet_order: int):
    """Newton-Hermite from J-jets at distinct nodes; divisions are unit differences."""
    jet_order = natural(jet_order, 'jet order')
    nodes = list(nodes)
    if len(jets) != len(nodes):
        raise ValueError('every node needs one jet block')
    points = []
    for node, block in zip(nodes, jets):
        if len(block) != jet_order:
            raise ValueError('every jet block must have J coefficients')
        points.extend([node] * jet_order)
    size = len(points)
    table = [[None] * size for _ in range(size)]
    # The diagonal entry at every point of a block is the value f(node);
    # confluent entries table[i][i+k] inside one block are the k-th jet.
    for i in range(size):
        table[i][i] = jets[i // jet_order][0]
    for k in range(1, size):
        for i in range(size - k):
            j = i + k
            if points[i] == points[j]:
                table[i][j] = jets[i // jet_order][k]
            else:
                difference = extension.sub(points[j], points[i])
                inverse = extension.inverse(difference)
                numerator = extension.sub(table[i + 1][j], table[i][j - 1])
                table[i][j] = extension.mul(numerator, inverse)
    polynomial = [extension.zero]
    basis = [extension.one]
    for k in range(size):
        if len(polynomial) < len(basis):
            polynomial.extend([extension.zero] * (len(basis) - len(polynomial)))
        for i, coefficient in enumerate(basis):
            polynomial[i] = extension.add(
                polynomial[i], extension.mul(table[0][k], coefficient)
            )
        basis = _mul_by_linear(extension, basis, points[k])
    while polynomial and not any(polynomial[-1]):
        polynomial.pop()
    return tuple(polynomial)


def _mul_by_linear(extension: CountingExtension, polynomial, root):
    output = [extension.zero] * (len(polynomial) + 1)
    neg_root = extension.neg(root)
    for i, coefficient in enumerate(polynomial):
        output[i + 1] = extension.add(output[i + 1], coefficient)
        output[i] = extension.add(output[i], extension.mul(coefficient, neg_root))
    return output
