"""Finite coefficient algebra from SqrtLogSpace/TruncatedRing.lean.

The shared marker U and the Taylor variable Z have separate cutoffs. Marker
readout is additive and is applied after the complete product. The ternary
indicator follows IdempotentCore, IdempotentPolynomial and IdempotentBase:
Newton lifting modulo (Y(Y-1)(Y-2))**t, followed by monic digit division.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass, field, fields
from typing import Sequence


class LeanAlgebraRefusal(ValueError):
    """An explicit algebra limit prevents returning a completed count."""


def natural(value: int, name: str) -> int:
    value = operator.index(value)
    if value < 0:
        raise ValueError(f'{name} must be nonnegative')
    return value


@dataclass(frozen=True)
class LeanAlgebraBudgets:
    max_variables: int = 100_000
    max_clauses: int = 100_000
    max_literals: int = 1_000_000
    max_query_rows: int = 100_000
    max_row_terms: int = 1_000_000
    max_formula_bits: int = 16_000_000
    max_precision: int = 64
    max_taylor_length: int = 256
    max_dense_slots: int = 65_536
    max_polynomial_terms: int = 16_384
    max_live_coefficients: int = 262_144
    max_algebra_bits: int = 16_000_000
    max_monomial_products: int = 2_000_000
    max_scalar_operations: int = 8_000_000

    def __post_init__(self) -> None:
        for item in fields(self):
            object.__setattr__(self, item.name, natural(getattr(self, item.name), item.name))


@dataclass
class AlgebraMeter:
    limits: LeanAlgebraBudgets
    scalar_operations: int = 0
    monomial_products: int = 0
    peak_live_coefficients: int = 0
    peak_algebra_bits: int = 0
    peak_polynomial_terms: int = 0

    def operations(self, number: int) -> None:
        if self.scalar_operations + number > self.limits.max_scalar_operations:
            raise LeanAlgebraRefusal('scalar operation budget exceeded')
        self.scalar_operations += number

    def storage(self, coefficients: int, bits: int) -> None:
        if coefficients > self.limits.max_live_coefficients:
            raise LeanAlgebraRefusal('live coefficient budget exceeded')
        if bits > self.limits.max_algebra_bits:
            raise LeanAlgebraRefusal('algebra encoding bit budget exceeded')
        self.peak_live_coefficients = max(self.peak_live_coefficients, coefficients)
        self.peak_algebra_bits = max(self.peak_algebra_bits, bits)


@dataclass(frozen=True)
class TruncatedRing:
    """Dense Z/modulus[U,Z]/(U**t,Z**J), in marker-major order.

    Zero-length cutoffs give the zero ring, as in Lean's Fin 0 arrays.
    The storage budget charges the two operands, result and one scalar scratch.
    These logical coefficient charges are not Python heap or tape measurements.
    """

    t: int
    J: int
    modulus: int
    budgets: LeanAlgebraBudgets = field(default_factory=LeanAlgebraBudgets, compare=False)

    def __post_init__(self) -> None:
        for name in ('t', 'J', 'modulus'):
            object.__setattr__(self, name, natural(getattr(self, name), name))
        if self.modulus < 2:
            raise ValueError('modulus must be at least two')
        if self.t > self.budgets.max_precision or self.J > self.budgets.max_taylor_length:
            raise LeanAlgebraRefusal('dense cutoff exceeds its dimension budget')
        if self.slots > self.budgets.max_dense_slots:
            raise LeanAlgebraRefusal('dense coefficient slot budget exceeded')
        coefficients = 3 * self.slots + 1
        if coefficients > self.budgets.max_live_coefficients:
            raise LeanAlgebraRefusal('dense live coefficient budget exceeded')
        if coefficients * (self.modulus - 1).bit_length() > self.budgets.max_algebra_bits:
            raise LeanAlgebraRefusal('dense coefficient bit budget exceeded')

    @property
    def slots(self) -> int:
        return self.t * self.J

    def element(self, coefficients: Sequence[int]) -> DenseElement:
        if len(coefficients) != self.slots:
            raise ValueError('the dense array must retain exactly t*J coefficients')
        return DenseElement(self, tuple(operator.index(c) % self.modulus for c in coefficients))

    def scalar(self, value: int) -> DenseElement:
        if not self.slots:
            return self.element(())
        return self.element((operator.index(value),) + (0,) * (self.slots - 1))

    def monomial(self, u: int, z: int, value: int = 1) -> DenseElement:
        u, z = natural(u, 'u'), natural(z, 'z')
        data = [0] * self.slots
        if u < self.t and z < self.J:
            data[u * self.J + z] = operator.index(value) % self.modulus
        return self.element(data)


@dataclass(frozen=True)
class DenseElement:
    ring: TruncatedRing
    coefficients: tuple[int, ...]

    def __post_init__(self) -> None:
        if len(self.coefficients) != self.ring.slots:
            raise ValueError('the dense array has the wrong length')
        object.__setattr__(
            self,
            'coefficients',
            tuple(operator.index(c) % self.ring.modulus for c in self.coefficients),
        )

    def _same(self, other: DenseElement) -> None:
        if not isinstance(other, DenseElement) or self.ring != other.ring:
            raise ValueError('dense operands must belong to the same truncated ring')

    def __add__(self, other: DenseElement) -> DenseElement:
        self._same(other)
        if self.ring.slots > self.ring.budgets.max_scalar_operations:
            raise LeanAlgebraRefusal('dense addition operation budget exceeded')
        return self.ring.element(
            tuple(a + b for a, b in zip(self.coefficients, other.coefficients))
        )

    def __neg__(self) -> DenseElement:
        if self.ring.slots > self.ring.budgets.max_scalar_operations:
            raise LeanAlgebraRefusal('dense negation operation budget exceeded')
        return self.ring.element(tuple(-c for c in self.coefficients))

    def __sub__(self, other: DenseElement) -> DenseElement:
        self._same(other)
        if self.ring.slots > self.ring.budgets.max_scalar_operations:
            raise LeanAlgebraRefusal('dense subtraction operation budget exceeded')
        return self.ring.element(
            tuple(a - b for a, b in zip(self.coefficients, other.coefficients))
        )

    def __mul__(self, other: DenseElement) -> DenseElement:
        self._same(other)
        t, J, modulus = self.ring.t, self.ring.J, self.ring.modulus
        pairs = t * (t + 1) // 2 * (J * (J + 1) // 2)
        if 2 * pairs > self.ring.budgets.max_scalar_operations:
            raise LeanAlgebraRefusal('dense convolution operation budget exceeded')
        output = [0] * self.ring.slots
        for u in range(t):
            for z in range(J):
                for a in range(u + 1):
                    for b in range(z + 1):
                        output[u * J + z] += (
                            self.coefficients[a * J + b] * other.coefficients[(u - a) * J + z - b]
                        )
                        output[u * J + z] %= modulus
        return self.ring.element(output)

    def __pow__(self, exponent: int) -> DenseElement:
        exponent = natural(exponent, 'exponent')
        pairs = self.ring.t * (self.ring.t + 1) // 2 * (self.ring.J * (self.ring.J + 1) // 2)
        multiplications = exponent.bit_count() + max(0, exponent.bit_length() - 1)
        if 2 * pairs * multiplications > self.ring.budgets.max_scalar_operations:
            raise LeanAlgebraRefusal('dense power operation budget exceeded')
        result, base = self.ring.scalar(1), self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent >>= 1
            if exponent:
                base = base * base
        return result

    def readout(self) -> tuple[int, ...]:
        """Sum U coefficients only after products have been truncated."""
        if self.ring.slots > self.ring.budgets.max_scalar_operations:
            raise LeanAlgebraRefusal('dense readout operation budget exceeded')
        return tuple(
            sum(self.coefficients[u * self.ring.J + z] for u in range(self.ring.t))
            % self.ring.modulus
            for z in range(self.ring.J)
        )


class _ScalarPolynomials:
    def __init__(self, modulus: int, meter: AlgebraMeter):
        self.modulus, self.meter = modulus, meter

    def trim(self, values: Sequence[int]) -> tuple[int, ...]:
        result = [c % self.modulus for c in values]
        while result and not result[-1]:
            result.pop()
        return tuple(result)

    def add(self, left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
        length = max(len(left), len(right))
        self.meter.operations(length)
        return self.trim(
            tuple(
                (left[i] if i < len(left) else 0) + (right[i] if i < len(right) else 0)
                for i in range(length)
            )
        )

    def scale(self, scalar: int, values: Sequence[int]) -> tuple[int, ...]:
        self.meter.operations(len(values))
        return self.trim(tuple(scalar * c for c in values))

    def mul(self, left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
        if not left or not right:
            return ()
        self.meter.operations(2 * len(left) * len(right))
        result = [0] * (len(left) + len(right) - 1)
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                result[i + j] = (result[i + j] + a * b) % self.modulus
        return self.trim(result)

    def divmod(self, numerator: Sequence[int], denominator: Sequence[int]):
        if not denominator or denominator[-1] != 1:
            raise ValueError('digit division requires a monic divisor')
        degree = len(denominator) - 1
        steps = max(0, len(numerator) - degree)
        self.meter.operations(2 * len(denominator) * steps)
        remainder, quotient = list(numerator), [0] * steps
        for i in reversed(range(steps)):
            coefficient = remainder[i + degree] % self.modulus
            quotient[i] = coefficient
            for j, d in enumerate(denominator):
                remainder[i + j] = (remainder[i + j] - coefficient * d) % self.modulus
        return self.trim(quotient), self.trim(remainder[:degree])


@dataclass(frozen=True)
class TernaryIndicator:
    t: int
    modulus: int
    digits: tuple[tuple[int, int, int], ...]
    polynomial: tuple[int, ...]
    newton_steps: int
    scalar_operations: int
    workspace_coefficients: int

    @property
    def coefficient_word_bits(self) -> int:
        return 2 * self.t

    def encode(self) -> bytes:
        """Lean's 3*t base coefficients, each with baseWidth 3 t = 2*t bits."""
        return bytes(
            (coefficient >> bit) & 1
            for row in self.digits
            for coefficient in row
            for bit in range(self.coefficient_word_bits)
        )

    def row(self, value: int) -> tuple[int, ...]:
        value = operator.index(value) % self.modulus
        defect = value * (value - 1) * (value - 2) % self.modulus
        power, output = 1, []
        for a, b, c in self.digits:
            output.append((a + b * value + c * value * value) * power % self.modulus)
            power = power * defect % self.modulus
        return tuple(output)

    def evaluate(self, value: int) -> int:
        return sum(self.row(value)) % self.modulus


def construct_ternary_indicator(
    t: int,
    *,
    budgets: LeanAlgebraBudgets | None = None,
    meter: AlgebraMeter | None = None,
) -> TernaryIndicator:
    """Compute the code-derived coefficients; no supplied indicator or table oracle."""
    t = natural(t, 't')
    if not t:
        raise ValueError('indicator precision must be positive')
    limits = budgets or (meter.limits if meter is not None else LeanAlgebraBudgets())
    if meter is not None and meter.limits != limits:
        raise ValueError('the indicator and its meter must use the same budgets')
    if t > limits.max_precision:
        raise LeanAlgebraRefusal('ternary precision budget exceeded')
    if 3 * t > limits.max_dense_slots:
        raise LeanAlgebraRefusal('ternary digit slot budget exceeded')
    # Covers modulus, Newton operands/products/division scratch and saved digits.
    # Normalized degree is <3*t and every Newton product has degree <9*t.
    workspace = 60 * t + 32
    used = meter or AlgebraMeter(limits)
    used.storage(workspace, workspace * (2 * t))
    start = used.scalar_operations
    modulus = 3**t
    algebra = _ScalarPolynomials(modulus, used)
    H = (0, 2 % modulus, (-3) % modulus, 1)
    divisor = (1,)
    for _ in range(t):
        divisor = algebra.mul(divisor, H)
    half = (modulus + 1) // 2
    current = algebra.trim((1, -3 * half, half))
    current = algebra.divmod(current, divisor)[1]
    steps = (t - 1).bit_length()
    for _ in range(steps):
        square = algebra.mul(current, current)
        cube = algebra.mul(square, current)
        next_value = algebra.add(algebra.scale(3, square), algebra.scale(-2, cube))
        current = algebra.divmod(next_value, divisor)[1]
    polynomial = current + (0,) * (3 * t - len(current))
    digits = []
    for _ in range(t):
        current, remainder = algebra.divmod(current, H)
        digits.append(tuple(remainder) + (0,) * (3 - len(remainder)))
    if current:
        raise ArithmeticError('the normalized ternary polynomial has a nonzero high digit')
    return TernaryIndicator(
        t,
        modulus,
        tuple(digits),
        polynomial,
        steps,
        used.scalar_operations - start,
        workspace,
    )
