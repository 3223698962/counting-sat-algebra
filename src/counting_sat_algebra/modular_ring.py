"""Exact finite coefficient rings from Appendix D.4.1.

An element of A_p = ((Z/2^p)[U]/f)[V]/(V^ell) is a tuple of integer
residues, in V-major order: index j*e+i represents U^i V^j.  The binary
polynomial f has an irreducible reduction modulo two.  With e=1 this is
the base ring B_{p,ell}.  No floating point or division by a nonunit is
used.  The finite primitives here do not prove Theorem D.1.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property, lru_cache
from typing import Sequence


Element = tuple[int, ...]


class ModularBudgetRefusal(ValueError):
    """A finite Appendix D operation exceeds an explicit resource budget."""


def _positive(value: int, name: str) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _gf2_remainder(value: int, divisor: int) -> int:
    while value.bit_length() >= divisor.bit_length():
        value ^= divisor << (value.bit_length() - divisor.bit_length())
    return value


def _gf2_gcd(left: int, right: int) -> int:
    while right:
        left, right = right, _gf2_remainder(left, right)
    return left


def _gf2_mul(left: int, right: int, polynomial: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left.bit_length() == polynomial.bit_length():
            left ^= polynomial
    return result


def _gf2_pow(value: int, exponent: int, polynomial: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = _gf2_mul(result, value, polynomial)
        exponent >>= 1
        value = _gf2_mul(value, value, polynomial)
    return result


@lru_cache(maxsize=128)
def is_irreducible_mod2(polynomial: int) -> bool:
    """Test a monic binary polynomial by exact Frobenius/gcd computations."""
    if not isinstance(polynomial, int) or polynomial < 2:
        return False
    degree = polynomial.bit_length() - 1
    if degree == 1:
        return True
    if not polynomial & 1:
        return False
    power = 2  # U
    for _ in range(1, degree // 2 + 1):
        power = _gf2_mul(power, power, polynomial)
        if _gf2_gcd(power ^ 2, polynomial) != 1:
            return False
    # A reducible degree-e polynomial has a factor of degree at most e/2.
    return True


def find_irreducible_polynomial(
    degree: int, *, max_degree: int = 32, max_candidates: int = 4096,
) -> int:
    """Deterministic search, without allocating the candidate list."""
    _positive(degree, "degree")
    _positive(max_degree, "max_degree")
    _positive(max_candidates, "max_candidates")
    if degree > max_degree:
        raise ModularBudgetRefusal("extension degree exceeds budget; search not started")
    for trial, lower in enumerate(range(1, 1 << degree, 2)):
        if trial >= max_candidates:
            raise ModularBudgetRefusal("irreducible-polynomial candidate budget exhausted")
        candidate = (1 << degree) | lower
        if is_irreducible_mod2(candidate):
            return candidate
    raise ArithmeticError("no irreducible polynomial of the requested degree")


def extension_ring_for_size(
    s: int, precision: int, ell: int = 1, *, max_degree: int = 32,
    max_candidates: int = 4096, max_precision_bits: int = 4096, max_coordinates: int = 4096,
) -> ModularRing:
    """D.4.1's e=4*ceil(log2(s)), with guards before search or ring allocation."""
    _positive(s, "s")
    _positive(precision, "precision")
    _positive(ell, "ell")
    if s < 2:
        raise ValueError("the size parameter s must be at least two")
    degree = 4 * (s - 1).bit_length()
    if precision > max_precision_bits or degree * ell > max_coordinates:
        raise ModularBudgetRefusal("coefficient ring exceeds budget; search not started")
    polynomial = find_irreducible_polynomial(degree, max_degree=max_degree, max_candidates=max_candidates)
    return ModularRing(precision, ell, polynomial)


@dataclass(frozen=True)
class ModularRing:
    precision: int
    ell: int = 1
    polynomial: int = 0b11

    def __post_init__(self) -> None:
        _positive(self.precision, "precision")
        _positive(self.ell, "ell")
        if not is_irreducible_mod2(self.polynomial):
            raise ValueError("f must be monic and irreducible modulo two")

    @property
    def extension_degree(self) -> int:
        return self.polynomial.bit_length() - 1

    @property
    def size(self) -> int:
        return self.extension_degree * self.ell

    @cached_property
    def modulus(self) -> int:
        return 1 << self.precision

    @cached_property
    def zero(self) -> Element:
        return (0,) * self.size

    @cached_property
    def one(self) -> Element:
        return (1,) + (0,) * (self.size - 1)

    def at_precision(self, precision: int) -> ModularRing:
        return ModularRing(precision, self.ell, self.polynomial)

    def element(self, coordinates: Sequence[int]) -> Element:
        """Normalize and zero-pad V-major coordinates; never truncate input."""
        if len(coordinates) > self.size:
            raise ValueError("too many ring coordinates")
        if any(not isinstance(x, int) for x in coordinates):
            raise ValueError("ring coordinates must be integers")
        return tuple(x % self.modulus for x in coordinates) + (0,) * (self.size - len(coordinates))

    def scalar(self, value: int) -> Element:
        return self.element((value,))

    def from_base(self, coefficients: Sequence[int]) -> Element:
        """Embed coefficients of 1,V,... from B_{p,ell}."""
        if len(coefficients) > self.ell:
            raise ValueError("too many V coefficients")
        coordinates = [0] * self.size
        for j, value in enumerate(coefficients):
            coordinates[j * self.extension_degree] = value
        return self.element(coordinates)

    def base_coefficients(self, value: Element) -> tuple[int, ...]:
        self._check(value)
        e = self.extension_degree
        if any(value[j] % self.modulus for j in range(self.size) if j % e):
            raise ValueError("element is outside the base coefficient ring")
        return tuple(value[j * e] % self.modulus for j in range(self.ell))

    def from_fraction(self, value: Fraction | int) -> Element:
        value = Fraction(value)
        if value.denominator % 2 == 0:
            raise ValueError("a two-adically integral scalar must have odd denominator")
        return self.scalar(value.numerator * pow(value.denominator, -1, self.modulus))

    def _check(self, value: Element) -> None:
        if len(value) != self.size:
            raise ValueError("ring element has the wrong number of coordinates")

    def add(self, left: Element, right: Element) -> Element:
        self._check(left)
        self._check(right)
        return tuple((x + y) % self.modulus for x, y in zip(left, right))

    def sub(self, left: Element, right: Element) -> Element:
        self._check(left)
        self._check(right)
        return tuple((x - y) % self.modulus for x, y in zip(left, right))

    def scale(self, value: Element, scalar: int) -> Element:
        self._check(value)
        return tuple((scalar * x) % self.modulus for x in value)

    def _u_mul(self, left: Sequence[int], right: Sequence[int]) -> list[int]:
        e = self.extension_degree
        coefficients = [0] * (2 * e - 1)
        for i, x in enumerate(left):
            if x:
                for j, y in enumerate(right):
                    coefficients[i + j] += x * y
        for degree in range(2 * e - 2, e - 1, -1):
            value = coefficients[degree] % self.modulus
            if value:
                for i in range(e):
                    if self.polynomial >> i & 1:
                        coefficients[degree - e + i] -= value
        return [x % self.modulus for x in coefficients[:e]]

    def mul(self, left: Element, right: Element) -> Element:
        self._check(left)
        self._check(right)
        e = self.extension_degree
        result = [0] * self.size
        for i in range(self.ell):
            row = left[i * e:(i + 1) * e]
            if not any(row):
                continue
            for j in range(self.ell - i):
                column = right[j * e:(j + 1) * e]
                if not any(column):
                    continue
                for k, value in enumerate(self._u_mul(row, column)):
                    result[(i + j) * e + k] += value
        return tuple(x % self.modulus for x in result)

    def pow(self, value: Element, exponent: int) -> Element:
        self._check(value)
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("exponent must be a nonnegative integer")
        result = self.one
        while exponent:
            if exponent & 1:
                result = self.mul(result, value)
            exponent >>= 1
            if exponent:
                value = self.mul(value, value)
        return result

    def inverse(self, value: Element) -> Element:
        """Invert units by residue-field inversion, Newton lifting, then V-series."""
        self._check(value)
        e = self.extension_degree
        constant = value[:e]
        residue = sum((x & 1) << i for i, x in enumerate(constant))
        if not residue:
            raise ValueError("cannot invert a nonunit")
        binary = _gf2_pow(residue, (1 << e) - 2, self.polynomial)
        inverse = [(binary >> i) & 1 for i in range(e)]
        correct_bits = 1
        while correct_bits < self.precision:
            product = self._u_mul(constant, inverse)
            correction = [-x for x in product]
            correction[0] += 2
            inverse = self._u_mul(inverse, correction)
            correct_bits *= 2
        rows = [inverse]
        for j in range(1, self.ell):
            total = [0] * e
            for i in range(1, j + 1):
                product = self._u_mul(value[i * e:(i + 1) * e], rows[j - i])
                total = [x + y for x, y in zip(total, product)]
            rows.append(self._u_mul(inverse, [-x for x in total]))
        return tuple(x for row in rows for x in row)

    def canonical_lift(self, value: Element, target: ModularRing) -> Element:
        """Lift each coordinate separately; this is not a ring homomorphism."""
        self._check(value)
        if target.ell != self.ell or target.polynomial != self.polynomial:
            raise ValueError("canonical lifts require the same coefficient basis")
        if target.precision < self.precision:
            raise ValueError("canonical lifting cannot lower precision")
        return tuple(x % self.modulus for x in value)

    def residue_point(self, index: int) -> Element:
        """Binary-coordinate lift of a residue-field element."""
        e = self.extension_degree
        if not isinstance(index, int) or not 0 <= index < 1 << e:
            raise ValueError("residue point index outside the field")
        return self.element(tuple((index >> i) & 1 for i in range(e)))


def evaluate_polynomial(ring: ModularRing, coefficients: Sequence[Element], point: Element) -> Element:
    value = ring.zero
    for coefficient in reversed(coefficients):
        value = ring.add(ring.mul(value, point), coefficient)
    return value


def interpolate_polynomial(
    ring: ModularRing, points: Sequence[Element], values: Sequence[Element],
) -> tuple[Element, ...]:
    """Newton interpolation over A_p; every point difference must be a unit."""
    if not points or len(points) != len(values):
        raise ValueError("need equally many nonempty points and values")
    for value in (*points, *values):
        ring._check(value)
    divided = list(values)
    for width in range(1, len(points)):
        for j in range(len(points) - 1, width - 1, -1):
            divided[j] = ring.mul(
                ring.sub(divided[j], divided[j - 1]),
                ring.inverse(ring.sub(points[j], points[j - width])),
            )
    coefficients = [divided[-1]]
    for j in range(len(points) - 2, -1, -1):
        updated = [ring.zero] * (len(coefficients) + 1)
        for k, value in enumerate(coefficients):
            updated[k] = ring.sub(updated[k], ring.mul(points[j], value))
            updated[k + 1] = ring.add(updated[k + 1], value)
        updated[0] = ring.add(updated[0], divided[j])
        coefficients = updated
    return tuple(coefficients)
