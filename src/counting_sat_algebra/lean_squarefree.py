"""Budgeted support-union algebra from SqrtLog/SquarefreeBasic.lean.

X_v**2 = X_v. Boolean summation uses the monomial identity
sum_x X_S = 2**(n-|S|), without an assignment table. This finite sparse
representation can still have exponentially many terms; it does not supply
the streamed packing and degree bounds of the complete Lean leaf program.
"""

from __future__ import annotations

import operator
import weakref
from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping, Sequence

from .lean_ring import AlgebraMeter, LeanAlgebraRefusal, natural

Support = frozenset[int]


@dataclass(eq=False, frozen=True)
class SquarefreePolynomial:
    algebra: SquarefreeAlgebra
    terms: Mapping[Support, int]

    def __add__(self, other: SquarefreePolynomial) -> SquarefreePolynomial:
        return self.algebra.add(self, other)

    def __sub__(self, other: SquarefreePolynomial) -> SquarefreePolynomial:
        return self.algebra.add(self, self.algebra.scale(-1, other))

    def __mul__(self, other: SquarefreePolynomial) -> SquarefreePolynomial:
        return self.algebra.multiply(self, other)

    def scale(self, coefficient: int) -> SquarefreePolynomial:
        return self.algebra.scale(coefficient, self)

    def evaluate(self, assignment: Sequence[int]) -> int:
        if len(assignment) != self.algebra.n or any(bit not in (0, 1) for bit in assignment):
            raise ValueError('a Boolean assignment must have exactly n bits')
        return (
            sum(
                coefficient
                for support, coefficient in self.terms.items()
                if all(assignment[v] for v in support)
            )
            % self.algebra.modulus
        )


class SquarefreeAlgebra:
    """Account for live immutable tables, scratch tables and prepared data.

    A table entry is charged an n-bit support mask, one fixed-width residue
    and eight framing bits. The receipt records logical encoded storage;
    Python objects and visited machine tapes are not being measured. Weak
    references only release charges after a polynomial is no longer retained.
    """

    def __init__(
        self,
        n: int,
        modulus: int,
        meter: AlgebraMeter,
        *,
        coefficient_bits: int | None = None,
        parked_coefficients: int = 0,
        parked_bits: int = 0,
    ):
        self.n, self.modulus = natural(n, 'n'), natural(modulus, 'modulus')
        if self.modulus < 2:
            raise ValueError('modulus must be at least two')
        self.meter = meter
        if self.n > meter.limits.max_variables:
            raise LeanAlgebraRefusal('squarefree variable budget exceeded')
        self.coefficient_bits = (
            (modulus - 1).bit_length()
            if coefficient_bits is None
            else natural(
                coefficient_bits,
                'coefficient_bits',
            )
        )
        if (modulus - 1).bit_length() > self.coefficient_bits:
            raise ValueError('the coefficient word cannot encode every residue')
        self.parked_coefficients, self.parked_bits = parked_coefficients, parked_bits
        self._live: weakref.WeakSet[SquarefreePolynomial] = weakref.WeakSet()
        self._reserve(0)
        self.zero, self.one = self.constant(0), self.constant(1)

    def _reserve(self, scratch: int) -> None:
        if scratch > self.meter.limits.max_polynomial_terms:
            raise LeanAlgebraRefusal('squarefree polynomial term budget exceeded')
        live = sum(len(poly.terms) for poly in self._live)
        coefficients = self.parked_coefficients + live + scratch
        bits = self.parked_bits + (live + scratch) * (self.n + self.coefficient_bits + 8)
        self.meter.storage(coefficients, bits)
        self.meter.peak_polynomial_terms = max(self.meter.peak_polynomial_terms, scratch)

    def _make(self, terms: dict[Support, int]) -> SquarefreePolynomial:
        self._reserve(len(terms))
        result = SquarefreePolynomial(self, MappingProxyType(terms))
        self._live.add(result)
        return result

    def reserve_scalars(self, count: int) -> None:
        """Charge the output array and modular multiplication/readout scratch."""
        live = sum(len(poly.terms) for poly in self._live)
        self.meter.storage(
            self.parked_coefficients + live + count,
            self.parked_bits
            + live * (self.n + self.coefficient_bits + 8)
            + count * self.coefficient_bits,
        )

    def _same(self, *polynomials: SquarefreePolynomial) -> None:
        if any(poly.algebra is not self for poly in polynomials):
            raise ValueError('polynomial operands must belong to the same squarefree algebra')

    def polynomial(self, terms: Mapping[Support, int]) -> SquarefreePolynomial:
        if len(terms) > self.meter.limits.max_polynomial_terms:
            raise LeanAlgebraRefusal('squarefree input term budget exceeded')
        output = {}
        for support, coefficient in terms.items():
            support = frozenset(natural(v, 'variable') for v in support)
            if any(v >= self.n for v in support):
                raise ValueError('squarefree support exceeds the declared n')
            coefficient = operator.index(coefficient) % self.modulus
            if coefficient:
                self._reserve(len(output) + 1)
                output[support] = coefficient
        return self._make(output)

    def constant(self, coefficient: int) -> SquarefreePolynomial:
        coefficient = operator.index(coefficient) % self.modulus
        return self._make({frozenset(): coefficient} if coefficient else {})

    def variable(self, variable: int) -> SquarefreePolynomial:
        variable = natural(variable, 'variable')
        if variable >= self.n:
            raise ValueError('variable exceeds the declared n')
        return self._make({frozenset((variable,)): 1})

    def _work_fits(self, products: int, operations: int) -> None:
        if self.meter.monomial_products + products > self.meter.limits.max_monomial_products:
            raise LeanAlgebraRefusal('squarefree monomial product budget exceeded')
        if self.meter.scalar_operations + operations > self.meter.limits.max_scalar_operations:
            raise LeanAlgebraRefusal('scalar operation budget exceeded')

    def _put(self, output: dict[Support, int], support: Support, coefficient: int) -> None:
        if coefficient:
            if support not in output:
                self._reserve(len(output) + 1)
            output[support] = coefficient
        else:
            output.pop(support, None)

    def add(self, left: SquarefreePolynomial, right: SquarefreePolynomial):
        self._same(left, right)
        self._work_fits(0, len(right.terms))
        self._reserve(len(left.terms))
        output = dict(left.terms)
        for support, coefficient in right.terms.items():
            self.meter.operations(1)
            self._put(output, support, (output.get(support, 0) + coefficient) % self.modulus)
        return self._make(output)

    def scale(self, coefficient: int, value: SquarefreePolynomial):
        self._same(value)
        coefficient = operator.index(coefficient) % self.modulus
        if not coefficient:
            return self.zero
        self._work_fits(0, len(value.terms))
        output = {}
        for support, c in value.terms.items():
            self.meter.operations(1)
            self._put(output, support, coefficient * c % self.modulus)
        return self._make(output)

    def multiply(self, left: SquarefreePolynomial, right: SquarefreePolynomial):
        self._same(left, right)
        products = len(left.terms) * len(right.terms)
        self._work_fits(products, 2 * products)
        output: dict[Support, int] = {}
        for a, c in left.terms.items():
            for b, d in right.terms.items():
                support = a | b
                self.meter.operations(2)
                self.meter.monomial_products += 1
                self._put(output, support, (output.get(support, 0) + c * d) % self.modulus)
        return self._make(output)

    def marker_product(
        self,
        left: Sequence[SquarefreePolynomial],
        right: Sequence[SquarefreePolynomial],
    ) -> tuple[SquarefreePolynomial, ...]:
        if len(left) != len(right):
            raise ValueError('marker arrays must have the same cutoff')
        self._same(*left, *right)
        output = [self.zero] * len(left)
        for u, a in enumerate(left):
            if not a.terms:
                continue
            for v in range(len(right) - u):
                if right[v].terms:
                    output[u + v] = output[u + v] + a * right[v]
        return tuple(output)

    def sum_boolean(
        self, value: SquarefreePolynomial, *, fixed: frozenset[int] = frozenset()
    ) -> int:
        """Sum all free bits through monomial weights; no Boolean cube is visited."""
        self._same(value)
        if any(v < 0 or v >= self.n for v in fixed):
            raise ValueError('fixed variable exceeds the declared n')
        free = self.n - len(fixed)
        total = 0
        for support, coefficient in value.terms.items():
            if support & fixed:
                raise ValueError('a fixed coordinate was not substituted out')
            weight = self.power_two(free - len(support))
            self.meter.operations(2)
            total = (total + coefficient * weight) % self.modulus
        return total

    def power_two(self, exponent: int) -> int:
        exponent = natural(exponent, 'exponent')
        result, base = 1, 2 % self.modulus
        while exponent:
            if exponent & 1:
                self.meter.operations(1)
                result = result * base % self.modulus
            exponent >>= 1
            if exponent:
                self.meter.operations(1)
                base = base * base % self.modulus
        return result
