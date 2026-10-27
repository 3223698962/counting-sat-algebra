"""Coordinate access to Section 4's full auxiliary operator without a dense matrix.

An orbit of the commuting swaps is a Boolean cube of dimension k <= Q.
Its Fourier eigenvalues are f_I(j). Krawtchouk coefficients evaluate one
matrix entry in polynomial time in Q and R, including at R=32, Q=192.
This representation provides entries; it does not enumerate all coordinates,
measurement modes, or the support of an enormous row.
"""

from __future__ import annotations

import itertools
import math
import operator
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .replacement_matrix import (
    _mismatch_positions, _phi, binomial_exceptional_count,
    default_interval_start, scalar_f_coefficients,
)
from .replacements import eval_poly

Label = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass(frozen=True)
class StructuredAuxiliary:
    r: int
    Q: int
    R: int
    A0: int
    f_coefficients: tuple[Fraction, ...]
    f_values: tuple[Fraction, ...]

    @property
    def dimension(self) -> int:
        return self.r ** (2 * self.Q)

    @property
    def interval(self) -> tuple[int, ...]:
        return tuple(range(self.A0, self.A0 + self.R))

    @property
    def measurement_mode_count(self) -> int:
        return binomial_exceptional_count(self.r, self.Q, self.interval)

    def _label(self, label: Sequence[Sequence[int]]) -> Label:
        if len(label) != 2 or any(len(sequence) != self.Q for sequence in label):
            raise ValueError("a label must contain two length-Q sequences")
        parsed = tuple(tuple(operator.index(value) for value in sequence) for sequence in label)
        if any(value < 0 or value >= self.r for sequence in parsed for value in sequence):
            raise ValueError("label symbols must lie in {0, ..., r-1}")
        return parsed  # type: ignore[return-value]

    def label_at(self, index: int) -> Label:
        index = operator.index(index)
        if not 0 <= index < self.dimension:
            raise ValueError("coordinate index outside the auxiliary dimension")
        symbols = [0] * (2 * self.Q)
        for position in range(len(symbols) - 1, -1, -1):
            index, symbols[position] = divmod(index, self.r)
        return tuple(symbols[:self.Q]), tuple(symbols[self.Q:])

    def matrix_entry(self, row: Sequence[Sequence[int]], column: Sequence[Sequence[int]]) -> Fraction:
        """Exact M0[u,y] = (V0 f_I(W0))[u,y] without enumerating the orbit."""
        row, column = self._label(row), self._label(column)
        active = 0
        mismatches = 0
        for a, d, c, e in zip(*row, *column):
            if a == d:
                if (c, e) != (a, d):
                    return Fraction(0)
            else:
                active += 1
                if (c, e) == (a, d):
                    mismatches += 1
                elif (c, e) != (d, a):
                    return Fraction(0)
        # A polynomial in the swaps of degree R-1 cannot reach this distance.
        if mismatches >= self.R:
            return Fraction(0)
        value = Fraction(0)
        for j in range(active + 1):
            weight = sum(
                (-1) ** i * math.comb(mismatches, i) * math.comb(active - mismatches, j - i)
                for i in range(max(0, j - active + mismatches), min(j, mismatches) + 1)
            )
            value += weight * self.f_values[j]
        return value / (2 ** active)

    def replacement_entry(self, row: Sequence[Sequence[int]], column: Sequence[Sequence[int]]) -> tuple[Fraction, ...]:
        """The coefficient polynomial P_u(T)[y] from equation (4.10)."""
        row, column = self._label(row), self._label(column)
        matrix = self.matrix_entry(row, column)
        diagonal = Fraction(row == column)
        if not matrix:
            return (diagonal,)
        poly = [matrix * value for value in _phi(_mismatch_positions(row, column))]
        poly[0] += diagonal - matrix
        while len(poly) > 1 and not poly[-1]:
            poly.pop()
        return tuple(poly)

    def replacement_support_bound(self, row: Sequence[Sequence[int]]) -> int:
        """Upper bound on nonzero coordinates of P_u, without walking its orbit."""
        row = self._label(row)
        active = sum(a != d for a, d in zip(*row))
        return 1 + sum(math.comb(active, j) for j in range(1, min(active, self.R - 1) + 1))

    def replacement_support(self, row: Sequence[Sequence[int]], *, max_terms: int = 4096):
        """Enumerate only the bounded swap support of equation (4.10)."""
        row = self._label(row)
        max_terms = operator.index(max_terms)
        if max_terms < 0:
            raise ValueError("max_terms must be nonnegative")
        if self.replacement_support_bound(row) > max_terms:
            raise ValueError("replacement support exceeds max_terms; no support enumerated")
        diagonal = self.replacement_entry(row, row)
        if any(diagonal):
            yield row, diagonal
        active = tuple(s for s, (a, d) in enumerate(zip(*row)) if a != d)
        for size in range(1, min(len(active), self.R - 1) + 1):
            for positions in itertools.combinations(active, size):
                first, second = list(row[1]), list(row[0])
                for s in positions:
                    first[s], second[s] = second[s], first[s]
                column = tuple(first), tuple(second)
                if column == row:
                    continue
                poly = self.replacement_entry(row, column)
                if any(poly):
                    yield column, poly

    def measurement_entry(self, mode: Sequence[tuple[int, int, int]], label: Sequence[Sequence[int]]) -> int:
        """Entry of one tensor eigenvector; kind 0=diagonal, +1=plus, -1=minus."""
        label = self._label(label)
        if len(mode) != self.Q:
            raise ValueError("measurement mode must have length Q")
        result = 1
        negatives = 0
        for (kind, a, d), c, e in zip(mode, *label):
            if kind == 0 and a == d and 0 <= a < self.r:
                result *= int((c, e) == (a, a))
            elif kind in (-1, 1) and 0 <= a < d < self.r:
                negatives += int(kind == -1)
                result *= 1 if (c, e) == (a, d) else (kind if (c, e) == (d, a) else 0)
            else:
                raise ValueError("invalid local eigenvector descriptor")
        if negatives in self.interval:
            raise ValueError("this eigenvector is inside I and is not a measurement mode")
        return result

    def measurement_modes(self, *, max_modes: int = 4096):
        """Enumerate descriptors only when their exact number is under the cap."""
        max_modes = operator.index(max_modes)
        if max_modes < 0:
            raise ValueError("max_modes must be nonnegative")
        if self.measurement_mode_count > max_modes:
            raise ValueError(f"{self.measurement_mode_count} measurement modes exceed cap {max_modes}; no modes enumerated")
        if self.measurement_mode_count == 0:
            return
        plus = [(0, a, a) for a in range(self.r)]
        plus += [(1, a, d) for a in range(self.r) for d in range(a + 1, self.r)]
        minus = [(-1, a, d) for a in range(self.r) for d in range(a + 1, self.r)]
        for count in range(self.Q + 1):
            if count in self.interval:
                continue
            for positions in itertools.combinations(range(self.Q), count):
                negative = set(positions)
                yield from itertools.product(*(minus if s in negative else plus for s in range(self.Q)))


def build_structured_auxiliary(r: int, Q: int, R: int, A0: int | None = None) -> StructuredAuxiliary:
    """Build a polynomial-size descriptor; the dense slice cap does not apply."""
    r, Q, R = operator.index(r), operator.index(Q), operator.index(R)
    if r < 2 or Q < 1 or not 1 <= R <= Q + 1 or Q > 4096:
        raise ValueError("need r >= 2, 1 <= Q <= 4096, 1 <= R <= Q+1")
    A0 = default_interval_start(r, Q, R) if A0 is None else operator.index(A0)
    if A0 < 0 or A0 + R - 1 > Q:
        raise ValueError("the interpolation interval must lie inside {0, ..., Q}")
    coeffs = scalar_f_coefficients(A0, R)
    values = tuple(eval_poly(list(coeffs), Fraction(j)) for j in range(Q + 1))
    return StructuredAuxiliary(r, Q, R, A0, coeffs, values)
