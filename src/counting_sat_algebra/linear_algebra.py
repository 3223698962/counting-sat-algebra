"""Reusable exact row spans, with sparse elimination over the rationals."""

from __future__ import annotations

from fractions import Fraction
from typing import Mapping, Sequence
import operator


class RationalRowSpan:
    """Eliminate a matrix once and test many targets against its row span.

    A target in the row span annihilates the matrix's entire rational kernel.
    Consequently equal measurements imply equal target values on Boolean
    assignments as well; checking those assignments adds no condition.
    """

    def __init__(self, rows: Sequence[Sequence[Fraction]], width: int):
        if width < 0 or any(len(row) != width for row in rows):
            raise ValueError("row widths must match the nonnegative width")
        self.width = width
        basis: dict[int, dict[int, Fraction]] = {}
        independent = []
        for row in rows:
            residual = {i: Fraction(value) for i, value in enumerate(row) if value}
            for pivot in sorted(basis):
                self._subtract(residual, pivot, basis[pivot])
            if residual:
                pivot = min(residual)
                scale = residual[pivot]
                basis[pivot] = {i: value / scale for i, value in residual.items()}
                independent.append(tuple(Fraction(value) for value in row))
        self._basis = tuple((pivot, basis[pivot]) for pivot in sorted(basis))
        self.independent_rows = tuple(independent)

    @staticmethod
    def _subtract(residual: dict[int, Fraction], pivot: int, row: dict[int, Fraction]) -> None:
        scale = residual.get(pivot, Fraction(0))
        if not scale:
            return
        for i, value in row.items():
            updated = residual.get(i, Fraction(0)) - scale * value
            if updated:
                residual[i] = updated
            else:
                residual.pop(i, None)

    @property
    def rank(self) -> int:
        return len(self._basis)

    def contains(self, target: Sequence[Fraction]) -> bool:
        if len(target) != self.width:
            raise ValueError("target width must match the row span")
        residual = {i: Fraction(value) for i, value in enumerate(target) if value}
        for pivot, row in self._basis:
            self._subtract(residual, pivot, row)
            if not residual:
                return True
        return not residual


class SparseEliminationBudgetExceeded(ValueError):
    """Sparse elimination would exceed its stored-coefficient budget."""


class SparseRationalRowSpan:
    """Exact echelon basis whose rows and targets contain only nonzero entries.

    Pivot lookup follows the leading residual coordinate. No array of length
    ``width`` is allocated, and testing a target does not scan every pivot.
    The independent rows remain original integer rows when the input is integer.
    """

    def __init__(
        self, rows: Sequence[Mapping[int, Fraction]], width: int, *,
        max_entries: int = 2_000_000,
    ):
        self.width = operator.index(width)
        self.max_entries = operator.index(max_entries)
        if min(self.width, self.max_entries) < 0:
            raise ValueError("width and max_entries must be nonnegative integers")
        self._basis: dict[int, dict[int, Fraction]] = {}
        self.stored_entries = 0
        independent = []
        for row in rows:
            original = self._normalize(row)
            residual = dict(original)
            self._reduce(residual)
            if residual:
                pivot = min(residual)
                scale = residual[pivot]
                normalized = {i: value / scale for i, value in residual.items()}
                self.stored_entries += len(normalized)
                if self.stored_entries > self.max_entries:
                    raise SparseEliminationBudgetExceeded("sparse elimination exceeded its coefficient budget")
                self._basis[pivot] = normalized
                independent.append(original)
        self.independent_rows = tuple(independent)

    def _normalize(self, row: Mapping[int, Fraction]) -> dict[int, Fraction]:
        normalized = {}
        for raw_index, raw_value in row.items():
            index = operator.index(raw_index)
            if not 0 <= index < self.width:
                raise ValueError("sparse coordinate outside the declared width")
            value = Fraction(raw_value)
            if value:
                normalized[index] = value
        return normalized

    def _reduce(self, residual: dict[int, Fraction]) -> None:
        while residual:
            if self.stored_entries + len(residual) > self.max_entries:
                raise SparseEliminationBudgetExceeded("sparse elimination exceeded its coefficient budget")
            pivot = min(residual)
            row = self._basis.get(pivot)
            if row is None:
                return
            RationalRowSpan._subtract(residual, pivot, row)

    @property
    def rank(self) -> int:
        return len(self._basis)

    def contains(self, target: Mapping[int, Fraction]) -> bool:
        residual = self._normalize(target)
        self._reduce(residual)
        return not residual
