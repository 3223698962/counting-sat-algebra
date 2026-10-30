"""Exact dyadic search from SqrtLogSpace/DisjointReuse/SparseBudget.lean.

These numerical conditions select M. They are not a proof that the Python
leaf evaluator or the complete pipeline meets the formalized time bound.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass


def _positive(value: int, name: str) -> int:
    value = operator.index(value)
    if value < 1:
        raise ValueError(f'{name} must be positive')
    return value


def occurrence_cap(k: int, M: int) -> int:
    return 2 * sum(M ** j for j in range(k))


def depth_constant(k: int, M: int) -> int:
    return 4 * sum((3 * M) ** j for j in range(k - 1))


def right_constant(k: int) -> int:
    return 4 * sum(3 ** j for j in range(k - 1))


def parameter_test(k: int, numerator: int, denominator: int, exponent: int) -> bool:
    """Lean's parameterTest, with an integer ceiling log and no float oracle."""
    k = _positive(k, 'k')
    numerator = _positive(numerator, 'numerator')
    denominator = _positive(denominator, 'denominator')
    exponent = operator.index(exponent)
    if exponent < 0:
        raise ValueError('exponent must be nonnegative')
    M = 1 << exponent
    argument = 4 * depth_constant(k, M) * M
    ceiling_log = max(0, (argument - 1).bit_length())
    gamma = right_constant(k)
    return M >= 2 and gamma <= M and gamma * ceiling_log * denominator <= numerator * M


@dataclass(frozen=True)
class SparseParameter:
    k: int
    numerator: int
    denominator: int
    exponent: int
    M: int
    occurrence_bound: int
    depth_coefficient: int
    right_coefficient: int
    candidates_tested: int


def choose_sparse_parameter(
    k: int, numerator: int = 1, denominator: int = 1, *, max_exponent: int = 128,
) -> SparseParameter:
    """Choose the first passing power of two, as Lean's chooseExponent does."""
    k = _positive(k, 'k')
    numerator = _positive(numerator, 'numerator')
    denominator = _positive(denominator, 'denominator')
    limit = operator.index(max_exponent)
    if limit < 0:
        raise ValueError('max_exponent must be nonnegative')
    for exponent in range(limit + 1):
        if parameter_test(k, numerator, denominator, exponent):
            M = 1 << exponent
            return SparseParameter(
                k, numerator, denominator, exponent, M, occurrence_cap(k, M),
                depth_constant(k, M), right_constant(k), exponent + 1,
            )
    raise ValueError('dyadic parameter search exceeds the exponent budget')
