"""Exact numerical hypotheses of reference Lemma 4.7, including its tail."""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .replacement_matrix import default_interval_start
from .scheduler import BalanceCheck, check_balance


@dataclass(frozen=True)
class ReferenceParameters:
    balance: BalanceCheck
    tau: Fraction
    interval: tuple[int, ...]
    tail_probability: Fraction
    tail_limit: Fraction | None
    admitted: bool
    reasons: tuple[str, ...]


def check_reference_parameters(
    n: int, scopes: Sequence[Sequence[int]], *, r: int, Q: int, R: int,
    tau: Fraction, A0: int | None = None,
) -> ReferenceParameters:
    """Check Lemma 4.7 without allocating auxiliary coordinates or schedules.

    The exact tail denominator is (2r)**Q. Small parameters can pass this
    lemma even though the asymptotic choices in Lemma 4.8 are much larger.
    Admission is a check of numerical hypotheses, not an executed counter.
    """
    balance = check_balance(n, scopes, r, Q, R)
    r, Q, R = balance.r, balance.Q, balance.R
    assert R is not None
    tau = Fraction(tau)
    if A0 is None:
        A0 = default_interval_start(r, Q, R)
    else:
        import operator
        A0 = operator.index(A0)
    interval = tuple(range(A0, A0 + R))
    tail = Fraction(sum(
        math.comb(Q, j) * (r - 1) ** j * (r + 1) ** (Q - j)
        for j in range(Q + 1) if j not in interval
    ), (2 * r) ** Q)
    limit = tau / (r * balance.b) if balance.b else None
    reasons = list(balance.reasons)
    if not 0 < tau <= Fraction(1, 2):
        reasons.append("need 0 < tau <= 1/2")
    if not 0 <= A0 or A0 + R - 1 > Q:
        reasons.append("the interval must lie inside {0, ..., Q}")
    if limit is None or tail > limit:
        reasons.append(f"binomial tail {tail} exceeds tau/(rb) = {limit}")
    return ReferenceParameters(balance, tau, interval, tail, limit, not reasons, tuple(reasons))
