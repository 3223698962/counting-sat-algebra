"""Exact numerical hypotheses of reference Lemma 4.7, and Lemma 4.8's schedule.

Lemma 4.8 fixes one integer schedule in n. For every n whose sparse-placement
budgets can run, that closed form still fails its preliminary guard, so no
auxiliary dimension is formed. A separate search then keeps a budget-feasible
tuple only when ``check_reference_parameters`` accepts it. Selecting such a
finite tuple is not a proof of the asymptotic time or space bound.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass, fields
from fractions import Fraction
from typing import Sequence

from .replacement_matrix import default_interval_start
from .scheduler import BALANCE_MAX_Q, BalanceCheck, check_balance


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
        A0 = operator.index(A0)
    interval = tuple(range(A0, A0 + R))
    tail = binomial_tail_probability(r, Q, interval)
    limit = tau / (r * balance.b) if balance.b else None
    reasons = list(balance.reasons)
    if not 0 < tau <= Fraction(1, 2):
        reasons.append("need 0 < tau <= 1/2")
    if not 0 <= A0 or A0 + R - 1 > Q:
        reasons.append("the interval must lie inside {0, ..., Q}")
    if limit is None or tail > limit:
        reasons.append(f"binomial tail {tail} exceeds tau/(rb) = {limit}")
    return ReferenceParameters(balance, tau, interval, tail, limit, not reasons, tuple(reasons))


def binomial_tail_probability(r: int, Q: int, interval: Sequence[int]) -> Fraction:
    """Exact Pr[W not in interval] for W ~ Bin(Q, (r-1)/(2r)).

    The denominator is (2r)**Q, as in Lemma 4.7. Interval points outside
    {0, ..., Q} do not change the sum. This is not Hoeffding's estimate.
    """
    r = operator.index(r)
    Q = operator.index(Q)
    if r < 1 or Q < 0:
        raise ValueError("the binomial tail needs r >= 1 and Q >= 0")
    inside = {operator.index(point) for point in interval}
    numerator = sum(
        math.comb(Q, j) * (r - 1) ** j * (r + 1) ** (Q - j)
        for j in range(Q + 1) if j not in inside
    )
    return Fraction(numerator, (2 * r) ** Q)


def _ceil_log2(value: int) -> int:
    """ceil(log2(value)) for value >= 1, with ceil(log2(1)) = 0."""
    value = operator.index(value)
    if value <= 1:
        return 0
    return (value - 1).bit_length()


def _ceil_sqrt(value: int) -> int:
    value = operator.index(value)
    if value < 0:
        raise ValueError("square root is not defined for a negative integer")
    root = math.isqrt(value)
    return root if root * root == value else root + 1


def lemma_4_8_integers(lam: int, beta: int = 1) -> tuple[int, int, int, int]:
    """Return (nu, r, Q, R) from lambda and beta, without a formula or a matrix.

    This is the integer schedule in the proof of Lemma 4.8. It does not form
    D = r**(2Q). Callers still have to apply the preliminary guard r >= 2 and
    Q > R before building tau or D.
    """
    lam = operator.index(lam)
    beta = operator.index(beta)
    if lam < 0 or beta < 1:
        raise ValueError("lambda must be nonnegative and beta must be >= 1")
    nu = _ceil_log2(lam) if lam >= 1 else 0
    return nu, _largest_r(lam, nu), (lam // (128 * nu) if nu else 0), (beta + 3) * _ceil_sqrt(lam)


def _largest_r(lam: int, nu: int) -> int:
    """Largest integer r with r**2 * nu**4 <= lam. Zero when nu is zero."""
    if nu <= 0 or lam <= 0:
        return 0
    scale = nu ** 4
    r = math.isqrt(lam // scale)
    while (r + 1) * (r + 1) * scale <= lam:
        r += 1
    while r > 0 and r * r * scale > lam:
        r -= 1
    return r


def _smallest_beta(lam: int, b: int, w: int) -> int:
    """Smallest integer beta >= 1 with b, w <= lam**beta, when that power grows."""
    need = max(b, w, 1)
    if lam <= 1:
        return 1
    beta = 1
    limit = lam
    while limit < need:
        if beta >= 64:
            break
        beta += 1
        limit *= lam
    return beta


def _dimension_within_power(dimension: int, lam: int) -> bool:
    """Return whether dimension**128 <= 2**lam, without building a matrix."""
    if dimension <= 1:
        return True
    if (dimension.bit_length() - 1) * 128 > lam:
        return False
    if dimension.bit_length() * 128 <= lam:
        return True
    return pow(dimension, 128) <= (1 << lam)


@dataclass(frozen=True)
class ChoiceBudgets:
    """Bounds for the finite search. They mirror the sparse-placement gates.

    ``max_dimension`` limits the integer D = r**(2Q) the search will form.
    It is not the dense-slice cap of 256, and forming D does not allocate a
    matrix. Alphabet and trial caps are the same numbers ``schedule_sparse_placement``
    uses before it starts the conditional-expectation search.
    """

    max_r: int = 16
    max_Q: int = 8
    max_dimension: int = 1 << 20
    max_alphabet: int = 32
    max_trials: int = 300_000
    max_groups: int = 256

    def __post_init__(self) -> None:
        for field in fields(self):
            value = operator.index(getattr(self, field.name))
            if value < 1:
                raise ValueError("parameter-choice budgets must be positive integers")
            object.__setattr__(self, field.name, value)
        if self.max_Q > BALANCE_MAX_Q:
            raise ValueError(f"max_Q cannot exceed {BALANCE_MAX_Q}")


@dataclass(frozen=True)
class Lemma48Proposal:
    """The exact integer schedule (4.17), whether or not it is admitted.

    When the preliminary guard fails, tau, epsilon, D and the interval stay
    unset. That is the paper's order: those objects are built only after
    r >= 2 and Q > R.
    """

    n: int
    M: int
    b: int
    w: int
    beta: int
    log_n: int
    log_log_n: int
    r: int
    Q: int
    R: int
    tau: Fraction | None
    epsilon: Fraction | None
    dimension: int | None
    interval: tuple[int, ...] | None
    dimension_formed: bool
    lemma_4_7_admitted: bool
    extras_admitted: bool
    admitted: bool
    reasons: tuple[str, ...]
    reference: ReferenceParameters | None


@dataclass(frozen=True)
class ParameterTuple:
    """One concrete (r, Q, R, tau, interval) accepted by Lemma 4.7's checker."""

    r: int
    Q: int
    R: int
    tau: Fraction
    A0: int
    interval: tuple[int, ...]
    dimension: int
    tail_probability: Fraction
    tail_limit: Fraction
    alphabet: int
    trials: int
    placeable: bool
    meets_r_at_least_32: bool
    meets_six_r_at_most_q: bool
    reference: ReferenceParameters


@dataclass(frozen=True)
class AsymptoticChoice:
    """Lemma 4.8's proposal together with the tuple selected for execution.

    ``source`` is ``lemma-4.8`` only when that closed form itself passes every
    Lemma 4.7 hypothesis and fits the budgets. ``finite-budget-member`` is a
    smaller tuple the same checker accepted. ``none`` means no tuple was
    accepted inside the budgets. ``admitted`` refers to the selected tuple.
    """

    proposal: Lemma48Proposal
    selected: ParameterTuple | None
    source: str
    admitted: bool
    placeable: bool
    examined: int
    admitted_candidates: int
    obstacles: tuple[str, ...]
    reasons: tuple[str, ...]


def _incidence(n: int, scopes: Sequence[Sequence[int]]) -> BalanceCheck:
    """Read M, b and w through the same scope normalization as the balance check."""
    return check_balance(n, scopes, 2, 1)


def lemma_4_8_proposal(
    n: int,
    scopes: Sequence[Sequence[int]],
    *,
    beta: int | None = None,
) -> Lemma48Proposal:
    """Compute Lemma 4.8's schedule and check it without allocating a matrix.

    The construction follows the proof's order. r is the greatest integer
    with r**2 * nu**4 <= lambda, Q is floor(lambda / (128 nu)), and
    R is (beta + 3) * ceil(sqrt(lambda)). Tau, epsilon, D and the centered
    interval are formed only after r >= 2 and Q > R. A failure at that guard
    does not build D = r**(2Q).
    """
    probe = _incidence(n, scopes)
    n, M, b, w = probe.n, probe.M, probe.b, probe.w
    lam = _ceil_log2(n) if n >= 1 else 0
    nu = _ceil_log2(lam) if lam >= 1 else 0
    if beta is None:
        beta = _smallest_beta(lam, b, w)
    else:
        beta = operator.index(beta)
        if beta < 1:
            raise ValueError("beta must be an integer >= 1")
    r = _largest_r(lam, nu)
    Q = lam // (128 * nu) if nu > 0 else 0
    R = (beta + 3) * _ceil_sqrt(lam) if lam >= 0 else 0
    reasons: list[str] = []
    if n < 20:
        reasons.append("Lemma 4.8 is stated for n >= 20")
    if b < 1 or w < 1 or not probe.range_ok:
        reasons.append(f"Lemma 4.7 needs n <= M <= b*n with b, w >= 1; got n={n}, M={M}, b={b}, w={w}")
    if r < 2 or Q <= R:
        reasons.append("preliminary guard failed: tau, D and the interval are formed only after r >= 2 and Q > R")
        if r < 2:
            reasons.append(f"r={r} < 2")
        if Q <= R:
            reasons.append(f"Q={Q} <= R={R}")
        return Lemma48Proposal(
            n, M, b, w, beta, lam, nu, r, Q, R, None, None, None, None,
            False, False, False, False, tuple(reasons), None,
        )
    if Q > BALANCE_MAX_Q:
        reasons.append(
            f"Q={Q} is above {BALANCE_MAX_Q}. The balance integer r**(2Q) was not "
            "formed, and no matrix was allocated"
        )
        return Lemma48Proposal(
            n, M, b, w, beta, lam, nu, r, Q, R, None, None, None, None,
            False, False, False, False, tuple(reasons), None,
        )
    tau = Fraction(1, r * lam) if lam > 0 else None
    epsilon = Fraction(1, b * r * r * lam) if lam > 0 and b > 0 else None
    if tau is None or not 0 < tau <= Fraction(1, 2):
        reasons.append("Lemma 4.8 tau = 1/(r*lambda) is outside (0, 1/2]")
        return Lemma48Proposal(
            n, M, b, w, beta, lam, nu, r, Q, R, tau, epsilon, None, None,
            False, False, False, False, tuple(reasons), None,
        )
    dimension = r ** (2 * Q)
    reference = check_reference_parameters(n, scopes, r=r, Q=Q, R=R, tau=tau)
    extras: list[str] = []
    if R < 32:
        extras.append(f"Lemma 4.8 asks for R >= 32; got R={R}")
    if 6 * R > Q:
        extras.append(f"Lemma 4.8 asks for 6R <= Q; got 6R={6 * R}, Q={Q}")
    if not _dimension_within_power(dimension, lam):
        extras.append("Lemma 4.8 asks for D <= 2**(lambda/128)")
    extras_admitted = not extras
    reasons.extend(reference.reasons)
    reasons.extend(extras)
    return Lemma48Proposal(
        n, M, b, w, beta, lam, nu, r, Q, R, tau, epsilon, dimension, reference.interval,
        True, reference.admitted, extras_admitted, reference.admitted and extras_admitted,
        tuple(reasons), reference,
    )


def _paper_tau(r: int, log_n: int, tail: Fraction, b: int) -> Fraction | None:
    """Least tau in (0, 1/2] meeting the exact tail, if any such tau exists.

    Lemma 4.8's own value 1/(r lambda) is used when the tail is zero, since
    every positive tau then works and that is the reference's choice. Otherwise
    the bound tau >= tail * r * b is tight, and equality is allowed.
    """
    if tail == 0:
        if log_n <= 0:
            return Fraction(1, 2)
        candidate = Fraction(1, r * log_n)
        return candidate if candidate <= Fraction(1, 2) else Fraction(1, 2)
    if b <= 0:
        return None
    candidate = tail * r * b
    if 0 < candidate <= Fraction(1, 2):
        return candidate
    return None


def _rank(candidate: ParameterTuple) -> tuple:
    """Prefer the measurement term, then R/Q, then 1/r. Ties break by D, r, Q, R.

    These are the three summands in Lemma 4.8's exponent, compared exactly.
    The comparison does not prove the O() bound attached to that exponent.
    """
    return (
        candidate.tau,
        Fraction(candidate.R, candidate.Q),
        Fraction(1, candidate.r),
        candidate.dimension,
        candidate.r,
        candidate.Q,
        candidate.R,
    )


def _placeable(r: int, Q: int, n: int, factor_count: int, budgets: ChoiceBudgets) -> tuple[bool, int, int]:
    alphabet = r ** Q
    trials = (n + factor_count) * alphabet
    ok = (
        r <= budgets.max_groups
        and Q <= budgets.max_Q
        and alphabet <= budgets.max_alphabet
        and trials <= budgets.max_trials
    )
    return ok, alphabet, trials


def choose_asymptotic_parameters(
    n: int,
    scopes: Sequence[Sequence[int]],
    *,
    beta: int | None = None,
    budgets: ChoiceBudgets | None = None,
    factor_count: int | None = None,
) -> AsymptoticChoice:
    """Propose Lemma 4.8's tuple and, inside the budgets, one Lemma 4.7 tuple.

    The finite search enumerates r, Q and R inside ``ChoiceBudgets`` and keeps
    a candidate only after ``check_reference_parameters`` accepts it. The
    centered interval is the one in the proof of Lemma 4.8. Tau is the least
    value in (0, 1/2] that meets the exact binomial tail. No matrix is
    allocated and no schedule search is started.
    """
    budgets = budgets or ChoiceBudgets()
    proposal = lemma_4_8_proposal(n, scopes, beta=beta)
    if factor_count is None:
        factor_count = len(scopes)
    else:
        factor_count = operator.index(factor_count)
        if factor_count < 0:
            raise ValueError("factor_count must be nonnegative")
    examined = 0
    obstacles: list[str] = []
    admitted_rows: list[ParameterTuple] = []
    placeable_rows: list[ParameterTuple] = []

    def obstacle(text: str) -> None:
        if len(obstacles) < 8 and text not in obstacles:
            obstacles.append(text)

    searchable = proposal.b >= 1 and proposal.w >= 1 and proposal.n <= proposal.M <= proposal.b * proposal.n
    if not searchable:
        obstacle("the scope totals do not meet n <= M <= b*n with b, w >= 1")
    else:
        for r in range(2, budgets.max_r + 1):
            for Q in range(2, budgets.max_Q + 1):
                dimension = r ** (2 * Q)
                if dimension > budgets.max_dimension:
                    obstacle(f"r={r}, Q={Q}: D={dimension} exceeds the choice dimension budget {budgets.max_dimension}")
                    continue
                variable_bound = 2 * Q * r * r
                if variable_bound > proposal.n:
                    obstacle(f"r={r}, Q={Q}: 2Qr^2 = {variable_bound} > n = {proposal.n}")
                    continue
                mass_bound = 2 * (proposal.b + proposal.w) * dimension * dimension
                if mass_bound > proposal.M:
                    obstacle(
                        f"r={r}, Q={Q}: 2(b+w)D^2 = {mass_bound} > M = {proposal.M}"
                    )
                    continue
                fits, alphabet, trials = _placeable(r, Q, proposal.n, factor_count, budgets)
                for R in range(1, Q):
                    examined += 1
                    A0 = default_interval_start(r, Q, R)
                    if A0 < 0 or A0 + R - 1 > Q:
                        obstacle(f"r={r}, Q={Q}, R={R}: centered interval is not inside {{0, ..., Q}}")
                        continue
                    interval = tuple(range(A0, A0 + R))
                    tail = binomial_tail_probability(r, Q, interval)
                    tau = _paper_tau(r, proposal.log_n, tail, proposal.b)
                    if tau is None:
                        obstacle(
                            f"r={r}, Q={Q}, R={R}: binomial tail {tail} exceeds every tau <= 1/2"
                        )
                        continue
                    reference = check_reference_parameters(
                        proposal.n, scopes, r=r, Q=Q, R=R, tau=tau, A0=A0,
                    )
                    if not reference.admitted:
                        obstacle(f"r={r}, Q={Q}, R={R}: {reference.reasons[0]}")
                        continue
                    candidate = ParameterTuple(
                        r, Q, R, tau, A0, interval, dimension, reference.tail_probability,
                        reference.tail_limit if reference.tail_limit is not None else tau,
                        alphabet, trials, fits, R >= 32, 6 * R <= Q, reference,
                    )
                    admitted_rows.append(candidate)
                    if fits:
                        placeable_rows.append(candidate)

    selected: ParameterTuple | None = None
    source = "none"
    if proposal.reference is not None and proposal.lemma_4_7_admitted and proposal.dimension is not None and proposal.tau is not None and proposal.interval is not None:
        proposal_placeable, alphabet, trials = _placeable(
            proposal.r, proposal.Q, proposal.n, factor_count, budgets,
        )
        proposal_placeable = proposal_placeable and proposal.dimension <= budgets.max_dimension
        if proposal_placeable:
            selected = ParameterTuple(
                proposal.r, proposal.Q, proposal.R, proposal.tau, proposal.interval[0],
                proposal.interval, proposal.dimension, proposal.reference.tail_probability,
                proposal.reference.tail_limit if proposal.reference.tail_limit is not None else proposal.tau,
                alphabet, trials, True, proposal.R >= 32, 6 * proposal.R <= proposal.Q,
                proposal.reference,
            )
            source = "lemma-4.8"
    if selected is None and placeable_rows:
        selected = min(placeable_rows, key=_rank)
        source = "finite-budget-member"
    elif selected is None and admitted_rows:
        selected = min(admitted_rows, key=_rank)
        source = "finite-budget-member"

    reasons: list[str] = []
    if selected is None:
        reasons.append("no budget-bounded tuple meets every Lemma 4.7 numerical hypothesis")
    elif not selected.placeable:
        reasons.append(
            "Lemma 4.7 admits this tuple, but the sparse-placement budgets refuse it "
            "before any schedule search or matrix allocation"
        )
        reasons.append(f"alphabet {selected.alphabet} and candidate trials {selected.trials} were not started")
    return AsymptoticChoice(
        proposal, selected, source, selected is not None, bool(selected and selected.placeable),
        examined, len(admitted_rows), tuple(obstacles), tuple(reasons),
    )


def choose_asymptotic_parameters_for_cnf(cnf, **kwargs) -> AsymptoticChoice:
    """Choose from a CNF. Factor scopes are the clause variable sets."""
    from .scheduler import scopes_from_cnf

    return choose_asymptotic_parameters(
        cnf.n, scopes_from_cnf(cnf), factor_count=cnf.m, **kwargs,
    )
