"""
Balance checks and exact conditional-expectation schedules (Steps 7–9).

Lemma 4.5 of the October 2026 preprint asks for schedules whose label
multiplicity is at most 2M/D and whose groups have size at most 2n/r, provided

    2(b+w) D^2 <= M,    2 Q r^2 <= n,

with D = r**(2Q), n <= M <= b n, b the maximum number of scopes containing a
variable, and w the maximum scope size. Together with M <= b n those two
inequalities force n to be at least on the order of D^2: when b >= 1,

    n >= 2 (b+w) D^2 / b >= 2 D^2.

Lemma 4.7 adds 1 <= R < Q and a binomial tail bound. This module checks the
numerical balance inequalities and the degree gap. It does not check the tail
bound.

When the inequalities fail, or when the slice would sit above the cap, a call
that does not supply schedules refuses and does not sample a placement. When
the inequalities hold and the slice is under the cap, Step 8 runs a
successive-choice search over the schedule alphabet {0, ..., r-1}^Q. At each
variable and each clause it tries every alphabet symbol and keeps the one that
minimizes a quadratic load potential. Ties break toward the lexicographically
first symbol. Nothing is sampled. The schedule cube, of size (r^Q)^(n+m), is
not enumerated.

The search minimizes the exact conditional expectation of Phi/(M(b+w))
+ Xi/(Qn), including dependencies of incidences sharing variables or factors.
Its potential is checked with rational arithmetic at the initial and final
states, and each choice must preserve monotonicity. The default design check
uses a cached rational row span, so there is no Boolean-point walk. The
legacy ``design_method="boolean"`` enables that extra check and retains the
Step 8 Boolean budget. Placement coefficients have their own allocation
budget. A supplied R that fails 1 <= R < Q blocks balancer acceptance.

Caller-supplied schedules are placed with ``place_formula`` and checked with
``design_problems_at`` and the T = 0 restoration. That is how the hand-written
Step 3 placements, and the two Step 7 slices, are accepted. A supplied
placement is not presented as the output of a balancer. The multiplicity bound
2M/D is recorded on that placement. It is not required for acceptance of a
supplied placement, because the hand-written rows can exceed it precisely when
the balance hypotheses fail. The group-size bound is already inside
``design_problems_at``.

Nothing here is Theorem 1.1. This module does not refute ETH. Nothing it runs
is a 2^{o(n)} count. Refusals are not timed. A passing balance check is not,
by itself, a schedule, and a schedule this search returns is not a count.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .cnf import CNF
from .conditional_balance import choose_schedules
from .replacement_matrix import (
    FULL_M0_LOG2_DIMENSION_AT_R2,
    FULL_M0_MIN_Q,
    FULL_M0_MIN_R,
    SLICE_MAX_DIMENSION,
    AuxiliarySlice,
    PlacedFormula,
    auxiliary_problems,
    build_auxiliary_slice,
    degree_cap_report,
    binomial_exceptional_count,
    default_interval_start,
    design_problems_at,
    groups_at,
    place_formula,
    restoration_problems,
)

# Above this, the balance integer r**(2Q) is not formed. Full M0 uses Q = 192.
BALANCE_MAX_Q = 4096

# Successive choice tries each alphabet symbol once per variable and per clause.
# The schedule cube, size (r**Q)**(n + m), is never enumerated.
BALANCER_MAX_ALPHABET = 32
BALANCER_MAX_TRIALS = 100_000
# ``design_problems_at`` walks 2**(group size) Boolean points on an outside group.
BALANCER_DESIGN_BOOLEAN_BUDGET = 4096
BALANCER_PLACEMENT_ENTRY_BUDGET = 4_000_000


def _as_index(value: object, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{what} must be an int, got {value!r}")
    return value


def _fmt(value: int) -> str:
    """Decimal for small integers, an exact power-of-two form for large ones."""
    if value < 0:
        return "-" + _fmt(-value)
    if value.bit_length() <= 40:
        return str(value)
    if value > 0 and value & (value - 1) == 0:
        return f"2**{value.bit_length() - 1}"
    shift = (value & -value).bit_length() - 1
    odd = value >> shift
    if shift and odd.bit_length() <= 40:
        return f"{odd}*2**{shift}"
    return f"integer with bit length {value.bit_length()}"


def scopes_from_cnf(cnf: CNF) -> tuple[tuple[int, ...], ...]:
    """Variable sets read by the clauses, in clause order.

    A repeated variable inside one clause counts once. An empty clause is an
    empty scope. This is the factor-scope list the balance inequalities read.
    """
    scopes: list[tuple[int, ...]] = []
    for clause in cnf.clauses:
        seen: list[int] = []
        present: set[int] = set()
        for var, _neg in clause:
            if var in present:
                continue
            present.add(var)
            seen.append(var)
        scopes.append(tuple(seen))
    return tuple(scopes)


def _normalize_scopes(
    n: int,
    scopes: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    n = _as_index(n, "n")
    if n < 0:
        raise ValueError(f"n must be nonnegative, got {n}")
    parsed: list[tuple[int, ...]] = []
    for index, scope in enumerate(scopes):
        seen: list[int] = []
        present: set[int] = set()
        for entry in scope:
            var = _as_index(entry, f"scope[{index}]")
            if not (0 <= var < n):
                raise ValueError(f"scope[{index}] variable {var} is outside {{0, ..., {n - 1}}}")
            if var in present:
                continue
            present.add(var)
            seen.append(var)
        parsed.append(tuple(seen))
    return tuple(parsed)


@dataclass(frozen=True)
class BalanceCheck:
    """Numerical hypotheses of the Lemma 4.5 / 4.7 balance package.

    ``admitted`` is the Lemma 4.5 package: valid parameters, n <= M <= b n,
    2(b+w)D^2 <= M, and 2Qr^2 <= n. ``degree_gap_ok`` is the Lemma 4.7
    condition 1 <= R < Q when R was supplied. The binomial tail is not a field
    here. ``admitted`` does not mean a balancer ran.
    """

    n: int
    M: int
    b: int
    w: int
    r: int
    Q: int
    R: int | None
    D: int
    mass_bound: int
    variable_bound: int
    minimum_n: int | None
    range_ok: bool
    mass_ok: bool
    variable_ok: bool
    admitted: bool
    degree_gap_ok: bool | None
    reasons: tuple[str, ...]

    @property
    def order_bound(self) -> int:
        """The weak necessary size 2 D^2. Admitted instances with b >= 1 meet it."""
        return 2 * self.D * self.D


def check_balance(
    n: int,
    scopes: Sequence[Sequence[int]],
    r: int,
    Q: int,
    R: int | None = None,
) -> BalanceCheck:
    """Check the balance inequalities. Does not allocate a matrix or choose schedules.

    ``scopes`` are the factor variable-sets. A CNF is converted with
    ``scopes_from_cnf``. M is the sum of the scope sizes, w is the largest
    scope, and b is the largest number of scopes containing one variable.
    """
    parsed = _normalize_scopes(n, scopes)
    n = _as_index(n, "n")
    r = _as_index(r, "r")
    Q = _as_index(Q, "Q")
    if r < 2:
        raise ValueError(f"r must be an integer ≥ 2, got {r}")
    if Q < 1:
        raise ValueError(f"Q must be a positive integer, got {Q}")
    if Q > BALANCE_MAX_Q:
        raise ValueError(
            f"Q={Q} is above {BALANCE_MAX_Q}. The balance integer r**(2Q) was not "
            "formed, and no matrix was allocated"
        )
    if R is not None:
        R = _as_index(R, "R")
        if R < 1:
            raise ValueError(f"R must be a positive integer, got {R}")
    counts = [0 for _ in range(n)]
    for scope in parsed:
        for var in scope:
            counts[var] += 1
    M = sum(len(scope) for scope in parsed)
    w = max((len(scope) for scope in parsed), default=0)
    b = max(counts, default=0)
    D = r ** (2 * Q)
    mass_bound = 2 * (b + w) * D * D
    variable_bound = 2 * Q * r * r
    if b >= 1:
        minimum_n: int | None = (mass_bound + b - 1) // b
    else:
        minimum_n = None
    range_ok = n <= M <= b * n
    mass_ok = mass_bound <= M
    variable_ok = variable_bound <= n
    admitted = range_ok and mass_ok and variable_ok
    if R is None:
        degree_gap_ok: bool | None = None
    else:
        degree_gap_ok = 1 <= R < Q

    reasons: list[str] = []
    if not range_ok:
        reasons.append(f"need n <= M <= b*n, got n={n}, M={M}, b={b}")
    if not mass_ok:
        tighter = ""
        if minimum_n is not None:
            tighter = f"; with M <= b*n this needs n >= { _fmt(minimum_n) }"
        reasons.append(
            f"2(b+w)D^2 = {_fmt(mass_bound)} > M = {M} "
            f"(b={b}, w={w}, D={_fmt(D)}){tighter}"
        )
    if not variable_ok:
        reasons.append(
            f"2Qr^2 = {variable_bound} > n = {n} (Q={Q}, r={r})"
        )
    if b >= 1 and n < 2 * D * D:
        reasons.append(
            "these inequalities need n on the order of D^2; "
            f"n = {n} < 2*D^2 = {_fmt(2 * D * D)}"
        )
    if degree_gap_ok is False:
        reasons.append(f"Lemma 4.7 asks for 1 <= R < Q; got R={R}, Q={Q}")
    return BalanceCheck(
        n=n,
        M=M,
        b=b,
        w=w,
        r=r,
        Q=Q,
        R=R,
        D=D,
        mass_bound=mass_bound,
        variable_bound=variable_bound,
        minimum_n=minimum_n,
        range_ok=range_ok,
        mass_ok=mass_ok,
        variable_ok=variable_ok,
        admitted=admitted,
        degree_gap_ok=degree_gap_ok,
        reasons=tuple(reasons),
    )


def _allocation_refused(check: BalanceCheck) -> bool:
    return check.Q > 30 or check.D > SLICE_MAX_DIMENSION


def _allocation_sentence(check: BalanceCheck) -> str:
    if not _allocation_refused(check):
        return ""
    return (
        f"Side length D={_fmt(check.D)} for r={check.r}, Q={check.Q} is above the "
        f"slice cap {SLICE_MAX_DIMENSION}. build_auxiliary_slice was not called; "
        f"the matrix was not allocated. Full M0 at R>={FULL_M0_MIN_R} and "
        f"Q>={FULL_M0_MIN_Q} has side length 2**{FULL_M0_LOG2_DIMENSION_AT_R2} at r=2. "
    )


def _gap_sentence(check: BalanceCheck) -> str:
    """Lemma 4.7's degree gap, when this call cannot be a witness for it."""
    if check.degree_gap_ok is False:
        return f" Lemma 4.7 asks for 1 <= R < Q; R={check.R}, Q={check.Q}."
    if check.R is None and check.Q <= 1:
        return (
            f" No integer R satisfies 1 <= R < Q at Q={check.Q}, so this call "
            "cannot be a Lemma 4.7 witness."
        )
    return ""


def _degree_gap_blocks_acceptance(check: BalanceCheck) -> bool:
    """Acceptance needs 1 <= R < Q when a slice would be built at a concrete R.

    A supplied R that fails the gap blocks acceptance. Q <= 1 blocks it as
    well: every positive R is then at least Q, including the R = 1 slice this
    module would build when the caller omitted R.
    """
    if check.degree_gap_ok is False:
        return True
    return check.Q <= 1 and check.degree_gap_ok is not True


def _unplaced_reason(check: BalanceCheck) -> str:
    allocation = _allocation_sentence(check)
    if check.admitted:
        return (
            allocation
            + "The Lemma 4.5 balance inequalities hold: "
            + f"2(b+w)D^2 = {_fmt(check.mass_bound)} <= M = {check.M} and "
            + f"2Qr^2 = {check.variable_bound} <= n = {check.n}. "
            + "n is on the order of D^2. "
            + "The balancer was not run, because the matrix was not allocated."
            + _gap_sentence(check)
        )
    failed = "; ".join(check.reasons)
    return (
        allocation
        + "Refusing to choose schedules: the Lemma 4.5 balance inequalities fail "
        + f"({failed}). No placement was invented. The balancer was not run."
    )


@dataclass(frozen=True)
class ScheduleOutcome:
    """Result of a balance check, plus an optional caller-supplied placement.

    ``accepted`` means restoration and design checks passed. For balancer
    output it additionally requires the degree gap and load bounds. It does
    not certify Lemma 4.7's tail or Theorem 1.1. ``design_method`` records the
    checker used. ``design_boolean_points`` is a hypothetical legacy walk,
    not an executed count. ``built_slice`` tracks allocation by this call.
    """

    balance: BalanceCheck
    placed: PlacedFormula | None
    design_problems: tuple[tuple[int, str], ...]
    restoration_problems: tuple[str, ...]
    accepted: bool
    caller_supplied: bool
    built_slice: bool
    allocation_refused: bool
    multiplicity: int | None
    multiplicity_bound: Fraction | None
    multiplicity_ok: bool | None
    group_size_ok: bool | None
    reason: str
    balancer_attempted: bool = False
    balancer_built: bool = False
    design_checked: bool = False
    max_group_size: int | None = None
    outside_group_size: int | None = None
    design_boolean_points: int | None = None
    trials: int | None = None
    design_method: str | None = None
    initial_potential: Fraction | None = None
    final_potential: Fraction | None = None


def _blank_outcome(check: BalanceCheck, *, caller_supplied: bool) -> ScheduleOutcome:
    return ScheduleOutcome(
        balance=check,
        placed=None,
        design_problems=(),
        restoration_problems=(),
        accepted=False,
        caller_supplied=caller_supplied,
        built_slice=False,
        allocation_refused=_allocation_refused(check),
        multiplicity=None,
        multiplicity_bound=Fraction(2 * check.M, check.D),
        multiplicity_ok=None,
        group_size_ok=None,
        reason=_unplaced_reason(check),
    )


def _alphabet_size(r: int, Q: int) -> int:
    return r ** Q


def _search_block_reason(n: int, m: int, r: int, Q: int) -> str | None:
    """Refuse before the successive loop when the alphabet or the trial count is too large.

    Returning a sentence means the loop did not run. The schedule cube is not
    mentioned as something that was enumerated; it was not.
    """
    alphabet = _alphabet_size(r, Q)
    if alphabet > BALANCER_MAX_ALPHABET:
        return (
            f"The balancer was not run: the schedule alphabet has size r**Q = {alphabet}, "
            f"above {BALANCER_MAX_ALPHABET}. That alphabet was not enumerated, and the "
            "schedule cube was not enumerated. No placement was invented."
        )
    trials = (n + m) * alphabet
    if trials > BALANCER_MAX_TRIALS:
        return (
            "The balancer was not run: the successive search would consider "
            f"(n+m)*r**Q = {trials} candidate symbols, above the budget "
            f"{BALANCER_MAX_TRIALS}. The schedule cube was not enumerated. "
            "No placement was invented."
        )
    return None


@dataclass(frozen=True)
class _SearchResult:
    """Schedules from successive choice, or a refusal that did not place them."""

    attempted: bool
    chosen: bool
    variable_schedules: tuple[tuple[int, ...], ...] | None
    clause_schedules: tuple[tuple[int, ...], ...] | None
    multiplicity: int | None
    multiplicity_bound: Fraction
    multiplicity_ok: bool | None
    max_group_size: int | None
    group_size_ok: bool | None
    outside_group_size: int | None
    design_boolean_points: int | None
    trials: int | None
    reason: str
    initial_potential: Fraction | None = None
    final_potential: Fraction | None = None


def _search_schedules(cnf: CNF, check: BalanceCheck) -> _SearchResult:
    """Choose the exact conditional expectation from reference Lemma 4.5.

    ``choose_schedules`` aggregates collision probabilities and checks the
    potential's monotonicity and final value with exact rational arithmetic.
    The alphabet and trial budgets are checked before the search runs.
    """
    bound = Fraction(2 * check.M, check.D)
    blocked = _search_block_reason(cnf.n, cnf.m, check.r, check.Q)
    if blocked is not None:
        return _SearchResult(
            attempted=False,
            chosen=False,
            variable_schedules=None,
            clause_schedules=None,
            multiplicity=None,
            multiplicity_bound=bound,
            multiplicity_ok=None,
            max_group_size=None,
            group_size_ok=None,
            outside_group_size=None,
            design_boolean_points=None,
            trials=None,
            reason=blocked,
        )
    schedules = choose_schedules(cnf.n, scopes_from_cnf(cnf), check.r, check.Q)
    variable_schedules = schedules.variables
    clause_schedules = schedules.factors
    group_counts = [[0 for _ in range(check.r)] for _ in range(check.Q)]
    for sigma in variable_schedules:
        for t, symbol in enumerate(sigma):
            group_counts[t][symbol] += 1
    label_counts: dict[tuple[tuple[int, ...], tuple[int, ...]], int] = {}
    for clause, tau in zip(cnf.clauses, clause_schedules):
        for var, _neg in clause:
            label = (variable_schedules[var], tau)
            label_counts[label] = label_counts.get(label, 0) + 1

    multiplicity = max(label_counts.values(), default=0)
    multiplicity_ok = multiplicity <= bound
    max_group = max((max(row) for row in group_counts), default=0)
    group_ok = all(size * check.r <= 2 * check.n for row in group_counts for size in row)
    outside, points = _outside_load(
        cnf,
        variable_schedules,
        clause_schedules,
        check.r,
        check.Q,
    )
    trials = schedules.trials
    if not multiplicity_ok or not group_ok:
        reason = _bound_failure_reason(
            check,
            multiplicity=multiplicity,
            bound=bound,
            multiplicity_ok=multiplicity_ok,
            max_group=max_group,
            group_ok=group_ok,
        )
        return _SearchResult(
            attempted=True,
            chosen=False,
            variable_schedules=tuple(variable_schedules),
            clause_schedules=tuple(clause_schedules),
            multiplicity=multiplicity,
            multiplicity_bound=bound,
            multiplicity_ok=multiplicity_ok,
            max_group_size=max_group,
            group_size_ok=group_ok,
            outside_group_size=outside,
            design_boolean_points=points,
            trials=trials,
            reason=reason,
            initial_potential=schedules.initial_potential,
            final_potential=schedules.final_potential,
        )
    return _SearchResult(
        attempted=True,
        chosen=True,
        variable_schedules=tuple(variable_schedules),
        clause_schedules=tuple(clause_schedules),
        multiplicity=multiplicity,
        multiplicity_bound=bound,
        multiplicity_ok=True,
        max_group_size=max_group,
        group_size_ok=True,
        outside_group_size=outside,
        design_boolean_points=points,
        trials=trials,
        reason="",
        initial_potential=schedules.initial_potential,
        final_potential=schedules.final_potential,
    )


def _outside_load(
    cnf: CNF,
    variable_schedules: Sequence[Sequence[int]],
    clause_schedules: Sequence[Sequence[int]],
    r: int,
    Q: int,
) -> tuple[int, int | None]:
    """Largest outside group, and the Boolean-point total when it is small enough to form.

    The total matches the walk in ``design_problems_at``: one 2**(group size)
    product for each occurrence and each group other than the clause's symbol,
    at each t. A size above 20 is already past the budget; the total is then
    left unformed.
    """
    sizes = [[0 for _ in range(r)] for _ in range(Q)]
    for schedule in variable_schedules:
        for t, symbol in enumerate(schedule):
            sizes[t][symbol] += 1
    outside = 0
    points = 0
    formed = True
    for clause, clause_schedule in zip(cnf.clauses, clause_schedules):
        width = len(clause)
        if width == 0:
            continue
        for t, assigned in enumerate(clause_schedule):
            for symbol in range(r):
                if symbol == assigned:
                    continue
                size = sizes[t][symbol]
                if size > outside:
                    outside = size
                if size > 20:
                    formed = False
                    continue
                points += width * (1 << size)
    return outside, points if formed else None


def _bound_failure_reason(
    check: BalanceCheck,
    *,
    multiplicity: int,
    bound: Fraction,
    multiplicity_ok: bool,
    max_group: int,
    group_ok: bool,
) -> str:
    sentences = [
        "The successive-choice balancer stopped before place_formula. "
        "The schedule cube was not enumerated."
    ]
    if multiplicity_ok:
        sentences.append(
            f"Label multiplicity {multiplicity} is at most 2M/D = {bound}."
        )
    else:
        sentences.append(
            f"Label multiplicity {multiplicity} exceeds 2M/D = {bound}. "
            "No placement was invented."
        )
    if group_ok:
        sentences.append(
            f"The largest group has size {max_group}, and "
            f"{max_group} * r = {max_group * check.r} <= 2n = {2 * check.n}."
        )
    else:
        sentences.append(
            f"The largest group has size {max_group}, and "
            f"{max_group} * r = {max_group * check.r} > 2n = {2 * check.n}. "
            "No placement was invented."
        )
    return " ".join(sentences)


def _design_budget_reason(check: BalanceCheck, search: _SearchResult) -> str:
    outside = search.outside_group_size if search.outside_group_size is not None else 0
    points = (
        f" The Boolean-point total would be {search.design_boolean_points}."
        if search.design_boolean_points is not None
        else ""
    )
    assert search.max_group_size is not None
    trials = f" The search considered {search.trials} candidate symbols." if search.trials else ""
    return (
        "The successive-choice balancer chose schedules over the alphabet "
        f"{{0, ..., {check.r - 1}}}^{check.Q} "
        f"({_alphabet_size(check.r, check.Q)} symbols). "
        "The schedule cube was not enumerated."
        + trials
        + f" Label multiplicity {search.multiplicity} is at most "
        f"2M/D = {search.multiplicity_bound}. "
        f"The largest group has size {search.max_group_size}, and "
        f"{search.max_group_size} * r = {search.max_group_size * check.r} "
        f"<= 2n = {2 * check.n}. "
        "design_problems_at was not called: an outside group has size "
        f"{outside}, so that call enumerates 2**{outside} Boolean points for "
        f"that group on each occurrence that reads it, above the budget "
        f"{BALANCER_DESIGN_BOOLEAN_BUDGET}.{points} "
        "place_formula was not called. No placement was accepted."
        + _gap_sentence(check)
    )


def _search_outcome(
    check: BalanceCheck, search: _SearchResult, *, reason: str | None = None
) -> ScheduleOutcome:
    if reason is None and search.chosen:
        reason = _design_budget_reason(check, search)
    elif reason is None:
        reason = search.reason + _gap_sentence(check)
    return ScheduleOutcome(
        balance=check,
        placed=None,
        design_problems=(),
        restoration_problems=(),
        accepted=False,
        caller_supplied=False,
        built_slice=False,
        allocation_refused=False,
        multiplicity=search.multiplicity,
        multiplicity_bound=search.multiplicity_bound,
        multiplicity_ok=search.multiplicity_ok,
        group_size_ok=search.group_size_ok,
        reason=reason,
        balancer_attempted=search.attempted,
        balancer_built=False,
        design_checked=False,
        max_group_size=search.max_group_size,
        outside_group_size=search.outside_group_size,
        design_boolean_points=search.design_boolean_points,
        trials=search.trials,
        initial_potential=search.initial_potential,
        final_potential=search.final_potential,
    )


def _design_budget_exceeded(search: _SearchResult) -> bool:
    outside = search.outside_group_size or 0
    if outside > 20:
        return True
    if search.design_boolean_points is None:
        return True
    return search.design_boolean_points > BALANCER_DESIGN_BOOLEAN_BUDGET


def _finish_chosen_placement(
    cnf: CNF,
    check: BalanceCheck,
    search: _SearchResult,
    auxiliary: AuxiliarySlice | None,
    *,
    design_method: str = "linear",
) -> ScheduleOutcome:
    """Place chosen schedules and check their rational kernels and restoration.

    The Boolean comparison is optional and is budgeted by the calling gate.
    Dense placement entries are budgeted before a slice is allocated.
    """
    assert search.variable_schedules is not None and search.clause_schedules is not None
    R = auxiliary.R if auxiliary is not None else (check.R if check.R is not None else 1)
    measurements = (len(auxiliary.C_rows) if auxiliary is not None else
                    binomial_exceptional_count(check.r, check.Q,
                        range(default_interval_start(check.r, check.Q, R),
                              default_interval_start(check.r, check.Q, R) + R)))
    entries = cnf.n * (check.M * R + (search.multiplicity or 0) * measurements)
    if entries > BALANCER_PLACEMENT_ENTRY_BUDGET:
        return _search_outcome(
            check, search,
            reason=(f"Placement would store at most {entries} coefficient entries, above "
                    f"the budget {BALANCER_PLACEMENT_ENTRY_BUDGET}. place_formula was not called."),
        )
    built = auxiliary is None
    if auxiliary is None:
        auxiliary = build_auxiliary_slice(check.r, check.Q, R)
    placed = place_formula(
        cnf,
        auxiliary,
        search.variable_schedules,
        search.clause_schedules,
    )
    design: list[tuple[int, str]] = []
    for t in range(1, auxiliary.Q + 1):
        for problem in design_problems_at(placed, t, boolean_check=design_method == "boolean"):
            design.append((t, problem))
    restored = tuple(restoration_problems(placed))
    multiplicity = placed.copies
    bound = Fraction(2 * check.M, check.D)
    multiplicity_ok = multiplicity <= bound
    group_ok = _group_size_ok(placed)
    accepted = not design and not restored and multiplicity_ok and group_ok
    gap_blocked = _degree_gap_blocks_acceptance(check)
    if gap_blocked:
        accepted = False
    if accepted:
        design_sentence = (
            "design_problems_at and T=0 restoration passed at every t in {1, ..., Q}. "
        )
    else:
        if design or restored:
            detail = design[0][1] if design else restored[0]
            design_sentence = (
                "design_problems_at or T=0 restoration reported a problem: "
                f"{detail} "
            )
        elif gap_blocked:
            design_sentence = (
                "design_problems_at and T=0 restoration passed. "
                "Acceptance is blocked because the Lemma 4.7 degree gap does not hold. "
            )
        else:
            design_sentence = "design_problems_at or T=0 restoration reported a problem. "
    if multiplicity_ok:
        multiplicity_sentence = f"Label multiplicity {multiplicity} is at most 2M/D = {bound}. "
    else:
        multiplicity_sentence = f"Label multiplicity {multiplicity} exceeds 2M/D = {bound}. "
    if group_ok:
        group_sentence = "Every group has size at most 2n/r. "
    else:
        group_sentence = "Some group has size above 2n/r. "
    balance_sentence = (
        "The Lemma 4.5 balance inequalities hold. "
        if check.admitted
        else "The Lemma 4.5 balance inequalities were not the gate for this call. "
    )
    slice_sentence = ""
    if built and check.R is None:
        slice_sentence = f"The slice was built at R={auxiliary.R}. "
    return ScheduleOutcome(
        balance=check,
        placed=placed,
        design_problems=tuple(design),
        restoration_problems=restored,
        accepted=accepted,
        caller_supplied=False,
        built_slice=built,
        allocation_refused=False,
        multiplicity=multiplicity,
        multiplicity_bound=bound,
        multiplicity_ok=multiplicity_ok,
        group_size_ok=group_ok,
        reason=(
            "Schedules were chosen by the successive-choice balancer. "
            "They were not supplied by the caller. The schedule cube was not enumerated. "
            + balance_sentence
            + slice_sentence
            + design_sentence
            + group_sentence
            + multiplicity_sentence
            + _gap_sentence(check).lstrip()
        ).rstrip(),
        balancer_attempted=True,
        balancer_built=True,
        design_checked=True,
        max_group_size=search.max_group_size,
        outside_group_size=search.outside_group_size,
        design_boolean_points=search.design_boolean_points,
        trials=search.trials,
        design_method=design_method,
        initial_potential=search.initial_potential,
        final_potential=search.final_potential,
    )


def _group_size_ok(placed: PlacedFormula) -> bool:
    n = placed.cnf.n
    r = placed.auxiliary.r
    for t in range(1, placed.auxiliary.Q + 1):
        for group in groups_at(placed, t):
            if len(group) * r > 2 * n:
                return False
    return True


def schedule_placement(
    cnf: CNF,
    r: int | None = None,
    Q: int | None = None,
    R: int | None = None,
    auxiliary: AuxiliarySlice | None = None,
    variable_schedules: Sequence[Sequence[int]] | None = None,
    clause_schedules: Sequence[Sequence[int]] | None = None,
    *,
    design_method: str = "linear",
) -> ScheduleOutcome:
    """Check balance, then place supplied schedules or run the successive-choice search.

    Parameters come from ``auxiliary`` when it is given. Explicit ``r``, ``Q``,
    and ``R`` must match that slice. With no schedules, the search runs only
    when the Lemma 4.5 inequalities hold and the slice is under the cap. It
    tries one alphabet symbol at a time. It does not draw a placement, and it
    does not enumerate the schedule cube. A failed bound or an exhausted
    budget returns no placement. Linear design checks are the default; the
    Boolean comparison is an optional legacy mode. With schedules, a slice above the cap is
    refused before ``build_auxiliary_slice``. A slice at or under the cap is
    built only when the caller did not pass one.
    """
    if design_method not in ("linear", "boolean"):
        raise ValueError("design_method must be 'linear' or 'boolean'")
    supplied_var = variable_schedules is not None
    supplied_clause = clause_schedules is not None
    if supplied_var != supplied_clause:
        raise ValueError("variable schedules and clause schedules must be supplied together")
    if auxiliary is not None:
        if r is not None and r != auxiliary.r:
            raise ValueError(f"r={r} does not match the auxiliary slice r={auxiliary.r}")
        if Q is not None and Q != auxiliary.Q:
            raise ValueError(f"Q={Q} does not match the auxiliary slice Q={auxiliary.Q}")
        if R is not None and R != auxiliary.R:
            raise ValueError(f"R={R} does not match the auxiliary slice R={auxiliary.R}")
        r, Q, R = auxiliary.r, auxiliary.Q, auxiliary.R
    if r is None or Q is None:
        raise ValueError("pass an auxiliary slice, or pass r and Q")
    check = check_balance(cnf.n, scopes_from_cnf(cnf), r, Q, R)
    if not supplied_var:
        if not check.admitted or _allocation_refused(check):
            return _blank_outcome(check, caller_supplied=False)
        search = _search_schedules(cnf, check)
        if not search.chosen or (design_method == "boolean" and _design_budget_exceeded(search)):
            return _search_outcome(check, search)
        return _finish_chosen_placement(cnf, check, search, auxiliary, design_method=design_method)
    if auxiliary is None and R is None:
        raise ValueError("pass R when the schedules need a new slice, or pass an auxiliary slice")
    if auxiliary is None and _allocation_refused(check):
        outcome = _blank_outcome(check, caller_supplied=True)
        return ScheduleOutcome(
            balance=outcome.balance,
            placed=None,
            design_problems=(),
            restoration_problems=(),
            accepted=False,
            caller_supplied=True,
            built_slice=False,
            allocation_refused=True,
            multiplicity=None,
            multiplicity_bound=outcome.multiplicity_bound,
            multiplicity_ok=None,
            group_size_ok=None,
            reason=(
                outcome.reason
                + " Caller-supplied schedules were not placed, because the matrix "
                + "was not allocated."
            ),
        )
    built = False
    if auxiliary is None:
        assert R is not None
        auxiliary = build_auxiliary_slice(r, Q, R)
        built = True
    assert variable_schedules is not None and clause_schedules is not None
    placed = place_formula(cnf, auxiliary, variable_schedules, clause_schedules)
    design: list[tuple[int, str]] = []
    for t in range(1, auxiliary.Q + 1):
        for problem in design_problems_at(placed, t, boolean_check=design_method == "boolean"):
            design.append((t, problem))
    restored = tuple(restoration_problems(placed))
    multiplicity = placed.copies
    bound = Fraction(2 * check.M, check.D)
    multiplicity_ok = multiplicity <= bound
    group_ok = _group_size_ok(placed)
    accepted = not design and not restored
    if check.admitted:
        balance_sentence = "The Lemma 4.5 balance inequalities hold. "
    else:
        balance_sentence = (
            "The Lemma 4.5 balance inequalities fail, so this placement is not "
            "a balancer output. "
        )
    if accepted:
        design_sentence = (
            "design_problems_at and T=0 restoration passed at every t in {1, ..., Q}. "
        )
    else:
        design_sentence = "design_problems_at or T=0 restoration reported a problem. "
    if multiplicity_ok:
        multiplicity_sentence = "Label multiplicity is at most 2M/D. "
    else:
        multiplicity_sentence = (
            f"Label multiplicity {multiplicity} exceeds 2M/D = {bound}. "
        )
    if group_ok:
        group_sentence = "Every group has size at most 2n/r. "
    else:
        group_sentence = "Some group has size above 2n/r. "
    return ScheduleOutcome(
        balance=check,
        placed=placed,
        design_problems=tuple(design),
        restoration_problems=restored,
        accepted=accepted,
        caller_supplied=True,
        built_slice=built,
        allocation_refused=False,
        multiplicity=multiplicity,
        multiplicity_bound=bound,
        multiplicity_ok=multiplicity_ok,
        group_size_ok=group_ok,
        reason=(
            "Schedules were supplied by the caller. They were not chosen by a balancer. "
            + balance_sentence
            + design_sentence
            + group_sentence
            + multiplicity_sentence
        ),
        balancer_attempted=False,
        balancer_built=False,
        design_checked=True,
        design_method=design_method,
    )


def placed_identity_problems(placed: PlacedFormula) -> list[str]:
    """Step 3 checks on one placed slice: matrix, T=0, design, both degree caps.

    The samples are the full-cube sum. This does not speed up the sum, and it
    does not time a refused matrix.
    """
    problems = [f"matrix: {problem}" for problem in auxiliary_problems(placed.auxiliary)]
    problems.extend(f"restore: {problem}" for problem in restoration_problems(placed))
    for t in range(1, placed.auxiliary.Q + 1):
        problems.extend(
            f"design: {problem}" for problem in design_problems_at(placed, t)
        )
    report = degree_cap_report(placed)
    if not report.ok:
        problems.append(
            "degree: "
            f"H(0)={report.direct_at_zero} interpolated={report.h0_occurrences} "
            f"#F={report.brute} deg={report.observed_degree} "
            f"caps={report.cap_occurrences}/{report.cap_step1}"
        )
    return problems


@dataclass(frozen=True)
class Step7Example:
    """One explicit slice above the old cap of 64.

    ``heavy`` marks the side-length-256 point. ``verify --quick`` skips it.
    The schedules are frozen data. They are not the output of a balancer.
    """

    name: str
    notes: str
    placed: PlacedFormula
    sample_points: tuple[Fraction, ...]
    heavy: bool

    @property
    def auxiliary(self) -> AuxiliarySlice:
        return self.placed.auxiliary


def _positive_three_clause() -> CNF:
    return CNF(3, [[(0, False), (1, False), (2, True)]])


def step7_examples(*, heavy: bool = True) -> list[Step7Example]:
    """Slices that the cap of 64 refused.

    ``r3-q2`` has side length 81, R < Q, and degree 0. ``r2-q4`` has side
    length 256, R < Q, and a row of positive degree. Both use one width-3
    clause and explicit schedules. ``heavy=False`` drops the 256-point.
    """
    examples: list[Step7Example] = []
    small = build_auxiliary_slice(3, 2, 1)
    small_placed = place_formula(
        CNF(3, [[(0, False), (1, False), (2, False)]]),
        small,
        variable_schedules=((0, 0), (0, 0), (1, 1)),
        clause_schedules=((0, 0),),
    )
    examples.append(
        Step7Example(
            name="r3-q2",
            notes=(
                "r=3, Q=2, R=1, side length 81, R<Q. One positive 3-clause. "
                "Schedules are explicit data, not a balancer output"
            ),
            placed=small_placed,
            sample_points=(Fraction(0), Fraction(1), Fraction(2)),
            heavy=False,
        )
    )
    if heavy:
        wide = build_auxiliary_slice(2, 4, 3)
        wide_placed = place_formula(
            _positive_three_clause(),
            wide,
            variable_schedules=((0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 1)),
            clause_schedules=((0, 0, 0, 0),),
        )
        examples.append(
            Step7Example(
                name="r2-q4",
                notes=(
                    "r=2, Q=4, R=3, side length 256, R<Q. One 3-clause with a "
                    "negated literal. Schedules are explicit data, not a balancer output"
                ),
                placed=wide_placed,
                sample_points=(
                    Fraction(0),
                    Fraction(1),
                    Fraction(2),
                    Fraction(3),
                    Fraction(4),
                ),
                heavy=True,
            )
        )
    return examples

