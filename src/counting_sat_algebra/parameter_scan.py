"""
Parameter scan for the pieces that are actually in this tree (Step 6).

The grid is fixed. Along n, the popcount block design is compared with brute
force where both fit in the rails already used by Steps 2 and 5. Lagrange,
paired Hermite, the bounded-occurrence degree check, implication-cycle
copying, and the tiny replacement matrix are run on the same kind of budget.
A cell that would pass the brute-force wall, the assignment budget, or the
measurement-tuple cap is not started. Full M0, a crossing clause with no
Step 3 placement, and a Hermite jet that would remove the 2^n pass are
recorded as refusals. Those constructions are not in this tree, so the scan
does not time them.

Two profiles share one predicate:

- ``quick`` is the small grid. Pytest and ``verify --quick`` run it.
- ``full`` adds larger n that still finish on one machine. ``verify --step6``
  runs it. It is not the default pytest path.

Wall times are whatever the clock on that process did. They are not a
complexity result. The figures that travel are the operation counts: the cube
2^n, the popcount tuple count (r+1)^b, the local-evaluation count, and the
planned assignment budget. A refused pipeline has no timing.

This module does not implement Theorem 1.1. It does not refute ETH. Nothing
it runs is a 2^{o(n)} counting algorithm.
"""

from __future__ import annotations

import json
import random
import time
from dataclasses import dataclass, replace
from fractions import Fraction
from pathlib import Path
from typing import Callable

from .blocked_sum import (
    DEFAULT_MAX_BLOCK_SIZE,
    DEFAULT_MAX_MEASUREMENT_TUPLES,
    blocked_sum,
    popcount_group_instance,
    popcount_operation_counts,
)
from .bounded_occurrence import (
    check_bounded_occurrence,
    degree_bound_3n,
    random_bounded_3cnf,
    run_degree_bound_experiment,
)
from .cnf import CNF, brute_count, random_3cnf
from .hermite import eval_H_hasse_jet, hermite_eval, hermite_multiplicities
from .replacement_matrix import (
    FULL_M0_LOG2_DIMENSION_AT_R2,
    FULL_M0_MIN_Q,
    FULL_M0_MIN_R,
    SLICE_MAX_DIMENSION,
    build_auxiliary_slice,
    degree_cap_report,
    place_formula,
)
from .replacements import degree_bound, eval_H_at, interpolate_H0
from .sparsify import (
    BRUTE_N_MAX,
    PIPELINE_ASSIGNMENT_BUDGET,
    copying_lower_bound,
    count_agreement,
    cross_block_formula,
    run_step2_pipeline,
    sparsify,
)

# Re-exported so a caller can name the rails without importing five modules.
ASSIGNMENT_BUDGET = PIPELINE_ASSIGNMENT_BUDGET
MEASUREMENT_TUPLE_CAP = DEFAULT_MAX_MEASUREMENT_TUPLES

TIMING_DISCLAIMER = (
    "Wall times are machine-local. They are not a complexity result and they "
    "are not a 2^{o(n)} bound. The reproducible figures are the operation "
    "counts: the cube 2^n, the popcount tuple count (r+1)^b, the "
    "local-evaluation count, and the planned cube-assignment budget. A refused "
    "pipeline has no timing."
)

DISTANCE_TO_FULL_ALGORITHM: tuple[str, ...] = (
    "Full M0 dense allocation at R=32, Q=192, D=2^384 remains refused. "
    "Step 11 provides structured coordinate queries at those parameters. "
    "Step 12 connects them to finite sparse placements and the common count "
    "entry point: n=16384, Q=3, R=2 uses 59392 polynomial entries and passes "
    "all numerical hypotheses and design checks. Its shared-measurement sum is "
    "refused before enumeration. Step 14 recovers count 1 from the same global "
    "placement through disjoint replacement supports, scanning 96 local Boolean "
    "vectors with memoization and 49152 without it. Execution of the full "
    "parameter family remains open.",
    "The default scheduler now uses exact reference conditional expectations and "
    "rational kernels; the admitted n=1024 placement passes. The legacy Boolean "
    "mode preserves the 2**32 budget refusal. Lemma 4.7's exact tail is now checked: "
    "the n=1024, Q=2, R=1 point fails it and the n=16384, Q=3, R=2 point passes. "
    "Complete asymptotic parameter selection "
    "and execution on general formulas remain open.",
    "Sparsification that keeps n' = n + o(n) is not implemented. Implication-cycle "
    "copying can produce n'=36 from n=5; Section 6's disjoint sparsification and "
    "the sparse-polynomial equation construction remain open.",
    "A subexponential measurement space for general inputs is not implemented. "
    "The Step 2 popcount tuples remain 3^{n/2} at r=2. Step 10 streams the "
    "reference's measurement boxes; a completed parameter and reduction pipeline "
    "is needed for the asymptotic bound.",
    "The original Hermite oracle still walks 2^n and its full jet has quadratic "
    "arithmetic per assignment. Step 10 instead recovers jets through restricted "
    "local sums and fixed substitutions, checked on supplied finite designs. "
    "Step 13 uses it on independently placed CNF components and counts "
    "10923 disjoint three-literal clauses as 7**10923. Step 15 counts one "
    "connected 1024-variable path plus two free bits as 4*Fibonacci(1026) with "
    "elimination width 2, and refuses a width-13 factor before building tables. "
    "Independent-component counts and disjoint-support counts are not a general "
    "connected-instance bound. General large connected components still face "
    "resource limits; the general-input reduction remains open.",
    "A crossing clause still fails the Step 2 block factorization. Step 10 can "
    "count crossing clauses on supplied designs, and Step 13 places small "
    "components afresh. Copying does not "
    "construct those placements for a general reduced formula.",
    "Theorem 1.1 is not implemented. This repository does not refute ETH. "
    "Step 15's finite Taylor tensor elimination counts connected supports only "
    "when intermediate tables fit explicit budgets; it is not Appendix D's "
    "general modular contraction or its asymptotic bound. Full bit-complexity "
    "and space analysis, and the Lean machine proof, remain open.",
)

_PROFILES = ("quick", "full")


@dataclass(frozen=True)
class ScanSpec:
    """One grid cell and the qualitative outcome the suite locks."""

    name: str
    step: str
    kind: str
    expect: str
    note: str
    profiles: frozenset[str]
    times_brute: bool
    times_pipeline: bool
    n: int | None = None
    groups: int | None = None
    block_size: int | None = None
    R: int | None = None
    degree_m: int | None = None
    degree_seed: int | None = None
    replacement_seed: int | None = None

    def __post_init__(self) -> None:
        if self.expect not in {"agreed", "refused"}:
            raise ValueError(f"{self.name}: expect must be agreed or refused")
        if not self.profiles or not set(self.profiles) <= set(_PROFILES):
            raise ValueError(f"{self.name}: bad profiles {self.profiles!r}")
        if self.times_pipeline and self.expect != "agreed":
            raise ValueError(f"{self.name}: a refused pipeline is not timed")


@dataclass(frozen=True)
class ScanCell:
    """What one cell did, or the plan it refused to run."""

    name: str
    step: str
    status: str
    expect: str
    reason: str
    note: str
    measured: bool
    n: int | None = None
    brute: int | None = None
    compared: Fraction | None = None
    cube: int | None = None
    tuples: int | None = None
    local_evaluations: int | None = None
    assignments: int | None = None
    blocks: int | None = None
    block_size: int | None = None
    R: int | None = None
    degree_cap: int | None = None
    observed_degree: int | None = None
    n_after: int | None = None
    occurrence_after: int | None = None
    lower_bound_n: int | None = None
    cube_passes: int | None = None
    jet_width: int | None = None
    side_log2: int | None = None
    seconds: float | None = None
    seconds_brute: float | None = None
    seconds_full: float | None = None
    problems: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return self.status == self.expect and not self.problems


@dataclass(frozen=True)
class ScanReport:
    """A finished profile: cells, the timing disclaimer, and the distance list."""

    profile: str
    cells: tuple[ScanCell, ...]
    disclaimer: str
    distance: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return bool(self.cells) and all(cell.ok for cell in self.cells)


def _positive_clause() -> CNF:
    return CNF(3, [[(0, False), (1, False), (2, False)]])


def _units(count: int) -> CNF:
    return CNF(1, [[(0, False)] for _ in range(count)])


def _blowup_cnf() -> CNF:
    """The Step 5 row whose reduced length is already past the brute wall."""
    return random_3cnf(5, 12, rng=random.Random(0))


def _random_heavy_cnf() -> CNF:
    """A 3-CNF with n = 4 and m = 8. Copying pushes n' past 20."""
    return random_3cnf(4, 8, rng=random.Random(1))


def lagrange_assignments(groups: int, block_size: int, R: int) -> tuple[int, int, int]:
    """Return ``(n, degree cap, planned cube assignments)`` for the Lagrange cell.

    ``interpolate_H0`` evaluates H at t = 1, …, D+1. The cell then evaluates
    H(0) once more and compares it with brute force. The plan is
    ``(D + 2) · 2^n`` with D = b · min(r, 3) · R.
    """
    n = groups * block_size
    degree_cap = groups * min(block_size, 3) * R
    return n, degree_cap, (degree_cap + 2) * (1 << n)


def hermite_paired_plan(groups: int, block_size: int, R: int) -> tuple[int, int, int, int, int]:
    """Return ``(n, D, passes, jet width, assignments)`` for the paired scheme.

    Assignments count one cube pass per distinct t, plus one direct pass at
    t = 0. The single-point jet of order D is not included. It is a different
    cell, and it is not run.
    """
    n, degree_cap, _ignored = lagrange_assignments(groups, block_size, R)
    multiplicities = hermite_multiplicities(degree_cap, 2)
    passes = len(multiplicities)
    jet_width = max(multiplicities)
    assignments = (passes + 1) * (1 << n)
    return n, degree_cap, passes, jet_width, assignments


def degree_plan(n: int, m: int, R: int) -> tuple[int, int, int, int]:
    """Return ``(tight cap, loose cap, window, assignments)`` for Step 1's sampler.

    The experiment samples H(1), …, H(window) with
    ``window = max(m·3·R, 3·n·R) + 2``. That is the same window
    ``run_degree_bound_experiment`` builds. Brute force is a separate 2^n
    pass and is not included here.
    """
    tight = m * 3 * R
    loose = degree_bound_3n(n, R)
    window = max(tight, loose) + 2
    return tight, loose, window, window * (1 << n)


def _blocked_expect(groups: int, block_size: int) -> str:
    n = groups * block_size
    tuples, _local, _cube = popcount_operation_counts(groups, block_size)
    if (
        n > BRUTE_N_MAX
        or tuples > MEASUREMENT_TUPLE_CAP
        or block_size > DEFAULT_MAX_BLOCK_SIZE
    ):
        return "refused"
    return "agreed"


def _lagrange_expect(groups: int, block_size: int, R: int) -> str:
    n, _degree_cap, assignments = lagrange_assignments(groups, block_size, R)
    if n > BRUTE_N_MAX or assignments > ASSIGNMENT_BUDGET:
        return "refused"
    return "agreed"


def _hermite_expect(groups: int, block_size: int, R: int) -> str:
    n, _degree_cap, _passes, _width, assignments = hermite_paired_plan(groups, block_size, R)
    if n > BRUTE_N_MAX or assignments > ASSIGNMENT_BUDGET:
        return "refused"
    return "agreed"


def _degree_expect(n: int, m: int, R: int) -> str:
    _tight, _loose, _window, assignments = degree_plan(n, m, R)
    if n > BRUTE_N_MAX or assignments > ASSIGNMENT_BUDGET:
        return "refused"
    return "agreed"


def _round_time(value: float | None) -> float | None:
    if value is None:
        return None
    return round(value, 6)


def _timed(function: Callable[[], object]) -> tuple[object, float]:
    started = time.perf_counter()
    value = function()
    return value, _round_time(time.perf_counter() - started) or 0.0


def _cell(spec: ScanSpec, *, status: str, reason: str, problems: tuple[str, ...] = (), **kwargs: object) -> ScanCell:
    extra = list(problems)
    if status not in {"agreed", "refused", "disagreed"}:
        extra.append(f"unknown status {status!r}")
        status = "disagreed"
    if status != spec.expect:
        extra.append(f"status {status}, expected {spec.expect}")
    return ScanCell(
        name=spec.name,
        step=spec.step,
        status=status,
        expect=spec.expect,
        reason=reason,
        note=spec.note,
        problems=tuple(extra),
        **kwargs,  # type: ignore[arg-type]
    )


def _timing_problems(spec: ScanSpec, cell: ScanCell) -> tuple[str, ...]:
    problems: list[str] = []
    if spec.times_pipeline:
        if cell.seconds is None:
            problems.append("the pipeline was supposed to be timed")
    elif cell.seconds is not None:
        problems.append("a pipeline timing was recorded for a cell that does not run it")
    if cell.status == "refused" and cell.seconds is not None:
        problems.append("a refused pipeline has a timing")
    if spec.times_brute:
        if cell.seconds_brute is None:
            problems.append("brute force was supposed to be timed")
    elif cell.seconds_brute is not None:
        problems.append("a brute-force timing was recorded although that cube was not counted")
    return tuple(problems)


def _popcount_counts(spec: ScanSpec) -> tuple[int, int, int, int, int]:
    assert spec.groups is not None and spec.block_size is not None and spec.R is not None
    tuples, local, cube = popcount_operation_counts(spec.groups, spec.block_size)
    degree_cap = spec.groups * min(spec.block_size, 3) * spec.R
    n = spec.groups * spec.block_size
    return n, tuples, local, cube, degree_cap


def _run_blocked(spec: ScanSpec) -> ScanCell:
    n, tuples, local, cube, degree_cap = _popcount_counts(spec)
    assert spec.groups is not None and spec.block_size is not None and spec.R is not None
    base = dict(
        measured=False,
        n=n,
        brute=None,
        compared=None,
        cube=cube,
        tuples=tuples,
        local_evaluations=local,
        assignments=cube,
        blocks=spec.groups,
        block_size=spec.block_size,
        R=spec.R,
        degree_cap=degree_cap,
        seconds=None,
        seconds_brute=None,
        seconds_full=None,
    )
    if n > BRUTE_N_MAX:
        cnf, _structured = popcount_group_instance(spec.groups, spec.block_size, spec.R)
        try:
            brute_count(cnf)
        except ValueError as exc:
            brute_message = str(exc)
        else:
            return _cell(
                spec,
                status="disagreed",
                reason=f"brute_count accepted n={n}, above the wall {BRUTE_N_MAX}",
                **base,
            )
        tuple_note = ""
        if tuples > MEASUREMENT_TUPLE_CAP:
            tuple_note = (
                f" The popcount tuple count {tuples} = ({spec.block_size}+1)^{spec.groups} "
                f"is also above the measurement cap {MEASUREMENT_TUPLE_CAP}."
            )
        return _cell(
            spec,
            status="refused",
            reason=(
                f"{brute_message} The blocked sum was not started, and neither "
                f"the cube 2^{n} nor the {tuples} measurement tuples were enumerated."
                f"{tuple_note} No timing was recorded."
            ),
            **base,
        )

    cnf, structured = popcount_group_instance(spec.groups, spec.block_size, spec.R)
    if structured.R != spec.R or cnf.n != n:
        return _cell(
            spec,
            status="disagreed",
            reason="popcount instance drifted from the scan cell",
            **base,
        )

    if tuples > MEASUREMENT_TUPLE_CAP or spec.block_size > DEFAULT_MAX_BLOCK_SIZE:
        brute, seconds_brute = _timed(lambda: brute_count(cnf))
        try:
            blocked_sum(cnf, structured, 0)
        except ValueError as exc:
            message = str(exc)
        else:
            return _cell(
                spec,
                status="disagreed",
                reason="blocked_sum ran past its measurement cap",
                brute=brute if isinstance(brute, int) else None,
                seconds_brute=seconds_brute,
                **{key: value for key, value in base.items() if key not in {"brute", "seconds_brute"}},
            )
        if "above the cap" not in message and "enumeration cap" not in message:
            return _cell(
                spec,
                status="disagreed",
                reason=f"blocked_sum refused for an unexpected reason: {message}",
                brute=brute if isinstance(brute, int) else None,
                seconds_brute=seconds_brute,
                **{key: value for key, value in base.items() if key not in {"brute", "seconds_brute"}},
            )
        return _cell(
            spec,
            status="refused",
            reason=(
                f"Brute force counted the cube 2^{n} = {cube}. The blocked sum "
                f"refused before its measurement loop: {message} "
                f"The closed-form local-evaluation count {local} was not performed. "
                "The blocked sum was not timed."
            ),
            brute=brute if isinstance(brute, int) else None,
            seconds_brute=seconds_brute,
            **{key: value for key, value in base.items() if key not in {"brute", "seconds_brute"}},
        )

    brute, seconds_brute = _timed(lambda: brute_count(cnf))
    blocked, seconds_blocked = _timed(lambda: blocked_sum(cnf, structured, 0))
    full, seconds_full = _timed(lambda: eval_H_at(cnf, structured.family, Fraction(0)))
    problems: list[str] = []
    if blocked.measurement_tuples != tuples:
        problems.append(
            f"tuple count {blocked.measurement_tuples} != closed form {tuples}"
        )
    if blocked.local_evaluations != local:
        problems.append(
            f"local evaluations {blocked.local_evaluations} != closed form {local}"
        )
    if blocked.full_space != cube:
        problems.append(f"full space {blocked.full_space} != 2^{n}")
    if blocked.degree_cap != degree_cap:
        problems.append(f"design D {blocked.degree_cap} != {degree_cap}")
    matched = (
        blocked.value == full == Fraction(brute)
        and Fraction(brute).denominator == 1
    )
    status = "agreed" if matched and not problems else "disagreed"
    reason = (
        f"Blocked H(0) = full H(0) = #F = {brute}. "
        f"Tuples ({spec.block_size}+1)^{spec.groups} = {tuples}, "
        f"local evaluations {local}, cube 2^{n} = {cube}. "
        "D is stored on the design and is not read by the sum."
        if status == "agreed"
        else f"blocked {blocked.value}, full {full}, #F {brute}; " + "; ".join(problems)
    )
    return _cell(
        spec,
        status=status,
        reason=reason,
        problems=tuple(problems),
        measured=True,
        n=n,
        brute=brute if isinstance(brute, int) else None,
        compared=blocked.value,
        cube=cube,
        tuples=blocked.measurement_tuples,
        local_evaluations=blocked.local_evaluations,
        assignments=cube,
        blocks=spec.groups,
        block_size=spec.block_size,
        R=spec.R,
        degree_cap=degree_cap,
        seconds=seconds_blocked,
        seconds_brute=seconds_brute,
        seconds_full=seconds_full,
    )


def _run_lagrange(spec: ScanSpec) -> ScanCell:
    assert spec.groups is not None and spec.block_size is not None and spec.R is not None
    n, degree_cap, assignments = lagrange_assignments(spec.groups, spec.block_size, spec.R)
    width = min(spec.block_size, 3)
    base = dict(
        measured=False,
        n=n,
        brute=None,
        compared=None,
        cube=1 << n,
        assignments=assignments,
        blocks=spec.groups,
        block_size=spec.block_size,
        R=spec.R,
        degree_cap=degree_cap,
        cube_passes=degree_cap + 1,
        seconds=None,
        seconds_brute=None,
        seconds_full=None,
    )
    if spec.expect == "refused":
        return _cell(
            spec,
            status="refused",
            reason=(
                f"Lagrange on this popcount formula would enumerate {assignments} "
                f"cube assignments, (D+2)·2^{n} with D = {degree_cap}. "
                f"The budget is {ASSIGNMENT_BUDGET}. "
                + (
                    f"n={n} also exceeds the brute-force wall {BRUTE_N_MAX}. "
                    if n > BRUTE_N_MAX
                    else ""
                )
                + "H(0) was not computed. No timing was recorded. "
                "Each sample of H is still one pass over the cube."
            ),
            **base,
        )
    cnf, structured = popcount_group_instance(spec.groups, spec.block_size, spec.R)
    if degree_bound(cnf, spec.R, k=width) != degree_cap:
        return _cell(spec, status="disagreed", reason="degree cap drifted", **base)

    def _interpolate() -> tuple[Fraction, Fraction]:
        h0, cap, family = interpolate_H0(cnf, fam=structured.family, k=width)
        if cap != degree_cap or family is not structured.family:
            raise RuntimeError("interpolate_H0 drifted from the scan family")
        direct = eval_H_at(cnf, structured.family, Fraction(0))
        return h0, direct

    brute, seconds_brute = _timed(lambda: brute_count(cnf))
    interpolated, seconds = _timed(_interpolate)
    h0, direct = interpolated
    matched = h0 == direct == Fraction(brute) and h0.denominator == 1
    status = "agreed" if matched else "disagreed"
    reason = (
        f"Lagrange H(0) = direct H(0) = #F = {brute} at D = {degree_cap}. "
        f"The path enumerated {assignments} cube assignments, (D+2)·2^{n}. "
        "That is still one 2^n pass per sample."
        if status == "agreed"
        else f"H(0)={h0}, direct={direct}, #F={brute}"
    )
    return _cell(
        spec,
        status=status,
        reason=reason,
        measured=True,
        n=n,
        brute=brute if isinstance(brute, int) else None,
        compared=h0 if isinstance(h0, Fraction) else None,
        cube=1 << n,
        assignments=assignments,
        blocks=spec.groups,
        block_size=spec.block_size,
        R=spec.R,
        degree_cap=degree_cap,
        cube_passes=degree_cap + 1,
        seconds=seconds,
        seconds_brute=seconds_brute,
        seconds_full=None,
    )


def _run_hermite(spec: ScanSpec) -> ScanCell:
    assert spec.groups is not None and spec.block_size is not None and spec.R is not None
    n, degree_cap, passes, jet_width, assignments = hermite_paired_plan(
        spec.groups, spec.block_size, spec.R
    )
    width = min(spec.block_size, 3)
    base = dict(
        measured=False,
        n=n,
        brute=None,
        compared=None,
        cube=1 << n,
        assignments=assignments,
        blocks=spec.groups,
        block_size=spec.block_size,
        R=spec.R,
        degree_cap=degree_cap,
        cube_passes=passes,
        jet_width=jet_width,
        seconds=None,
        seconds_brute=None,
        seconds_full=None,
    )
    if spec.expect == "refused":
        return _cell(
            spec,
            status="refused",
            reason=(
                f"Paired Hermite would enumerate {assignments} cube assignments "
                f"({passes} points of jet width {jet_width}, plus H(0)), above the "
                f"budget {ASSIGNMENT_BUDGET} or the brute-force wall. "
                "The jets were not started. No timing was recorded. "
                "The condition count would still have been D+1 "
                f"= {degree_cap + 1}. Each distinct t is one pass over 2^{n}."
            ),
            **base,
        )
    if jet_width < 2 or passes >= degree_cap + 1:
        return _cell(
            spec,
            status="disagreed",
            reason="paired scheme did not pack D+1 conditions into fewer points",
            **base,
        )
    cnf, structured = popcount_group_instance(spec.groups, spec.block_size, spec.R)
    multiplicities = hermite_multiplicities(degree_cap, 2)

    def _recover() -> Fraction:
        samples: list[tuple[Fraction, tuple[Fraction, ...]]] = []
        for index, multiplicity in enumerate(multiplicities):
            jet = eval_H_hasse_jet(
                cnf,
                structured.family,
                Fraction(index + 1),
                multiplicity - 1,
            )
            if len(jet) != multiplicity:
                raise RuntimeError("Hasse jet length drifted")
            samples.append((Fraction(index + 1), jet))
        return hermite_eval(samples, Fraction(0))

    brute, seconds_brute = _timed(lambda: brute_count(cnf))
    h0, seconds = _timed(_recover)
    direct, seconds_full = _timed(lambda: eval_H_at(cnf, structured.family, Fraction(0)))
    conditions = sum(multiplicities)
    problems: list[str] = []
    if conditions != degree_cap + 1:
        problems.append(f"conditions {conditions} != D+1")
    if degree_bound(cnf, spec.R, k=width) != degree_cap:
        problems.append("degree cap drifted")
    matched = h0 == direct == Fraction(brute) and isinstance(h0, Fraction) and h0.denominator == 1
    status = "agreed" if matched and not problems else "disagreed"
    reason = (
        f"Paired Hermite H(0) = direct H(0) = #F = {brute}. "
        f"{passes} points, jet width {jet_width}, conditions {conditions} = D+1. "
        f"Cube assignments {assignments} = (passes + 1)·2^{n}. "
        "The single-point jet of order D was not run."
        if status == "agreed"
        else f"H(0)={h0}, direct={direct}, #F={brute}; " + "; ".join(problems)
    )
    return _cell(
        spec,
        status=status,
        reason=reason,
        problems=tuple(problems),
        measured=True,
        n=n,
        brute=brute if isinstance(brute, int) else None,
        compared=h0 if isinstance(h0, Fraction) else None,
        cube=1 << n,
        assignments=assignments,
        blocks=spec.groups,
        block_size=spec.block_size,
        R=spec.R,
        degree_cap=degree_cap,
        cube_passes=passes,
        jet_width=jet_width,
        seconds=seconds,
        seconds_brute=seconds_brute,
        seconds_full=seconds_full,
    )


def _run_hermite_cost(spec: ScanSpec) -> ScanCell:
    """Record the single-point jet as not a 2^n reduction, and do not run it.

    The reference formula is the popcount row b = 1, r = 2, where the paired
    cell already compares H(0) with brute force. One node of order D would
    still walk the cube once, with jet width D+1.
    """
    n, degree_cap, _passes, _width, _paired = hermite_paired_plan(1, 2, 1)
    jet_width = degree_cap + 1
    assignments = 1 << n
    return _cell(
        spec,
        status="refused",
        reason=(
            f"A single Hermite node of order D = {degree_cap} on the popcount "
            f"row n = {n} would enumerate the cube once ({assignments} assignments) "
            f"and do quadratic arithmetic in the jet width D+1 = {jet_width}. "
            "That pass was not started. Packing every coefficient into t = 1 "
            "does not evaluate H below exponential cost in n. No timing was recorded."
        ),
        measured=False,
        n=n,
        cube=assignments,
        assignments=assignments,
        blocks=1,
        block_size=2,
        R=1,
        degree_cap=degree_cap,
        cube_passes=1,
        jet_width=jet_width,
        seconds=None,
        seconds_brute=None,
        seconds_full=None,
    )


def _run_degree(spec: ScanSpec) -> ScanCell:
    assert spec.n is not None and spec.degree_m is not None and spec.R is not None
    n, m, R = spec.n, spec.degree_m, spec.R
    tight, loose, window, assignments = degree_plan(n, m, R)
    base = dict(
        measured=False,
        n=n,
        brute=None,
        compared=None,
        cube=1 << n,
        assignments=assignments,
        R=R,
        degree_cap=tight,
        seconds=None,
        seconds_brute=None,
        seconds_full=None,
    )
    if spec.kind == "explicit-clause":
        cnf = _positive_clause()
    else:
        cnf = random_bounded_3cnf(n, m, rng=random.Random(spec.degree_seed))
    report = check_bounded_occurrence(cnf)
    if not report.ok:
        return _cell(
            spec,
            status="disagreed",
            reason="bounded-occurrence checker rejected the scan formula: "
            + "; ".join(report.violations),
            **base,
        )
    if spec.expect == "refused":
        return _cell(
            spec,
            status="refused",
            reason=(
                f"The degree experiment would enumerate {assignments} cube "
                f"assignments, window {window} = max(m·3·R, 3·n·R)+2 with "
                f"caps {tight} and {loose}, above the budget {ASSIGNMENT_BUDGET}. "
                "The formula was built and checked for occurrence ≤ 3. "
                "H was not sampled. Brute force was not run. No timing was recorded."
            ),
            **base,
        )
    seed = 0 if spec.replacement_seed is None else spec.replacement_seed
    brute, seconds_brute = _timed(lambda: brute_count(cnf))
    result, seconds = _timed(
        lambda: run_degree_bound_experiment(cnf, R=R, rng=random.Random(seed))
    )
    ok = bool(result["ok_mkR"] and result["ok_3nR"] and result["brute"] == brute)
    status = "agreed" if ok else "disagreed"
    compared = result["H0_mkR"]
    reason = (
        f"Both caps recover H(0) = #F = {brute}. "
        f"Tight cap m·3·R = {tight}, loose cap 3·n·R = {loose}, "
        f"observed degree {result['observed_degree']}. "
        f"The sampler enumerated {assignments} cube assignments. "
        "Each of those samples is a full pass over the cube."
        if status == "agreed"
        else f"degree experiment disagreed: {result}"
    )
    return _cell(
        spec,
        status=status,
        reason=reason,
        measured=True,
        n=n,
        brute=brute if isinstance(brute, int) else None,
        compared=compared if isinstance(compared, Fraction) else None,
        cube=1 << n,
        assignments=assignments,
        R=R,
        degree_cap=tight,
        observed_degree=result["observed_degree"] if status == "agreed" else None,
        seconds=seconds,
        seconds_brute=seconds_brute,
        seconds_full=None,
    )


def _copy_source(spec: ScanSpec) -> CNF:
    if spec.kind == "copy-popcount":
        cnf, _structured = popcount_group_instance(1, 2, 1)
        return cnf
    if spec.kind == "copy-blowup":
        return _blowup_cnf()
    if spec.kind == "copy-units":
        return _units(21)
    if spec.kind == "copy-random":
        return _random_heavy_cnf()
    if spec.kind == "copy-cross":
        return cross_block_formula()
    raise ValueError(f"unknown copy kind {spec.kind}")


def _run_copy(spec: ScanSpec) -> ScanCell:
    source = _copy_source(spec)
    reduced_spec = sparsify(source)
    bound = copying_lower_bound(reduced_spec.literal_occurrences_before, source.n)
    base = dict(
        measured=False,
        n=source.n,
        n_after=reduced_spec.n_after,
        occurrence_after=reduced_spec.max_occurrence_after,
        lower_bound_n=bound,
        cube=1 << source.n,
        assignments=1 << reduced_spec.n_after,
        seconds=None,
        seconds_brute=None,
        seconds_full=None,
    )
    if reduced_spec.max_occurrence_after > reduced_spec.cap:
        return _cell(
            spec,
            status="disagreed",
            reason=f"occurrence {reduced_spec.max_occurrence_after} exceeds the cap",
            **base,
        )

    if spec.kind == "copy-cross":
        blocks = ((0, 1), (2,))
        brute_source, seconds_source = _timed(lambda: brute_count(reduced_spec.source))
        brute_reduced, seconds_reduced = _timed(lambda: brute_count(reduced_spec.reduced))
        agreement = count_agreement(reduced_spec)
        blocked = run_step2_pipeline(reduced_spec, blocks, R=1)
        problems: list[str] = []
        if agreement.status != "agreed" or brute_source != brute_reduced:
            problems.append(
                f"count status {agreement.status}, source {brute_source}, reduced {brute_reduced}"
            )
        if blocked.status != "refused" or "block" not in blocked.reason:
            problems.append(f"Step 2 status {blocked.status}: {blocked.reason}")
        if blocked.seconds != 0.0 or blocked.blocked_at_zero is not None:
            problems.append("Step 2 timed or evaluated a crossing clause")
        status = "refused" if not problems else "disagreed"
        return _cell(
            spec,
            status=status,
            reason=(
                f"#F' = #F = {brute_reduced} was counted (n {source.n}→{reduced_spec.n_after}). "
                f"The blocked sum was not run: {blocked.reason} "
                "Step 3 is not applied. Copying does not produce a schedule, and "
                "full M0 is not built. The blocked sum was not timed."
            ),
            problems=tuple(problems),
            brute=brute_reduced if isinstance(brute_reduced, int) else None,
            seconds_brute=_round_time((seconds_source or 0.0) + (seconds_reduced or 0.0)),
            **{key: value for key, value in base.items() if key != "seconds_brute"},
        )

    if reduced_spec.n_after > BRUTE_N_MAX:
        brute_source, seconds_brute = _timed(lambda: brute_count(source))
        agreement = count_agreement(reduced_spec)
        problems = []
        if agreement.status != "refused" or agreement.brute_reduced is not None:
            problems.append(f"expected a reduced-cube refusal, got {agreement.status}")
        if agreement.brute_source != brute_source:
            problems.append("source count drifted from count_agreement")
        try:
            brute_count(reduced_spec.reduced)
        except ValueError as exc:
            reduced_message = str(exc)
        else:
            problems.append("brute_count accepted the reduced formula")
            reduced_message = ""
        status = "refused" if not problems else "disagreed"
        return _cell(
            spec,
            status=status,
            reason=(
                f"Source n={source.n} was counted, #F = {brute_source}. "
                f"Reduced n'={reduced_spec.n_after} (lower bound {bound}) was not: "
                f"{reduced_message} Occurrence after copying is "
                f"{reduced_spec.max_occurrence_after}. "
                f"The recovery #F' = #F was not checked. The reduced cube was not timed."
            ),
            problems=tuple(problems),
            brute=brute_source if isinstance(brute_source, int) else None,
            seconds_brute=seconds_brute,
            **{key: value for key, value in base.items() if key != "seconds_brute"},
        )

    brute_source, seconds_source = _timed(lambda: brute_count(source))
    brute_reduced, seconds_reduced = _timed(lambda: brute_count(reduced_spec.reduced))
    agreement = count_agreement(reduced_spec)
    matched = (
        agreement.status == "agreed"
        and brute_source == brute_reduced == agreement.brute_reduced
        and reduced_spec.n_after == source.n
        and reduced_spec.cycle_clause_count == 0
    )
    status = "agreed" if matched else "disagreed"
    counted = (
        (1 << source.n) + (1 << reduced_spec.n_after)
        if isinstance(brute_source, int) and isinstance(brute_reduced, int)
        else base["assignments"]
    )
    return _cell(
        spec,
        status=status,
        reason=(
            f"#F' = #F = {brute_reduced}. The formula already had occurrence ≤ 3, "
            f"so n' = n = {source.n} and the cycle is empty."
            if status == "agreed"
            else f"copy comparison failed: source {brute_source}, reduced {brute_reduced}, {agreement.reason}"
        ),
        measured=status == "agreed",
        n=base["n"],
        brute=brute_reduced if isinstance(brute_reduced, int) else None,
        compared=Fraction(brute_reduced) if isinstance(brute_reduced, int) and status == "agreed" else None,
        cube=base["cube"],
        assignments=counted,
        n_after=base["n_after"],
        occurrence_after=base["occurrence_after"],
        lower_bound_n=base["lower_bound_n"],
        seconds=seconds_reduced,
        seconds_brute=seconds_source,
        seconds_full=None,
    )


def _run_matrix(spec: ScanSpec) -> ScanCell:
    if spec.kind == "matrix-full":
        try:
            build_auxiliary_slice(2, FULL_M0_MIN_Q, FULL_M0_MIN_R)
        except ValueError as exc:
            message = str(exc)
        else:
            return _cell(
                spec,
                status="disagreed",
                reason="build_auxiliary_slice built full M0",
                measured=False,
                side_log2=FULL_M0_LOG2_DIMENSION_AT_R2,
                R=FULL_M0_MIN_R,
                seconds=None,
                seconds_brute=None,
                seconds_full=None,
            )
        refused_ok = (
            "not built" in message
            and f"2**{FULL_M0_LOG2_DIMENSION_AT_R2}" in message
        )
        status = "refused" if refused_ok else "disagreed"
        return _cell(
            spec,
            status=status,
            reason=(
                f"build_auxiliary_slice(r=2, Q={FULL_M0_MIN_Q}, R={FULL_M0_MIN_R}) "
                f"refused before allocating the matrix. The side length is "
                f"2^{FULL_M0_LOG2_DIMENSION_AT_R2}, above the slice cap "
                f"{SLICE_MAX_DIMENSION}. No scheduler was invoked. "
                "No timing was recorded."
            ),
            problems=() if refused_ok else (message,),
            measured=False,
            side_log2=FULL_M0_LOG2_DIMENSION_AT_R2,
            R=FULL_M0_MIN_R,
            degree_cap=None,
            seconds=None,
            seconds_brute=None,
            seconds_full=None,
        )

    auxiliary = build_auxiliary_slice(2, 1, 1)
    placed = place_formula(
        CNF(2, [[(0, False)], [(1, False)]]),
        auxiliary,
        ((0,), (1,)),
        ((1,), (0,)),
    )
    brute, seconds_brute = _timed(lambda: brute_count(placed.cnf))
    report, seconds = _timed(lambda: degree_cap_report(placed))
    status = "agreed" if report.ok and report.brute == brute else "disagreed"
    window = max(report.cap_occurrences, report.cap_step1, 0) + 2
    assignments = (window + 1) * (1 << placed.cnf.n)
    return _cell(
        spec,
        status=status,
        reason=(
            f"Tiny M0 at side length {auxiliary.dimension} restores H(0) = #F = {brute}. "
            f"Observed degree {report.observed_degree}, caps "
            f"{report.cap_occurrences} and {report.cap_step1}. "
            f"This is the hand-written swap-4 placement, not full M0. "
            f"The degree report enumerated {assignments} cube assignments on n = {placed.cnf.n}."
            if status == "agreed"
            else f"matrix slice disagreed: H(0)={report.h0_occurrences}, #F={report.brute}"
        ),
        measured=status == "agreed",
        n=placed.cnf.n,
        brute=brute if isinstance(brute, int) else None,
        compared=report.direct_at_zero,
        cube=1 << placed.cnf.n,
        assignments=assignments,
        R=auxiliary.R,
        degree_cap=report.cap_step1,
        observed_degree=report.observed_degree,
        side_log2=None,
        seconds=seconds,
        seconds_brute=seconds_brute,
        seconds_full=None,
    )


_RUNNERS: dict[str, Callable[[ScanSpec], ScanCell]] = {
    "blocked": _run_blocked,
    "lagrange": _run_lagrange,
    "hermite": _run_hermite,
    "hermite-cost": _run_hermite_cost,
    "degree": _run_degree,
    "copy": _run_copy,
    "matrix": _run_matrix,
}


def _spec(
    name: str,
    step: str,
    kind: str,
    note: str,
    profiles: frozenset[str],
    *,
    expect: str,
    times_brute: bool,
    times_pipeline: bool,
    **kwargs: object,
) -> ScanSpec:
    return ScanSpec(
        name=name,
        step=step,
        kind=kind,
        expect=expect,
        note=note,
        profiles=profiles,
        times_brute=times_brute,
        times_pipeline=times_pipeline,
        **kwargs,  # type: ignore[arg-type]
    )


_QUICK = frozenset({"quick", "full"})
_FULL = frozenset({"full"})


def _blocked_spec(
    name: str,
    groups: int,
    block_size: int,
    R: int,
    profiles: frozenset[str],
    note: str,
) -> ScanSpec:
    expect = _blocked_expect(groups, block_size)
    n = groups * block_size
    times_brute = n <= BRUTE_N_MAX
    return _spec(
        name,
        "blocked",
        "blocked",
        note,
        profiles,
        expect=expect,
        times_brute=times_brute,
        times_pipeline=expect == "agreed",
        n=n,
        groups=groups,
        block_size=block_size,
        R=R,
    )


def _lagrange_spec(
    name: str,
    groups: int,
    block_size: int,
    R: int,
    profiles: frozenset[str],
) -> ScanSpec:
    expect = _lagrange_expect(groups, block_size, R)
    n, degree_cap, assignments = lagrange_assignments(groups, block_size, R)
    return _spec(
        name,
        "lagrange",
        "lagrange",
        (
            f"Lagrange on the popcount formula, D = {degree_cap}, "
            f"planned assignments {assignments}."
        ),
        profiles,
        expect=expect,
        times_brute=expect == "agreed",
        times_pipeline=expect == "agreed",
        n=n,
        groups=groups,
        block_size=block_size,
        R=R,
    )


def _hermite_spec(
    name: str,
    groups: int,
    block_size: int,
    R: int,
    profiles: frozenset[str],
) -> ScanSpec:
    expect = _hermite_expect(groups, block_size, R)
    n, degree_cap, passes, _width, assignments = hermite_paired_plan(groups, block_size, R)
    return _spec(
        name,
        "hermite",
        "hermite",
        (
            f"Paired Hermite on the popcount formula, D = {degree_cap}, "
            f"{passes} points, planned assignments {assignments}."
        ),
        profiles,
        expect=expect,
        times_brute=expect == "agreed",
        times_pipeline=expect == "agreed",
        n=n,
        groups=groups,
        block_size=block_size,
        R=R,
    )


def _degree_spec(
    name: str,
    n: int,
    m: int,
    R: int,
    kind: str,
    profiles: frozenset[str],
    note: str,
    *,
    degree_seed: int | None = None,
    replacement_seed: int | None = None,
) -> ScanSpec:
    expect = _degree_expect(n, m, R)
    return _spec(
        name,
        "degree",
        kind,
        note,
        profiles,
        expect=expect,
        times_brute=expect == "agreed",
        times_pipeline=expect == "agreed",
        n=n,
        degree_m=m,
        R=R,
        degree_seed=degree_seed,
        replacement_seed=replacement_seed,
    )


def _build_specs() -> tuple[ScanSpec, ...]:
    blowup = sparsify(_blowup_cnf())
    units = sparsify(_units(21))
    heavy = sparsify(_random_heavy_cnf())
    specs: list[ScanSpec] = [
        _blocked_spec(
            "blocked-r2-b1",
            1,
            2,
            1,
            _QUICK,
            "One width-2 clause on one block. Tuple count 3, cube 4, #F = 3.",
        ),
        _blocked_spec(
            "blocked-r2-b2",
            2,
            2,
            1,
            _QUICK,
            "Two blocks of size 2. Local evaluations 24 already exceed the cube 16.",
        ),
        _blocked_spec(
            "blocked-r2-b3",
            3,
            2,
            1,
            _QUICK,
            "Three blocks of size 2. Tuples 27, local evaluations 108, cube 64.",
        ),
        _blocked_spec(
            "blocked-r2-b2-R2",
            2,
            2,
            2,
            _QUICK,
            "Same blocks as blocked-r2-b2 at replacement degree R = 2. "
            "The tuple count does not depend on R. H(0) is still #F.",
        ),
        _blocked_spec(
            "blocked-r1-b4",
            4,
            1,
            1,
            _QUICK,
            "Blocks of size 1. The tuple count equals the cube, 2^n.",
        ),
        _blocked_spec(
            "blocked-r1-b17-tuple-cap",
            17,
            1,
            1,
            _QUICK,
            "n = 17 is inside the brute-force wall. The popcount tuple count "
            "2^17 exceeds the measurement cap, so the blocked sum refuses and "
            "brute force still runs.",
        ),
        _blocked_spec(
            "blocked-r2-b11-brute-wall",
            11,
            2,
            1,
            _QUICK,
            "n = 22 is past the brute-force wall, and the tuple count 3^11 is "
            "past the measurement cap. Neither enumeration is started.",
        ),
        _lagrange_spec("lagrange-r2-b1", 1, 2, 1, _QUICK),
        _lagrange_spec("lagrange-r2-b2", 2, 2, 1, _QUICK),
        _lagrange_spec("lagrange-r2-b7-budget", 7, 2, 1, _QUICK),
        _hermite_spec("hermite-paired-r2-b1", 1, 2, 1, _QUICK),
        _hermite_spec("hermite-paired-r2-b7-budget", 7, 2, 1, _QUICK),
        _spec(
            "hermite-single-point-not-run",
            "hermite",
            "hermite-cost",
            "The single-point jet does not remove the 2^n pass. It is not run.",
            _QUICK,
            expect="refused",
            times_brute=False,
            times_pipeline=False,
            n=2,
            groups=1,
            block_size=2,
            R=1,
        ),
        _degree_spec(
            "degree-clause-n3",
            3,
            1,
            1,
            "explicit-clause",
            _QUICK,
            "One positive 3-clause. Both Step 1 caps are checked against brute force.",
            replacement_seed=0,
        ),
        _degree_spec(
            "degree-bounded-n6-m4",
            6,
            4,
            1,
            "bounded-random",
            _QUICK,
            "Step 1 bench shape n = 6, m = 4, formula seed 42, replacement seed 42.",
            degree_seed=42,
            replacement_seed=42,
        ),
        _degree_spec(
            "degree-bounded-n12-m12-budget",
            12,
            12,
            1,
            "bounded-random",
            _QUICK,
            "A legal occurrence-≤3 3-CNF whose degree window exceeds the assignment budget.",
            degree_seed=0,
            replacement_seed=0,
        ),
        _spec(
            "copy-popcount-b1",
            "copy",
            "copy-popcount",
            "The popcount row already has occurrence 1. Copying is the identity.",
            _QUICK,
            expect="agreed",
            times_brute=True,
            times_pipeline=True,
            n=2,
        ),
        _spec(
            "copy-blowup-n5-m12",
            "copy",
            "copy-blowup",
            "Seeded 3-CNF, n = 5, m = 12. The cap-3 bound already exceeds n = 20.",
            _QUICK,
            expect="refused",
            times_brute=True,
            times_pipeline=False,
            n=5,
        ),
        _spec(
            "copy-units-n21",
            "copy",
            "copy-units",
            "Twenty-one positive units. n' = 21, just past the brute-force wall.",
            _QUICK,
            expect="refused",
            times_brute=True,
            times_pipeline=False,
            n=1,
        ),
        _spec(
            "copy-cross-block",
            "copy",
            "copy-cross",
            "A 3-clause meets two blocks. Copying does not repair it, and Step 3 is not applied.",
            _QUICK,
            expect="refused",
            times_brute=True,
            times_pipeline=False,
            n=3,
        ),
        _spec(
            "matrix-swap-4",
            "matrix",
            "matrix-swap",
            "Hand-written 4×4 exchange. H(0) matches brute force. This is not full M0.",
            _QUICK,
            expect="agreed",
            times_brute=True,
            times_pipeline=True,
            n=2,
            R=1,
        ),
        _spec(
            "matrix-full-m0",
            "matrix",
            "matrix-full",
            "R = 32, Q = 192, side length 2^384. Refused before allocation.",
            _QUICK,
            expect="refused",
            times_brute=False,
            times_pipeline=False,
            R=FULL_M0_MIN_R,
        ),
        _blocked_spec(
            "blocked-r2-b4",
            4,
            2,
            1,
            _FULL,
            "Four blocks of size 2. n = 8. Local evaluations exceed the cube.",
        ),
        _blocked_spec(
            "blocked-r2-b5",
            5,
            2,
            1,
            _FULL,
            "Five blocks of size 2. n = 10.",
        ),
        _blocked_spec(
            "blocked-r2-b6",
            6,
            2,
            1,
            _FULL,
            "Six blocks of size 2. n = 12.",
        ),
        _blocked_spec(
            "blocked-r2-b7",
            7,
            2,
            1,
            _FULL,
            "Seven blocks of size 2. n = 14. The blocked sum still runs; "
            "Lagrange and paired Hermite on this same n are over the assignment budget.",
        ),
        _blocked_spec(
            "blocked-r3-b3",
            3,
            3,
            1,
            _FULL,
            "Three blocks of size 3. Tuples 64, local evaluations 384, cube 512.",
        ),
        _lagrange_spec("lagrange-r2-b3", 3, 2, 1, _FULL),
        _lagrange_spec("lagrange-r2-b4", 4, 2, 1, _FULL),
        _lagrange_spec("lagrange-r2-b5", 5, 2, 1, _FULL),
        _hermite_spec("hermite-paired-r2-b2", 2, 2, 1, _FULL),
        _hermite_spec("hermite-paired-r2-b3", 3, 2, 1, _FULL),
        _spec(
            "copy-random-n4-m8",
            "copy",
            "copy-random",
            "Seeded 3-CNF, n = 4, m = 8. Copying pushes n' past 20.",
            _FULL,
            expect="refused",
            times_brute=True,
            times_pipeline=False,
            n=4,
        ),
    ]
    # Touch the sparsifications so a drift in n' fails at import, not only at run time.
    if blowup.n_after != 36 or units.n_after != 21 or heavy.n_after != 24:
        raise RuntimeError(
            "copying lengths drifted: "
            f"blowup {blowup.n_after}, units {units.n_after}, heavy {heavy.n_after}"
        )
    names = [spec.name for spec in specs]
    if len(names) != len(set(names)):
        raise RuntimeError("duplicate scan cell names")
    return tuple(specs)


_ALL_SPECS = _build_specs()


def scan_specs(profile: str = "quick") -> tuple[ScanSpec, ...]:
    """Cells for ``quick`` or ``full``. ``full`` includes every quick cell."""
    if profile not in _PROFILES:
        raise ValueError(f"profile must be one of {_PROFILES}, got {profile!r}")
    return tuple(spec for spec in _ALL_SPECS if profile in spec.profiles)


def _runner_for(spec: ScanSpec) -> Callable[[ScanSpec], ScanCell]:
    if spec.kind == "hermite-cost":
        return _run_hermite_cost
    return _RUNNERS[spec.step]


def run_parameter_scan(profile: str = "quick") -> ScanReport:
    """Run one profile. A cell whose status is not the one in its spec fails the report."""
    specs = scan_specs(profile)
    cells: list[ScanCell] = []
    for spec in specs:
        cell = _runner_for(spec)(spec)
        extra = _timing_problems(spec, cell)
        if extra:
            cell = replace(cell, problems=cell.problems + extra)
        cells.append(cell)
    return ScanReport(
        profile=profile,
        cells=tuple(cells),
        disclaimer=TIMING_DISCLAIMER,
        distance=DISTANCE_TO_FULL_ALGORITHM,
    )


def _json_value(value: object) -> object:
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    return value


def cell_record(cell: ScanCell) -> dict:
    """JSON-ready record. Fractions become strings. Refused pipelines keep null timings."""
    raw = {
        "name": cell.name,
        "step": cell.step,
        "status": cell.status,
        "expect": cell.expect,
        "measured": cell.measured,
        "reason": cell.reason,
        "note": cell.note,
        "n": cell.n,
        "brute": cell.brute,
        "compared": cell.compared,
        "cube": cell.cube,
        "tuples": cell.tuples,
        "local_evaluations": cell.local_evaluations,
        "assignments": cell.assignments,
        "blocks": cell.blocks,
        "block_size": cell.block_size,
        "R": cell.R,
        "degree_cap": cell.degree_cap,
        "observed_degree": cell.observed_degree,
        "n_after": cell.n_after,
        "occurrence_after": cell.occurrence_after,
        "lower_bound_n": cell.lower_bound_n,
        "cube_passes": cell.cube_passes,
        "jet_width": cell.jet_width,
        "side_log2": cell.side_log2,
        "seconds": cell.seconds,
        "seconds_brute": cell.seconds_brute,
        "seconds_full": cell.seconds_full,
        "problems": list(cell.problems),
    }
    return {key: _json_value(value) for key, value in raw.items()}


def scan_document(report: ScanReport) -> dict:
    """Document for ``benches/parameter_scan.json``. Timings are machine-local."""
    return {
        "profile": report.profile,
        "disclaimer": report.disclaimer,
        "distance": list(report.distance),
        "description": (
            "Step 6 parameter scan. Statuses and operation counts are the "
            "reproducible record. seconds, seconds_brute, and seconds_full are "
            "one machine's wall times and are not a complexity result. "
            "Regenerate with counting_sat_algebra.parameter_scan.write_scan_document."
        ),
        "cells": [cell_record(cell) for cell in report.cells],
    }


def scan_bench_path(root: Path | None = None) -> Path:
    base = root or Path(__file__).resolve().parents[2]
    return base / "benches" / "parameter_scan.json"


def write_scan_document(path: Path | None = None, profile: str = "full") -> dict:
    """Run ``profile`` and write the bench document. Default profile is ``full``."""
    document = scan_document(run_parameter_scan(profile))
    target = path or scan_bench_path()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return document


def load_scan_document(path: Path | None = None) -> dict:
    target = path or scan_bench_path()
    return json.loads(target.read_text(encoding="utf-8"))


def _format_seconds(value: float) -> str:
    """Render a finished timing. Sub-millisecond runs stay distinct from ``not timed``."""
    if value < 0.0001:
        return "<0.0001s"
    return f"{value:.4f}s"


def format_scan_line(index: int, total: int, cell: ScanCell) -> str:
    """One CLI row. Refused pipelines print ``not timed`` instead of a zero."""
    if cell.brute is None:
        brute_text = "—"
    else:
        brute_text = str(cell.brute)
    if cell.compared is None:
        compared_text = "—"
    else:
        compared_text = str(cell.compared)
    n_text = "—" if cell.n is None else str(cell.n)
    enumeration: list[str] = []
    if cell.tuples is not None:
        enumeration.append(f"tuples={cell.tuples}")
    if cell.local_evaluations is not None:
        enumeration.append(f"local={cell.local_evaluations}")
    if cell.cube is not None:
        enumeration.append(f"cube={cell.cube}")
    if cell.assignments is not None:
        enumeration.append(f"assignments={cell.assignments}")
    facts: list[str] = []
    if cell.n_after is not None:
        facts.append(f"n'={cell.n_after}")
    if cell.side_log2 is not None:
        facts.append(f"side=2^{cell.side_log2}")
    enum_text = " ".join(enumeration)
    if enum_text and not cell.measured:
        enum_text = "plan " + enum_text
    count_text = " ".join(part for part in (enum_text, " ".join(facts)) if part)
    if cell.seconds is None and cell.seconds_brute is None:
        clock = "not timed"
    else:
        parts = []
        if cell.seconds_brute is not None:
            parts.append(f"brute {_format_seconds(cell.seconds_brute)}")
        if cell.seconds is not None:
            parts.append(f"pipeline {_format_seconds(cell.seconds)}")
        if cell.seconds_full is not None:
            parts.append(f"full {_format_seconds(cell.seconds_full)}")
        clock = ", ".join(parts)
    label = {"agreed": "AGREED", "refused": "REFUSED", "disagreed": "FAIL"}.get(
        cell.status, cell.status.upper()
    )
    return (
        f"  [H {index:02d}/{total}] {cell.name}  step={cell.step}  n={n_text}  "
        f"#F={brute_text}  compared={compared_text}  {count_text}  {label}  ({clock})"
    )


def format_distance(distance: tuple[str, ...] | None = None) -> str:
    lines = ["Distance to a full counting algorithm:"]
    for item in distance if distance is not None else DISTANCE_TO_FULL_ALGORITHM:
        lines.append(f"  - {item}")
    return "\n".join(lines)
