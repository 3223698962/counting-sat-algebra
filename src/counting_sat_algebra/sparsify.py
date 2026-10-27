"""
Counting-preserving variable copying for a small CNF (Step 5).

The aim is an occurrence bound, not a shorter clause list. A variable that
already occurs in at most B clauses is kept as one variable. A variable that
occurs d > B times is replaced by

    c = ceil(d / (B - 2))

copies. Occurrences are dealt out in clause order, at most B − 2 per copy,
and the copies are tied together by the implication cycle

    (¬y_1 ∨ y_2) ∧ (¬y_2 ∨ y_3) ∧ … ∧ (¬y_c ∨ y_1).

Each copy sits in two cycle clauses and at most B − 2 original clauses, so it
occurs in at most B clauses. The smallest cap this gadget can meet is B = 3:
two cycle clauses leave one slot, which is exactly one original occurrence.
That is also the cap in the Step 1 testbed, and it is the default. A smaller
cap is refused.

The cycle is satisfied exactly when every copy carries the same bit. Writing
that bit equal to the original variable is a bijection between source
assignments and reduced assignments that satisfy the cycles, and it preserves
the truth of every original clause. Model counts are therefore equal:

    #F' = #F.

There is no extra factor of 2. A copy is not a free variable.

The output is a mixed-width CNF. Original clauses keep their width, and every
cycle clause has width 2. ``check_bounded_occurrence`` still requires exact
width 3, so it rejects a formula that actually used the cycle. The occurrence
cap is the same number 3. The degree pipeline calls ``interpolate_H0`` with k
equal to the maximum clause width and, for B = 3, also checks the loose cap
from ``degree_bound_3n``. Both bounds apply because each replacement has
degree at most R and, at B = 3, the number of literal occurrences is at most
3n'.

Step 2's blocked sum is run only when the caller supplies a partition of the
source variables that already contains every source clause. Copies are placed
in the block of their source variable, so cycle clauses stay inside that
block. A source clause that meets two blocks still meets two blocks after
copying; that call is recorded as a refusal. Step 3 is not used. Its schedules
are hand-written for specific placements, and this reduction does not produce
one.

Cost, for B = 3. A variable of occurrence d > 3 becomes d variables. With M
the number of literal occurrences, the reduced length satisfies

    n' ≥ n    and    n' ≥ M − 2n.

A 3-CNF with 3m − 2n > 20 therefore leaves the brute-force window n ≤ 20 no
matter how the clauses are arranged. Those instances are recorded and not
evaluated. Interpolation is also refused when the number of cube assignments
it would enumerate exceeds a fixed budget. The identity on the instances that
were evaluated is exact. It is not a subexponential algorithm. Theorem 1.1 is
not implemented.
"""

from __future__ import annotations

import itertools
import random
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .blocked_sum import compare_blocked_to_full, constant_alpha, hand_built_replacement
from .bounded_occurrence import degree_bound_3n, equispaced_degree, occurrence_counts
from .cnf import CNF, brute_count, is_satisfied, random_3cnf
from .design import popcount_design, validate_design
from .replacements import (
    degree_bound,
    eval_H_at,
    interpolate_H0,
    lagrange_eval,
    lagrange_eval_at_zero,
    random_replacements,
)

# Two cycle clauses per copy. One remaining slot is the smallest positive
# number of original occurrences a copy can absorb, so B = 3 is the minimum.
MIN_OCCURRENCE_CAP = 3
DEFAULT_OCCURRENCE_CAP = 3

# Same wall as ``brute_count`` in cnf.py. Reduced formulas past this wall are
# recorded and not enumerated.
BRUTE_N_MAX = 20

# Bijection check enumerates both cubes. Above this it is left to the
# cardinality check, which itself stops at BRUTE_N_MAX.
EXTENSION_N_MAX = 12

# ``(window + degree cap + 2) * 2^n`` assignments. The cap is a safety rail
# for the test suite, not a property of the identity.
PIPELINE_ASSIGNMENT_BUDGET = 100_000

RECOVERY_FORMULA = "#F' = #F"


def _as_positive_cap(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"occurrence cap must be an int, got {value!r}")
    if value < MIN_OCCURRENCE_CAP:
        raise ValueError(
            f"occurrence cap {value} is below {MIN_OCCURRENCE_CAP}. "
            "Each copy sits in two cycle clauses, so this gadget needs a cap "
            "of at least 3 to leave a slot for an original clause."
        )
    return value


def copies_for_occurrence(occurrences: int, cap: int = DEFAULT_OCCURRENCE_CAP) -> int:
    """How many variables replace one source variable of the given occurrence.

    Unsplit variables, including those that occur zero times, contribute one
    coordinate. Split variables contribute ``ceil(occurrences / (cap - 2))``.
    """
    cap = _as_positive_cap(cap)
    if isinstance(occurrences, bool) or not isinstance(occurrences, int) or occurrences < 0:
        raise ValueError(f"occurrences must be a nonnegative int, got {occurrences!r}")
    if occurrences <= cap:
        return 1
    slots = cap - 2
    return (occurrences + slots - 1) // slots


def copying_lower_bound(total_occurrences: int, n: int) -> int:
    """Lower bound n' ≥ M − 2n for the cap-3 gadget.

    Each source variable contributes f(d) copies, with f(d) = d when d > 3 and
    f(d) = 1 otherwise. That satisfies f(d) ≥ d − 2, and summing gives the
    bound. It can be negative; the construction also has n' ≥ n.
    """
    if isinstance(total_occurrences, bool) or not isinstance(total_occurrences, int):
        raise ValueError(f"M must be an int, got {total_occurrences!r}")
    if isinstance(n, bool) or not isinstance(n, int):
        raise ValueError(f"n must be an int, got {n!r}")
    return total_occurrences - 2 * n


@dataclass(frozen=True)
class VariableCopies:
    """Images of one source variable under the copying map."""

    original: int
    occurrences: int
    copies: tuple[int, ...]

    @property
    def split(self) -> bool:
        return len(self.copies) > 1


@dataclass(frozen=True)
class Sparsification:
    """Source formula, reduced formula, and the copying table that relates them."""

    source: CNF
    reduced: CNF
    cap: int
    variables: tuple[VariableCopies, ...]
    cycle_clause_count: int
    max_occurrence_before: int
    max_occurrence_after: int
    literal_occurrences_before: int
    literal_occurrences_after: int

    @property
    def n_before(self) -> int:
        return self.source.n

    @property
    def m_before(self) -> int:
        return self.source.m

    @property
    def n_after(self) -> int:
        return self.reduced.n

    @property
    def m_after(self) -> int:
        return self.reduced.m

    @property
    def recovery(self) -> str:
        return RECOVERY_FORMULA


@dataclass(frozen=True)
class CountAgreement:
    """Brute-force comparison of #F and #F', or a recorded refusal."""

    status: str
    reason: str
    brute_source: int | None
    brute_reduced: int | None


@dataclass(frozen=True)
class Step1PipelineResult:
    """``interpolate_H0`` on the reduced formula, plus the occurrence caps."""

    status: str
    reason: str
    R: int
    k: int | None
    degree_cap: int | None
    degree_cap_mass: int | None
    degree_cap_loose: int | None
    observed_degree: int | None
    h0: Fraction | None
    h0_loose: Fraction | None
    direct_at_zero: Fraction | None
    brute: int | None
    assignments: int
    seconds: float


@dataclass(frozen=True)
class Step2PipelineResult:
    """Blocked sum on the reduced formula for a source block partition."""

    status: str
    reason: str
    r: int | None
    blocks: int | None
    tuples: int | None
    local_evaluations: int | None
    full_space: int | None
    brute: int | None
    blocked_at_zero: Fraction | None
    full_at_zero: Fraction | None
    blocked_at_one: Fraction | None
    seconds: float


def _profile(cnf: CNF) -> tuple[list[int], int, int, list[str]]:
    problems: list[str] = []
    for index, clause in enumerate(cnf.clauses):
        variables = [var for var, _neg in clause]
        if len(set(variables)) != len(variables):
            problems.append(f"clause {index} repeats a variable")
    counts = occurrence_counts(cnf)
    maximum = max(counts) if counts else 0
    mass = sum(len(clause) for clause in cnf.clauses)
    return counts, maximum, mass, problems


def _max_width(cnf: CNF) -> int:
    if not cnf.clauses:
        return 0
    return max(len(clause) for clause in cnf.clauses)


def sparsify(cnf: CNF, max_occurrence: int = DEFAULT_OCCURRENCE_CAP) -> Sparsification:
    """Copy variables so each occurs in at most ``max_occurrence`` clauses.

    The returned formula satisfies #F' = #F whenever both sides are counted.
    This function builds the formula either way. Enumeration is a separate
    step and stops at ``BRUTE_N_MAX``.
    """
    cap = _as_positive_cap(max_occurrence)
    counts, max_before, mass_before, problems = _profile(cnf)
    if problems:
        raise ValueError(
            "sparsify refuses a clause that repeats a variable "
            f"({problems[0]})"
        )
    slots = cap - 2
    variables: list[VariableCopies] = []
    cursor = 0
    for var, occ in enumerate(counts):
        n_copies = copies_for_occurrence(occ, cap)
        copies = tuple(range(cursor, cursor + n_copies))
        cursor += n_copies
        variables.append(VariableCopies(var, occ, copies))

    next_occurrence = [0] * cnf.n

    def take(var: int) -> int:
        info = variables[var]
        if not info.split:
            return info.copies[0]
        index = next_occurrence[var]
        next_occurrence[var] = index + 1
        copy_index = index // slots
        if copy_index >= len(info.copies):
            raise RuntimeError(
                f"variable {var} dealt occurrence {index} past {len(info.copies)} copies"
            )
        return info.copies[copy_index]

    rewritten: list[list[tuple[int, bool]]] = []
    for clause in cnf.clauses:
        rewritten.append([(take(var), neg) for var, neg in clause])

    cycles: list[list[tuple[int, bool]]] = []
    for info in variables:
        if not info.split:
            continue
        width = len(info.copies)
        for index, copy in enumerate(info.copies):
            nxt = info.copies[(index + 1) % width]
            cycles.append([(copy, True), (nxt, False)])

    reduced = CNF(cursor, rewritten + cycles)
    _check_copy_invariant(variables, rewritten, cycles, cap, slots)
    _counts_after, max_after, mass_after, reduced_problems = _profile(reduced)
    if reduced_problems or max_after > cap:
        detail = reduced_problems or [f"max occurrence {max_after} exceeds {cap}"]
        raise RuntimeError(f"copying invariant failed: {detail[0]}")
    if reduced.n < cnf.n:
        raise RuntimeError(f"reduced n={reduced.n} is below source n={cnf.n}")
    if cap == MIN_OCCURRENCE_CAP:
        lower = max(cnf.n, copying_lower_bound(mass_before, cnf.n))
        if reduced.n < lower:
            raise RuntimeError(
                f"reduced n={reduced.n} is below the cap-3 lower bound {lower}"
            )
    if mass_after > cap * reduced.n:
        raise RuntimeError(
            f"literal occurrences {mass_after} exceed {cap} per each of {reduced.n} variables"
        )
    return Sparsification(
        source=cnf,
        reduced=reduced,
        cap=cap,
        variables=tuple(variables),
        cycle_clause_count=len(cycles),
        max_occurrence_before=max_before,
        max_occurrence_after=max_after,
        literal_occurrences_before=mass_before,
        literal_occurrences_after=mass_after,
    )


def _check_copy_invariant(
    variables: Sequence[VariableCopies],
    rewritten: Sequence[Sequence[tuple[int, bool]]],
    cycles: Sequence[Sequence[tuple[int, bool]]],
    cap: int,
    slots: int,
) -> None:
    n_after = 0
    for info in variables:
        for copy in info.copies:
            n_after = max(n_after, copy + 1)
    original_uses = [0] * n_after
    cycle_uses = [0] * n_after
    for clause in rewritten:
        for var, _neg in clause:
            original_uses[var] += 1
    for clause in cycles:
        if len(clause) != 2:
            raise RuntimeError("cycle clause is not binary")
        for var, _neg in clause:
            cycle_uses[var] += 1
    for info in variables:
        if info.split:
            for copy in info.copies:
                if cycle_uses[copy] != 2:
                    raise RuntimeError(f"copy {copy} is in {cycle_uses[copy]} cycle clauses")
                if not 1 <= original_uses[copy] <= slots:
                    raise RuntimeError(
                        f"copy {copy} absorbed {original_uses[copy]} original clauses "
                        f"(slots {slots})"
                    )
                if original_uses[copy] + cycle_uses[copy] > cap:
                    raise RuntimeError(f"copy {copy} exceeds the occurrence cap")
        else:
            copy = info.copies[0]
            if cycle_uses[copy] != 0:
                raise RuntimeError(f"unsplit variable {info.original} received a cycle clause")
            if original_uses[copy] != info.occurrences:
                raise RuntimeError(
                    f"unsplit variable {info.original} occurrence changed under renaming"
                )


def rewritten_source_clauses(spec: Sparsification) -> list[list[tuple[int, bool]]]:
    """Reduced clauses that come from the source, before the cycle suffix."""
    if spec.cycle_clause_count == 0:
        return [list(clause) for clause in spec.reduced.clauses]
    stop = len(spec.reduced.clauses) - spec.cycle_clause_count
    return [list(clause) for clause in spec.reduced.clauses[:stop]]


def cycle_clauses(spec: Sparsification) -> list[list[tuple[int, bool]]]:
    """Implication-cycle clauses appended by ``sparsify``."""
    if spec.cycle_clause_count == 0:
        return []
    return [list(clause) for clause in spec.reduced.clauses[-spec.cycle_clause_count :]]


def extend_assignment(spec: Sparsification, bits: Sequence[int]) -> list[int]:
    """The unique copy-constant extension of a source assignment."""
    if len(bits) != spec.source.n:
        raise ValueError(
            f"assignment has length {len(bits)}, source formula has n={spec.source.n}"
        )
    out = [0] * spec.reduced.n
    for info in spec.variables:
        bit = bits[info.original]
        if bit not in (0, 1):
            raise ValueError(f"bit for variable {info.original} must be 0 or 1, got {bit!r}")
        for copy in info.copies:
            out[copy] = int(bit)
    return out


def project_assignment(spec: Sparsification, bits: Sequence[int]) -> list[int] | None:
    """Source assignment read off a reduced assignment whose copies agree.

    Returns ``None`` when some source variable's copies are not all equal.
    Those assignments falsify the cycle, so they are not models of the reduced
    formula.
    """
    if len(bits) != spec.reduced.n:
        raise ValueError(
            f"assignment has length {len(bits)}, reduced formula has n={spec.reduced.n}"
        )
    projected = [0] * spec.source.n
    for info in spec.variables:
        values = {bits[copy] for copy in info.copies}
        if len(values) != 1:
            return None
        bit = next(iter(values))
        if bit not in (0, 1):
            raise ValueError(f"bit for copy of variable {info.original} must be 0 or 1")
        projected[info.original] = int(bit)
    return projected


def extension_problems(spec: Sparsification, max_n: int = EXTENSION_N_MAX) -> list[str]:
    """Check the copy-constant bijection by enumerating both cubes.

    Raises ``ValueError`` above ``max_n``. An empty list means every source
    assignment satisfies the source formula exactly when its extension
    satisfies the reduced formula, and every reduced model is that extension
    of a source model.
    """
    if spec.source.n > max_n or spec.reduced.n > max_n:
        raise ValueError(
            f"extension check refuses n {spec.source.n}→{spec.reduced.n} above {max_n}"
        )
    problems: list[str] = []
    source_models = 0
    for bits in itertools.product((0, 1), repeat=spec.source.n):
        extended = extend_assignment(spec, bits)
        source_ok = is_satisfied(bits, spec.source)
        reduced_ok = is_satisfied(extended, spec.reduced)
        if source_ok != reduced_ok:
            problems.append(
                f"assignment {bits} satisfies the source ({source_ok}) but its "
                f"extension satisfies the reduced formula ({reduced_ok})"
            )
        if source_ok:
            source_models += 1
    reduced_models = 0
    for bits in itertools.product((0, 1), repeat=spec.reduced.n):
        if not is_satisfied(bits, spec.reduced):
            continue
        reduced_models += 1
        projected = project_assignment(spec, bits)
        if projected is None:
            problems.append(f"reduced model {bits} has unequal copies")
            continue
        if extend_assignment(spec, projected) != list(bits):
            problems.append(f"reduced model {bits} is not the constant extension of its projection")
        if not is_satisfied(projected, spec.source):
            problems.append(f"projection {projected} of reduced model {bits} is not a source model")
    if source_models != reduced_models and len(problems) < 8:
        problems.append(f"model counts differ: source {source_models}, reduced {reduced_models}")
    return problems


def count_agreement(spec: Sparsification) -> CountAgreement:
    """Compare brute-force model counts, or record that a side is past the wall.

    ``brute_count`` refuses n > 20 before it enumerates. The same cutoff is
    applied here. When the reduced side is past the wall, the source side is
    still counted if it fits, and the equality is reported as not checked.
    """
    source_count: int | None
    reduced_count: int | None
    if spec.source.n > BRUTE_N_MAX:
        source_count = None
    else:
        source_count = brute_count(spec.source)
    if spec.reduced.n > BRUTE_N_MAX:
        reduced_count = None
    else:
        reduced_count = brute_count(spec.reduced)
    if source_count is None or reduced_count is None:
        parts = []
        if source_count is None:
            parts.append(f"source n={spec.source.n}")
        if reduced_count is None:
            parts.append(f"reduced n={spec.reduced.n}")
        reason = (
            f"brute force stops at n={BRUTE_N_MAX}; "
            + " and ".join(parts)
            + f" (n {spec.n_before}→{spec.n_after}, m {spec.m_before}→{spec.m_after}). "
            f"The recovery {RECOVERY_FORMULA} was not checked on this instance."
        )
        return CountAgreement("refused", reason, source_count, reduced_count)
    if source_count != reduced_count:
        return CountAgreement(
            "disagreed",
            f"#F = {source_count} and #F' = {reduced_count} differ",
            source_count,
            reduced_count,
        )
    return CountAgreement(
        "agreed",
        RECOVERY_FORMULA,
        source_count,
        reduced_count,
    )


def _pipeline_assignments(n: int, degree_cap: int, loose_cap: int) -> int:
    """Cube assignments ``run_step1_pipeline`` will enumerate if it runs.

    The window is H(1), …, H(max(D, D_loose)+2). ``interpolate_H0`` then
    repeats H(1), …, H(D+1), and one more pass evaluates H(0) directly.
    """
    window = max(degree_cap, loose_cap) + 2
    return (window + degree_cap + 2) * (1 << n)


def _recover_at_cap(
    samples: Sequence[Fraction],
    degree_cap: int,
) -> tuple[Fraction, bool]:
    xs = [Fraction(t) for t in range(1, degree_cap + 2)]
    ys = list(samples[: degree_cap + 1])
    h0 = lagrange_eval_at_zero(xs, ys)
    predicted = lagrange_eval(xs, ys, Fraction(degree_cap + 2))
    return h0, predicted == samples[degree_cap + 1]


def run_step1_pipeline(
    spec: Sparsification,
    R: int = 1,
    rng: random.Random | None = None,
    budget: int = PIPELINE_ASSIGNMENT_BUDGET,
) -> Step1PipelineResult:
    """Recover H(0) of the reduced formula and compare it with #F'.

    Uses ``interpolate_H0`` at D = m·k·R, with k the maximum clause width, and
    a second interpolation at the loose occurrence cap B·n·R. For B = 3 that
    loose cap is ``degree_bound_3n``. The call refuses, without enumerating,
    when n' > 20 or when the planned assignment count exceeds ``budget``.
    """
    if isinstance(R, bool) or not isinstance(R, int) or R < 0:
        raise ValueError(f"R must be a nonnegative int, got {R!r}")
    reduced = spec.reduced
    k = _max_width(reduced)
    degree_cap = degree_bound(reduced, R=R, k=k)
    mass_cap = spec.literal_occurrences_after * R
    loose_cap = spec.cap * reduced.n * R
    if spec.cap == MIN_OCCURRENCE_CAP and loose_cap != degree_bound_3n(reduced.n, R):
        raise RuntimeError("loose cap drifted from degree_bound_3n")
    if mass_cap > loose_cap:
        raise RuntimeError(
            f"literal-occurrence cap {mass_cap} exceeds the per-variable cap {loose_cap}"
        )
    planned = _pipeline_assignments(reduced.n, degree_cap, loose_cap) if reduced.n <= 30 else budget + 1
    refused = _refuse_pipeline(spec, planned, budget, R, k, degree_cap, mass_cap, loose_cap)
    if refused is not None:
        return refused

    started = time.perf_counter()
    family = random_replacements(reduced, R=R, rng=rng or random.Random(0))
    window = max(degree_cap, loose_cap) + 2
    samples = [eval_H_at(reduced, family, Fraction(t)) for t in range(1, window + 1)]
    h0_width, extra_width = _recover_at_cap(samples, degree_cap)
    h0_loose, extra_loose = _recover_at_cap(samples, loose_cap)
    observed = equispaced_degree(samples)
    direct = eval_H_at(reduced, family, Fraction(0))
    interpolated, interpolated_cap, _family = interpolate_H0(reduced, fam=family, k=k)
    elapsed = time.perf_counter() - started
    brute = brute_count(reduced)
    problems: list[str] = []
    if interpolated_cap != degree_cap:
        problems.append(f"interpolate_H0 used cap {interpolated_cap}, expected {degree_cap}")
    if interpolated != h0_width:
        problems.append(f"interpolate_H0 returned {interpolated}, window returned {h0_width}")
    if h0_width != h0_loose:
        problems.append(f"width cap gave H(0)={h0_width}, loose cap gave H(0)={h0_loose}")
    if direct != h0_width:
        problems.append(f"direct H(0)={direct}, interpolated H(0)={h0_width}")
    if h0_width != Fraction(brute):
        problems.append(f"H(0)={h0_width} disagrees with #F'={brute}")
    if h0_width.denominator != 1:
        problems.append(f"H(0)={h0_width} is not an integer")
    if not extra_width:
        problems.append(f"extra node at the width cap {degree_cap} disagreed")
    if not extra_loose:
        problems.append(f"extra node at the loose cap {loose_cap} disagreed")
    if observed > mass_cap:
        problems.append(f"observed degree {observed} exceeds M·R = {mass_cap}")
    if observed > degree_cap or observed > loose_cap:
        problems.append(
            f"observed degree {observed} exceeds caps {degree_cap} and {loose_cap}"
        )
    status = "agreed" if not problems else "disagreed"
    reason = RECOVERY_FORMULA if not problems else "; ".join(problems)
    return Step1PipelineResult(
        status=status,
        reason=reason,
        R=R,
        k=k,
        degree_cap=degree_cap,
        degree_cap_mass=mass_cap,
        degree_cap_loose=loose_cap,
        observed_degree=observed,
        h0=h0_width,
        h0_loose=h0_loose,
        direct_at_zero=direct,
        brute=brute,
        assignments=planned,
        seconds=elapsed,
    )


def _refuse_pipeline(
    spec: Sparsification,
    planned: int,
    budget: int,
    R: int,
    k: int,
    degree_cap: int,
    mass_cap: int,
    loose_cap: int,
) -> Step1PipelineResult | None:
    if spec.reduced.n > BRUTE_N_MAX:
        return Step1PipelineResult(
            status="refused",
            reason=(
                f"reduced n={spec.reduced.n} exceeds the brute-force wall "
                f"{BRUTE_N_MAX}; interpolation was not started"
            ),
            R=R,
            k=k,
            degree_cap=degree_cap,
            degree_cap_mass=mass_cap,
            degree_cap_loose=loose_cap,
            observed_degree=None,
            h0=None,
            h0_loose=None,
            direct_at_zero=None,
            brute=None,
            assignments=0,
            seconds=0.0,
        )
    if planned > budget:
        return Step1PipelineResult(
            status="refused",
            reason=(
                f"interpolation would enumerate {planned} cube assignments, "
                f"above the budget {budget}; H(0) was not computed. "
                f"n'={spec.reduced.n}, D=m·k·R={degree_cap}, loose cap={loose_cap}."
            ),
            R=R,
            k=k,
            degree_cap=degree_cap,
            degree_cap_mass=mass_cap,
            degree_cap_loose=loose_cap,
            observed_degree=None,
            h0=None,
            h0_loose=None,
            direct_at_zero=None,
            brute=None,
            assignments=planned,
            seconds=0.0,
        )
    return None


def expand_blocks(
    spec: Sparsification,
    blocks: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    """Lift a partition of the source variables to the copied variables.

    Copies of a source variable stay in that variable's block, in the order
    the copies were allocated.
    """
    if not blocks:
        raise ValueError("need at least one block")
    owner = [-1] * spec.source.n
    for group, block in enumerate(blocks):
        if len(block) == 0:
            raise ValueError(f"block {group} is empty")
        for var in block:
            if isinstance(var, bool) or not isinstance(var, int):
                raise ValueError(f"block {group} contains non-int variable {var!r}")
            if not 0 <= var < spec.source.n:
                raise ValueError(
                    f"block {group} variable {var} is outside [0, {spec.source.n})"
                )
            if owner[var] != -1:
                raise ValueError(f"variable {var} is repeated in the partition")
            owner[var] = group
    missing = tuple(var for var, group in enumerate(owner) if group < 0)
    if missing:
        raise ValueError(
            f"blocks do not partition the source variables; missing {missing}"
        )
    expanded = []
    for block in blocks:
        copies: list[int] = []
        for var in block:
            copies.extend(spec.variables[var].copies)
        expanded.append(tuple(copies))
    return tuple(expanded)


def run_step2_pipeline(
    spec: Sparsification,
    blocks: Sequence[Sequence[int]],
    R: int = 1,
    budget: int = PIPELINE_ASSIGNMENT_BUDGET,
) -> Step2PipelineResult:
    """Blocked sum of the reduced formula when ``blocks`` partition the source.

    The design is the Step 2 popcount design on the expanded blocks, with
    every higher coefficient equal to 1. Sample points are t = 0 and t = 1.
    A clause that meets two blocks, a block above the Step 2 size rail, or a
    cube past the budget is recorded as a refusal and is not summed.
    """
    if isinstance(R, bool) or not isinstance(R, int) or R < 0:
        raise ValueError(f"R must be a nonnegative int, got {R!r}")
    expanded = expand_blocks(spec, blocks)
    reduced = spec.reduced
    if reduced.n > BRUTE_N_MAX:
        return _step2_refused(
            f"reduced n={reduced.n} exceeds the brute-force wall {BRUTE_N_MAX}; "
            "the blocked sum was not started"
        )
    # Two full-cube evaluations, plus the blocked enumeration beside them.
    planned_cube = 2 * (1 << reduced.n)
    if planned_cube > budget:
        return _step2_refused(
            f"blocked comparison would enumerate {planned_cube} cube assignments, "
            f"above the budget {budget}"
        )
    k = _max_width(reduced)
    design_cap = degree_bound(reduced, R=R, k=k)
    try:
        design = popcount_design(reduced.n, expanded, D=design_cap)
        validate_design(design, reduced)
    except ValueError as exc:
        return _step2_refused(str(exc))
    alpha = constant_alpha(reduced, R, design.tau, Fraction(1))
    try:
        structured = hand_built_replacement(reduced, design, R, alpha=alpha)
        started = time.perf_counter()
        comparisons = compare_blocked_to_full(reduced, structured, (Fraction(0), Fraction(1)))
        elapsed = time.perf_counter() - started
    except ValueError as exc:
        return _step2_refused(str(exc))
    problems: list[str] = []
    for t, full, blocked in comparisons:
        if blocked.value != full:
            problems.append(f"at t={t}, blocked {blocked.value} != full {full}")
    at_zero = comparisons[0]
    at_one = comparisons[1]
    brute = brute_count(reduced)
    if at_zero[0] != 0 or at_zero[2].value != Fraction(brute):
        problems.append(
            f"blocked H(0)={at_zero[2].value} disagrees with #F'={brute}"
        )
    sample = at_zero[2]
    status = "agreed" if not problems else "disagreed"
    reason = "blocked sum = full H(t), and H(0) = #F'" if not problems else "; ".join(problems)
    return Step2PipelineResult(
        status=status,
        reason=reason,
        r=sample.r,
        blocks=sample.blocks,
        tuples=sample.measurement_tuples,
        local_evaluations=sample.local_evaluations,
        full_space=sample.full_space,
        brute=brute,
        blocked_at_zero=at_zero[2].value,
        full_at_zero=at_zero[1],
        blocked_at_one=at_one[2].value,
        seconds=elapsed,
    )


def _step2_refused(reason: str) -> Step2PipelineResult:
    return Step2PipelineResult(
        status="refused",
        reason=reason,
        r=None,
        blocks=None,
        tuples=None,
        local_evaluations=None,
        full_space=None,
        brute=None,
        blocked_at_zero=None,
        full_at_zero=None,
        blocked_at_one=None,
        seconds=0.0,
    )


def _units(n_source_clauses: int, n_vars: int = 1, variable: int = 0) -> CNF:
    if n_vars < 1 or not 0 <= variable < n_vars:
        raise ValueError("unit formula dimensions are inconsistent")
    return CNF(n_vars, [[(variable, False)] for _ in range(n_source_clauses)])


def heavy_3cnf() -> CNF:
    """One variable in four width-3 clauses; the others occur at most three times."""
    return CNF(
        4,
        [
            [(0, False), (1, False), (2, False)],
            [(0, False), (1, False), (3, False)],
            [(0, False), (2, False), (3, False)],
            [(0, True), (1, True), (2, True)],
        ],
    )


def two_block_formula() -> CNF:
    """Width ≤ 2 formula whose clauses respect blocks {0, 1} and {2}.

    Variables 0 and 2 each occur four times, so both are copied. Variable 1
    occurs twice and is not copied.
    """
    return CNF(
        3,
        [
            [(0, False), (1, False)],
            [(0, False), (1, True)],
            [(0, False)],
            [(0, False)],
            [(2, False)],
            [(2, False)],
            [(2, False)],
            [(2, False)],
        ],
    )


def cross_block_formula() -> CNF:
    """A 3-clause on variables 0, 1, 2, plus three units on variable 0.

    Blocks {0, 1} and {2} split the 3-clause. Copying variable 0 does not
    pull variable 2 into its block.
    """
    return CNF(
        3,
        [
            [(0, False), (1, False), (2, False)],
            [(0, False)],
            [(0, False)],
            [(0, False)],
        ],
    )


@dataclass(frozen=True)
class Step5Case:
    """One source formula and the outcome the suite expects from it."""

    name: str
    cnf: CNF
    note: str
    expect_count: str
    expect_pipeline: str | None
    expect_blocks: str | None = None
    blocks: tuple[tuple[int, ...], ...] | None = None
    replacement_seed: int = 0

    def __post_init__(self) -> None:
        if self.expect_count not in {"agree", "refuse"}:
            raise ValueError(f"expect_count must be agree or refuse, got {self.expect_count!r}")
        if self.expect_pipeline not in {None, "agree", "refuse-budget", "refuse-brute"}:
            raise ValueError(f"unexpected pipeline expectation {self.expect_pipeline!r}")
        if self.expect_blocks not in {None, "agree", "refuse"}:
            raise ValueError(f"unexpected block expectation {self.expect_blocks!r}")
        if (self.expect_blocks is None) != (self.blocks is None):
            raise ValueError("block expectation and block partition must be set together")


@dataclass(frozen=True)
class Step5Report:
    """Count check, optional Step 1 pipeline, optional Step 2 pipeline."""

    name: str
    note: str
    sparsification: Sparsification
    count: CountAgreement
    pipeline: Step1PipelineResult | None
    blocked: Step2PipelineResult | None
    problems: tuple[str, ...]
    seconds: float

    @property
    def ok(self) -> bool:
        return not self.problems


def step5_cases() -> tuple[Step5Case, ...]:
    """Hand formulas, one seeded 3-CNF past the brute wall, and the block rows.

    Every case uses the default cap B = 3. Rows whose reduced n exceeds 20, or
    whose interpolation budget would, are expected refusals. They are part of
    the check: the suite passes when the refusal is the one recorded here.
    """
    positive_clause = CNF(3, [[(0, False), (1, False), (2, False)]])
    return (
        Step5Case(
            name="units-n1-m4",
            cnf=_units(4),
            note="Four positive units on one variable. Copies 4, cycle length 4, #F' = #F = 1.",
            expect_count="agree",
            expect_pipeline="agree",
            replacement_seed=701,
        ),
        Step5Case(
            name="free-bit-n2",
            cnf=_units(4, n_vars=2, variable=0),
            note="Four units on variable 0 and a free variable 1. Both sides count 2 models.",
            expect_count="agree",
            expect_pipeline="agree",
            replacement_seed=702,
        ),
        Step5Case(
            name="heavy-3cnf-n4",
            cnf=heavy_3cnf(),
            note="Width-3 formula, one variable in all four clauses. n 4→7, then interpolate_H0.",
            expect_count="agree",
            expect_pipeline="agree",
            replacement_seed=703,
        ),
        Step5Case(
            name="bounded-clause-n3",
            cnf=positive_clause,
            note="Already occurrence ≤ 3. The reduction renames nothing and adds no cycle.",
            expect_count="agree",
            expect_pipeline="agree",
            replacement_seed=704,
        ),
        Step5Case(
            name="empty-clause-n2",
            cnf=CNF(2, [[]]),
            note="An empty clause stays empty. Both counts are 0, and H is the zero polynomial.",
            expect_count="agree",
            expect_pipeline="agree",
            replacement_seed=705,
        ),
        Step5Case(
            name="two-block-copies",
            cnf=two_block_formula(),
            note=(
                "Clauses respect two source blocks. Copies stay in the owner's block. "
                "The blocked sum checks H(0) here. Degree-cap interpolation of copied "
                "formulas is the rows with n' ≤ 7."
            ),
            expect_count="agree",
            expect_pipeline=None,
            expect_blocks="agree",
            blocks=((0, 1), (2,)),
            replacement_seed=706,
        ),
        Step5Case(
            name="cross-block-copies",
            cnf=cross_block_formula(),
            note="The 3-clause meets both blocks. Copying does not repair that; Step 2 refuses.",
            expect_count="agree",
            expect_pipeline="agree",
            expect_blocks="refuse",
            blocks=((0, 1), (2,)),
            replacement_seed=707,
        ),
        Step5Case(
            name="units-n1-m16",
            cnf=_units(16),
            note="n' = 16 is inside the brute wall and outside the interpolation budget.",
            expect_count="agree",
            expect_pipeline="refuse-budget",
            replacement_seed=708,
        ),
        Step5Case(
            name="units-n1-m21",
            cnf=_units(21),
            note="n' = 21. brute_count refuses before enumerating. Equality is not checked.",
            expect_count="refuse",
            expect_pipeline="refuse-brute",
            replacement_seed=709,
        ),
        Step5Case(
            name="random-3cnf-n5-m12",
            cnf=random_3cnf(5, 12, rng=random.Random(0)),
            note="Exact 3-CNF. The cap-3 bound n' ≥ 3m − 2n = 26 already exceeds n = 20.",
            expect_count="refuse",
            expect_pipeline="refuse-brute",
            replacement_seed=710,
        ),
    )


def analyze_step5_case(
    case: Step5Case,
    budget: int = PIPELINE_ASSIGNMENT_BUDGET,
) -> Step5Report:
    """Run the expectations for one Step 5 case.

    An expected refusal is a passing report. A disagreement, or a refusal of
    a different kind than the case asked for, is a problem.
    """
    started = time.perf_counter()
    spec = sparsify(case.cnf)
    count = count_agreement(spec)
    pipeline: Step1PipelineResult | None = None
    if case.expect_pipeline is not None:
        pipeline = run_step1_pipeline(
            spec,
            R=1,
            rng=random.Random(case.replacement_seed),
            budget=budget,
        )
    blocked: Step2PipelineResult | None = None
    if case.blocks is not None:
        blocked = run_step2_pipeline(spec, case.blocks, R=1, budget=budget)
    problems: list[str] = []
    if spec.max_occurrence_after > spec.cap:
        problems.append(
            f"occurrence {spec.max_occurrence_after} exceeds cap {spec.cap}"
        )
    if (
        spec.source.n <= EXTENSION_N_MAX
        and spec.reduced.n <= EXTENSION_N_MAX
    ):
        problems.extend(extension_problems(spec))
    if case.expect_count == "agree" and count.status != "agreed":
        problems.append(f"count status {count.status}: {count.reason}")
    if case.expect_count == "refuse" and count.status != "refused":
        problems.append(f"expected a brute-force refusal, got {count.status}: {count.reason}")
    if case.expect_pipeline == "agree":
        if pipeline is None or pipeline.status != "agreed":
            reason = None if pipeline is None else pipeline.reason
            problems.append(f"pipeline status {None if pipeline is None else pipeline.status}: {reason}")
        elif count.brute_reduced is not None and pipeline.h0 != Fraction(count.brute_reduced):
            problems.append(f"pipeline H(0)={pipeline.h0} disagrees with #F'={count.brute_reduced}")
    elif case.expect_pipeline == "refuse-budget":
        if pipeline is None or pipeline.status != "refused" or "budget" not in pipeline.reason:
            problems.append(
                "expected an interpolation-budget refusal, got "
                f"{None if pipeline is None else pipeline.status}: "
                f"{None if pipeline is None else pipeline.reason}"
            )
    elif case.expect_pipeline == "refuse-brute":
        if pipeline is None or pipeline.status != "refused" or "brute-force wall" not in pipeline.reason:
            problems.append(
                "expected a brute-wall pipeline refusal, got "
                f"{None if pipeline is None else pipeline.status}"
            )
    if case.expect_blocks == "agree":
        if blocked is None or blocked.status != "agreed":
            problems.append(
                f"blocked status {None if blocked is None else blocked.status}: "
                f"{None if blocked is None else blocked.reason}"
            )
    elif case.expect_blocks == "refuse":
        if blocked is None or blocked.status != "refused" or "block" not in blocked.reason:
            problems.append(
                "expected Step 2 to refuse a clause that meets two blocks, got "
                f"{None if blocked is None else blocked.reason}"
            )
    elapsed = time.perf_counter() - started
    return Step5Report(
        name=case.name,
        note=case.note,
        sparsification=spec,
        count=count,
        pipeline=pipeline,
        blocked=blocked,
        problems=tuple(problems),
        seconds=elapsed,
    )
