"""
Blocked sum for a hand-built replacement family on a simplified (r, τ, D)-design.

The full sum in ``replacements.eval_H_at`` is

    H(t) = Σ_{x in {0,1}^n}  Π_e Φ_e( ⟨P_{e,v}(t), x⟩ ).

Here each replacement is hand-built so that its non-constant part is a
polynomial in the design's shared measurements:

    ⟨P_{e,v}(t), X⟩ = X_v + Σ_{j=1}^{R} t^j · ⟨α_{e,v,j}, μ(X)⟩,

with μ(X) in Q^τ the measurement vector of the design. The constant term is
the standard basis vector e_v, so the form still restores X_v at t = 0.

When every nonempty clause sits inside one block of the design, the summand
depends on a block's own bits and on the global measurement only. Grouping
Boolean assignments by the local measurement λ_g of each block gives the identity

    H(t) = Σ_{λ_1, ..., λ_b}  Π_g A_g(λ_g, λ_1 + ... + λ_b),

    A_g(λ, m) = Σ_{y in {0,1}^{B_g} : μ(y)=λ}  Π_{e in block g} Φ_e(y, m).

The sum on the right is the blocked sum. It is an exact regrouping of a finite
sum. It is not a subexponential algorithm.

Cost, for this implementation
-----------------------------
Every realizable local measurement is enumerated. The number of tuples is
Π_g |Λ_g|, where Λ_g is the set of measurement values attained on block g.
That quantity is at most 2^n and, for the popcount design, exactly
Π_g (|B_g| + 1). With b blocks of size r the popcount tuple count is (r + 1)^b,
which is exponential in the number of groups. Local clause products are
evaluated b · 2^r · (r + 1)^{b-1} times in that same case. Both counts are
returned on ``BlockedSumResult`` so a caller can see the growth. Wall-clock
time on small n is also recorded; it is one machine's measurement, not a
complexity bound.

D is carried on the design and is not consulted. The matrix M0 is not built.
Theorem 1.1 is not implemented. Nothing here runs in 2^{o(n)}.
"""

from __future__ import annotations

import itertools
import random
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .cnf import CNF, brute_count
from .design import (
    BlockDesign,
    clause_block_index,
    make_design,
    popcount_design,
    validate_design,
)
from .replacements import (
    ReplacementFamily,
    degree_bound,
    eval_H_at,
    phi_clause_numeric,
)

# Small rationals, same set as the random replacements in Step 0.
_ALPHA_CHOICES: tuple[Fraction, ...] = (
    Fraction(-2),
    Fraction(-1),
    Fraction(-1, 2),
    Fraction(0),
    Fraction(1, 2),
    Fraction(1),
    Fraction(2),
)

# Guards so a mistaken call cannot start an enormous measurement loop.
# The identity is unchanged below the cap; the cap is only a safety rail.
DEFAULT_MAX_BLOCK_SIZE = 12
DEFAULT_MAX_MEASUREMENT_TUPLES = 100_000


def _as_fraction(value: object, what: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise ValueError(f"{what} must be an int or Fraction, got {value!r}")
    return Fraction(value)


# alpha[clause][occurrence][power index j][measurement q]
# multiplies t^{j+1} and measurement q.
AlphaTable = tuple[tuple[tuple[tuple[Fraction, ...], ...], ...], ...]


@dataclass
class StructuredReplacement:
    """Hand-built family whose higher T-coefficients lie in the measurement span.

    ``family`` is the coefficient representation consumed by ``eval_H_at``.
    ``alpha`` is the same data in measurement coordinates, which is what the
    blocked sum evaluates. The constructor checks that the two agree on a few
    assignments, including the restoration ⟨P(0), x⟩ = x_v.
    """

    family: ReplacementFamily
    design: BlockDesign
    alpha: AlphaTable
    R: int


@dataclass(frozen=True)
class BlockedSumResult:
    """One blocked evaluation, together with the enumeration counts that dominate its cost."""

    value: Fraction
    n: int
    blocks: int
    r: int
    tau: int
    degree_cap: int
    measurement_tuples: int
    local_evaluations: int
    full_space: int
    seconds: float


@dataclass(frozen=True)
class BlockedExample:
    """Named hand-built instance used by tests and by the verify CLI."""

    name: str
    cnf: CNF
    structured: StructuredReplacement
    sample_points: tuple[Fraction, ...]
    notes: str


@dataclass(frozen=True)
class GrowthRow:
    """Blocked sum against the full sum for one group-count in the popcount scan."""

    groups: int
    blocked: BlockedSumResult
    full: Fraction

    @property
    def ok(self) -> bool:
        return self.blocked.value == self.full


def measurement_linear_form(
    structured: StructuredReplacement,
    clause_index: int,
    occurrence: int,
    bit: int,
    measurement: Sequence[Fraction],
    t: Fraction,
) -> Fraction:
    """⟨P, X⟩ = bit + Σ_{j≥1} t^j ⟨α_j, μ⟩, evaluated by Horner in the powers of t."""
    alpha_occ = structured.alpha[clause_index][occurrence]
    if len(measurement) != structured.design.tau:
        raise ValueError(
            f"measurement has length {len(measurement)}, design has tau={structured.design.tau}"
        )
    extra = Fraction(0)
    for power in reversed(range(structured.R)):
        dot = Fraction(0)
        row = alpha_occ[power]
        for q, coeff in enumerate(row):
            dot += coeff * measurement[q]
        extra = extra * t + dot
    return Fraction(bit) + extra * t


def _normalize_alpha(alpha: Sequence[Sequence[Sequence[Sequence[object]]]], cnf: CNF, R: int, tau: int) -> AlphaTable:
    if R < 0:
        raise ValueError(f"R must be nonnegative, got {R}")
    if len(alpha) != cnf.m:
        raise ValueError(f"alpha has {len(alpha)} clauses, formula has {cnf.m}")
    normalized_clauses = []
    for clause_index, clause in enumerate(cnf.clauses):
        clause_alpha = alpha[clause_index]
        if len(clause_alpha) != len(clause):
            raise ValueError(
                f"alpha[{clause_index}] has {len(clause_alpha)} occurrences, "
                f"clause has {len(clause)}"
            )
        normalized_occs = []
        for occ in range(len(clause)):
            powers = clause_alpha[occ]
            if len(powers) != R:
                raise ValueError(
                    f"alpha[{clause_index}][{occ}] has {len(powers)} powers, expected R={R}"
                )
            normalized_powers = []
            for power in range(R):
                row = powers[power]
                if len(row) != tau:
                    raise ValueError(
                        f"alpha[{clause_index}][{occ}][{power}] has length {len(row)}, "
                        f"expected tau={tau}"
                    )
                normalized_powers.append(
                    tuple(
                        _as_fraction(entry, f"alpha[{clause_index}][{occ}][{power}][{q}]")
                        for q, entry in enumerate(row)
                    )
                )
            normalized_occs.append(tuple(normalized_powers))
        normalized_clauses.append(tuple(normalized_occs))
    return tuple(normalized_clauses)


def _random_alpha(cnf: CNF, R: int, tau: int, rng: random.Random) -> AlphaTable:
    raw = [
        [
            [[rng.choice(_ALPHA_CHOICES) for _q in range(tau)] for _power in range(R)]
            for _occ in clause
        ]
        for clause in cnf.clauses
    ]
    return _normalize_alpha(raw, cnf, R, tau)


def _coefficient_rows(
    cnf: CNF,
    design: BlockDesign,
    R: int,
    alpha: AlphaTable,
) -> tuple[list[list[list[list[Fraction]]]], list[list[int]], list[list[bool]]]:
    rows: list[list[list[list[Fraction]]]] = []
    occ_vars: list[list[int]] = []
    occ_neg: list[list[bool]] = []
    for clause_index, clause in enumerate(cnf.clauses):
        clause_rows: list[list[list[Fraction]]] = []
        variables: list[int] = []
        negated: list[bool] = []
        for occ, (var, neg) in enumerate(clause):
            row: list[list[Fraction]] = []
            for i in range(cnf.n):
                coeffs = [Fraction(0) for _ in range(R + 1)]
                if i == var:
                    coeffs[0] = Fraction(1)
                for power in range(R):
                    coeff = Fraction(0)
                    for q in range(design.tau):
                        coeff += alpha[clause_index][occ][power][q] * design.weights[q][i]
                    coeffs[power + 1] = coeff
                row.append(coeffs)
            clause_rows.append(row)
            variables.append(var)
            negated.append(neg)
        rows.append(clause_rows)
        occ_vars.append(variables)
        occ_neg.append(negated)
    return rows, occ_vars, occ_neg


def _representations_agree(structured: StructuredReplacement) -> None:
    """Lock the coefficient rows to the measurement formula on a few assignments."""
    n = structured.design.n
    probes = [tuple(0 for _ in range(n)), tuple(1 for _ in range(n))]
    if n >= 2:
        probes.append(tuple(i % 2 for i in range(n)))
    points = (Fraction(0), Fraction(1), Fraction(2))
    for bits in probes:
        measurement = structured.design.measure(bits)
        for clause_index, clause in enumerate(structured.family.occ_vars):
            for occ, var in enumerate(clause):
                for t in points:
                    from_rows = structured.family.linear_form(clause_index, occ, bits, t)
                    from_measurements = measurement_linear_form(
                        structured,
                        clause_index,
                        occ,
                        bits[var],
                        measurement,
                        t,
                    )
                    if from_rows != from_measurements:
                        raise RuntimeError(
                            "replacement coefficients do not match the measurement formula "
                            f"at clause {clause_index} occurrence {occ}"
                        )
                    if t == 0 and from_rows != Fraction(bits[var]):
                        raise RuntimeError(
                            f"T=0 failed to restore X_{var} "
                            f"(clause {clause_index} occurrence {occ})"
                        )


def hand_built_replacement(
    cnf: CNF,
    design: BlockDesign,
    R: int,
    alpha: Sequence[Sequence[Sequence[Sequence[object]]]] | None = None,
    rng: random.Random | None = None,
) -> StructuredReplacement:
    """Replacement family with ⟨P, X⟩ = X_v + Σ_{j≥1} t^j ⟨α_j, μ(X)⟩.

    Pass ``alpha`` for a fully specified family. When ``alpha`` is omitted the
    coefficients are drawn from the Step 0 rational set; the default generator
    is ``Random(0)`` so the draw is reproducible.
    """
    if isinstance(R, bool) or not isinstance(R, int) or R < 0:
        raise ValueError(f"R must be a nonnegative int, got {R!r}")
    validate_design(design, cnf)
    if alpha is None:
        alpha_table = _random_alpha(cnf, R, design.tau, rng or random.Random(0))
    else:
        alpha_table = _normalize_alpha(alpha, cnf, R, design.tau)
    rows, occ_vars, occ_neg = _coefficient_rows(cnf, design, R, alpha_table)
    family = ReplacementFamily(
        n=cnf.n,
        R=R,
        rows=rows,
        occ_vars=occ_vars,
        occ_neg=occ_neg,
    )
    structured = StructuredReplacement(
        family=family,
        design=design,
        alpha=alpha_table,
        R=R,
    )
    _representations_agree(structured)
    return structured


def constant_alpha(cnf: CNF, R: int, tau: int, value: Fraction = Fraction(1)) -> AlphaTable:
    """Every measurement coefficient equal to ``value``. Still restores X_v at t = 0."""
    if isinstance(R, bool) or not isinstance(R, int) or R < 0:
        raise ValueError(f"R must be a nonnegative int, got {R!r}")
    if isinstance(tau, bool) or not isinstance(tau, int) or tau < 0:
        raise ValueError(f"tau must be a nonnegative int, got {tau!r}")
    raw = [
        [
            [[value for _q in range(tau)] for _power in range(R)]
            for _occ in clause
        ]
        for clause in cnf.clauses
    ]
    return _normalize_alpha(raw, cnf, R, tau)


def _local_product(
    cnf: CNF,
    structured: StructuredReplacement,
    clause_ids: Sequence[int],
    assignment: dict[int, int],
    measurement: tuple[Fraction, ...],
    t: Fraction,
) -> Fraction:
    product = Fraction(1)
    for clause_index in clause_ids:
        clause = cnf.clauses[clause_index]
        values = [
            measurement_linear_form(
                structured,
                clause_index,
                occ,
                assignment[var],
                measurement,
                t,
            )
            for occ, (var, _neg) in enumerate(clause)
        ]
        product *= phi_clause_numeric(values, structured.family.occ_neg[clause_index])
        if product == 0:
            break
    return product


def blocked_sum(
    cnf: CNF,
    structured: StructuredReplacement,
    t: Fraction | int,
    max_block_size: int = DEFAULT_MAX_BLOCK_SIZE,
    max_measurement_tuples: int = DEFAULT_MAX_MEASUREMENT_TUPLES,
) -> BlockedSumResult:
    """Evaluate the blocked-sum regrouping of H(t).

    Raises ``ValueError`` when the formula is not legal for the design, when a
    block is larger than ``max_block_size``, or when the measurement tuple
    count would exceed ``max_measurement_tuples``. The last two limits are
    safety rails. Hitting them means this exponential enumeration got too big
    for the caller, which is the cost this step is willing to report rather
    than hide.
    """
    if structured.family.n != cnf.n or structured.design.n != cnf.n:
        raise ValueError("formula, replacement family, and design disagree on n")
    if structured.family.R != structured.R:
        raise ValueError("replacement degree does not match the structured family")
    validate_design(structured.design, cnf)
    point = Fraction(t)
    design = structured.design
    for group, block in enumerate(design.blocks):
        if len(block) > max_block_size:
            raise ValueError(
                f"block {group} has size {len(block)}, above the enumeration cap "
                f"{max_block_size}"
            )

    if any(len(clause) == 0 for clause in cnf.clauses):
        # Φ of an empty clause is 0, so H(t) is 0. Skip the measurement loop.
        return BlockedSumResult(
            value=Fraction(0),
            n=cnf.n,
            blocks=len(design.blocks),
            r=design.r,
            tau=design.tau,
            degree_cap=design.D,
            measurement_tuples=0,
            local_evaluations=0,
            full_space=2**cnf.n,
            seconds=0.0,
        )

    clause_ids: list[list[int]] = [[] for _ in design.blocks]
    for clause_index, clause in enumerate(cnf.clauses):
        group = clause_block_index(design, clause)
        if group is None:
            continue
        clause_ids[group].append(clause_index)

    buckets: list[dict[tuple[Fraction, ...], list[dict[int, int]]]] = []
    for group, block in enumerate(design.blocks):
        bucket: dict[tuple[Fraction, ...], list[dict[int, int]]] = {}
        for bits in itertools.product((0, 1), repeat=len(block)):
            assignment = {var: bit for var, bit in zip(block, bits)}
            local = design.local_measurement(group, assignment)
            bucket.setdefault(local, []).append(assignment)
        buckets.append(bucket)

    tuple_count = 1
    for bucket in buckets:
        tuple_count *= len(bucket)
    if tuple_count > max_measurement_tuples:
        raise ValueError(
            f"measurement enumeration would cover {tuple_count} tuples, "
            f"above the cap {max_measurement_tuples}. This simplified blocked sum "
            "is exponential in the number of blocks; the call refuses this instance"
        )

    keys = [list(bucket.keys()) for bucket in buckets]
    started = time.perf_counter()
    total = Fraction(0)
    seen_tuples = 0
    local_evaluations = 0
    for local_tuple in itertools.product(*keys):
        seen_tuples += 1
        measurement = tuple(
            sum((local_tuple[group][q] for group in range(len(buckets))), start=Fraction(0))
            for q in range(design.tau)
        )
        term = Fraction(1)
        for group, local in enumerate(local_tuple):
            block_sum = Fraction(0)
            for assignment in buckets[group][local]:
                local_evaluations += 1
                block_sum += _local_product(
                    cnf,
                    structured,
                    clause_ids[group],
                    assignment,
                    measurement,
                    point,
                )
            term *= block_sum
        total += term
    elapsed = time.perf_counter() - started
    if seen_tuples != tuple_count:
        raise RuntimeError(
            f"enumerated {seen_tuples} measurement tuples, expected {tuple_count}"
        )
    return BlockedSumResult(
        value=total,
        n=cnf.n,
        blocks=len(design.blocks),
        r=design.r,
        tau=design.tau,
        degree_cap=design.D,
        measurement_tuples=seen_tuples,
        local_evaluations=local_evaluations,
        full_space=2**cnf.n,
        seconds=elapsed,
    )


def compare_blocked_to_full(
    cnf: CNF,
    structured: StructuredReplacement,
    points: Sequence[Fraction | int],
) -> list[tuple[Fraction, Fraction, BlockedSumResult]]:
    """Return ``(t, full H(t), blocked result)`` for each sample point."""
    rows: list[tuple[Fraction, Fraction, BlockedSumResult]] = []
    for point in points:
        t = Fraction(point)
        full = eval_H_at(cnf, structured.family, t)
        blocked = blocked_sum(cnf, structured, t)
        rows.append((t, full, blocked))
    return rows


def popcount_operation_counts(groups: int, block_size: int) -> tuple[int, int, int]:
    """Exact enumeration sizes for unit popcount weights and equal blocks.

    Returns ``(measurement_tuples, local_evaluations, full_space)``.

    Every popcount 0, ..., block_size occurs on a block, so each block
    contributes ``block_size + 1`` measurement values. The figures do not
    depend on the clauses or on α. They describe this implementation, which
    walks every realizable measurement even when a coefficient happens to be 0.
    """
    if isinstance(groups, bool) or not isinstance(groups, int) or groups < 1:
        raise ValueError(f"groups must be a positive int, got {groups!r}")
    if isinstance(block_size, bool) or not isinstance(block_size, int) or block_size < 1:
        raise ValueError(f"block_size must be a positive int, got {block_size!r}")
    tuples = (block_size + 1) ** groups
    local = groups * (2**block_size) * (block_size + 1) ** (groups - 1)
    full = 2 ** (groups * block_size)
    return tuples, local, full


def _group_clause(start: int, block_size: int) -> list[tuple[int, bool]]:
    width = min(block_size, 3)
    return [(start + offset, (start + offset) % 2 == 1) for offset in range(width)]


def popcount_group_instance(
    num_groups: int,
    block_size: int = 2,
    R: int = 1,
) -> tuple[CNF, StructuredReplacement]:
    """Deterministic popcount instance with one clause per block and α = 1.

    Group g owns variables ``[g * block_size, (g + 1) * block_size)``. The
    clause on that group uses the first ``min(block_size, 3)`` variables, with
    the odd variables negated. Degree cap stored on the design is
    ``m * width * R`` and is not used by the sum.
    """
    if num_groups < 1 or block_size < 1:
        raise ValueError("num_groups and block_size must be positive")
    if block_size > DEFAULT_MAX_BLOCK_SIZE:
        raise ValueError(f"block_size {block_size} exceeds the enumeration cap")
    n = num_groups * block_size
    blocks = tuple(
        tuple(range(group * block_size, (group + 1) * block_size))
        for group in range(num_groups)
    )
    clauses = [_group_clause(group * block_size, block_size) for group in range(num_groups)]
    cnf = CNF(n, clauses)
    width = min(block_size, 3)
    design = popcount_design(
        n,
        blocks,
        D=degree_bound(cnf, R=R, k=width),
        r=block_size,
    )
    structured = hand_built_replacement(
        cnf,
        design,
        R=R,
        alpha=constant_alpha(cnf, R, tau=1, value=Fraction(1)),
    )
    return cnf, structured


def popcount_growth(
    group_counts: Sequence[int],
    block_size: int = 2,
    R: int = 1,
    t: Fraction | int = 1,
) -> list[GrowthRow]:
    """Blocked sum versus full H(t) as the number of equal popcount blocks grows."""
    rows: list[GrowthRow] = []
    point = Fraction(t)
    for groups in group_counts:
        cnf, structured = popcount_group_instance(groups, block_size=block_size, R=R)
        blocked = blocked_sum(cnf, structured, point)
        full = eval_H_at(cnf, structured.family, point)
        rows.append(GrowthRow(groups=groups, blocked=blocked, full=full))
    return rows


def _varied_unit_alpha(cnf: CNF, values: Sequence[Fraction], tau: int, R: int) -> AlphaTable:
    """Cycle ``values`` across occurrences. Each power repeats that occurrence's row.

    Only measurement 0 receives the value; any further measurement coefficients are 0.
    """
    occ_index = 0
    raw_clauses = []
    for clause in cnf.clauses:
        raw_occs = []
        for _var, _neg in clause:
            value = values[occ_index % len(values)]
            occ_index += 1
            row = tuple(value if q == 0 else Fraction(0) for q in range(tau))
            raw_occs.append(tuple(tuple(row) for _power in range(R)))
        raw_clauses.append(raw_occs)
    return _normalize_alpha(raw_clauses, cnf, R, tau)


def _explicit_alpha_tau2(cnf: CNF) -> AlphaTable:
    """R = 2, τ = 2, coefficients written out per occurrence parity."""
    raw = []
    occurrence = 0
    for clause in cnf.clauses:
        clause_alpha = []
        for _var, _neg in clause:
            if occurrence % 2 == 0:
                powers = (
                    (Fraction(1), Fraction(0)),
                    (Fraction(0), Fraction(1)),
                )
            else:
                powers = (
                    (Fraction(-1, 2), Fraction(1, 2)),
                    (Fraction(1), Fraction(-1)),
                )
            clause_alpha.append(powers)
            occurrence += 1
        raw.append(clause_alpha)
    return _normalize_alpha(raw, cnf, R=2, tau=2)


def hand_built_examples() -> list[BlockedExample]:
    """Small hand-built families. Every blocked sum in the suite is checked against full H(t)."""
    one = CNF(3, [[(0, False), (1, False), (2, True)]])
    one_design = popcount_design(3, [(0, 1, 2)], D=degree_bound(one, R=1, k=3), r=3)
    one_alpha = (
        (
            ((Fraction(1),),),
            ((Fraction(1, 2),),),
            ((Fraction(-1),),),
        ),
    )
    one_structured = hand_built_replacement(one, one_design, R=1, alpha=one_alpha)

    two = CNF(
        6,
        [
            [(0, False), (1, False), (2, True)],
            [(0, True), (1, True), (2, False)],
            [(3, False), (4, True), (5, False)],
        ],
    )
    two_design = popcount_design(
        6,
        [(0, 1, 2), (3, 4, 5)],
        D=degree_bound(two, R=1, k=3),
        r=3,
    )
    two_structured = hand_built_replacement(
        two,
        two_design,
        R=1,
        alpha=constant_alpha(two, R=1, tau=1, value=Fraction(1)),
    )

    three_cnf, three_structured = popcount_group_instance(3, block_size=2, R=1)
    # Replace the constant-α family with an explicit alternating table so this
    # example is not identical to the growth scan's α = 1 family.
    three_design = three_structured.design
    three_alpha = _varied_unit_alpha(
        three_cnf,
        (Fraction(1), Fraction(-1, 2)),
        tau=1,
        R=1,
    )
    three_structured = hand_built_replacement(three_cnf, three_design, R=1, alpha=three_alpha)

    flagged_weights = (
        tuple(Fraction(1) for _ in range(6)),
        (Fraction(1), Fraction(0), Fraction(0), Fraction(1), Fraction(0), Fraction(0)),
    )
    flagged = make_design(
        6,
        [(0, 1, 2), (3, 4, 5)],
        flagged_weights,
        D=degree_bound(two, R=2, k=3),
        r=3,
    )
    flagged_structured = hand_built_replacement(two, flagged, R=2, alpha=_explicit_alpha_tau2(two))

    return [
        BlockedExample(
            name="one-block-n3",
            cnf=one,
            structured=one_structured,
            sample_points=(Fraction(0), Fraction(1), Fraction(2)),
            notes="single block, τ=1 popcount, R=1, occurrence coefficients 1, 1/2, -1",
        ),
        BlockedExample(
            name="two-block-3cnf-n6",
            cnf=two,
            structured=two_structured,
            sample_points=(Fraction(0), Fraction(1), Fraction(2), Fraction(3, 2)),
            notes="two blocks of size 3, τ=1 popcount, R=1, every α coefficient 1",
        ),
        BlockedExample(
            name="three-block-n6",
            cnf=three_cnf,
            structured=three_structured,
            sample_points=(Fraction(0), Fraction(1), Fraction(2)),
            notes="three blocks of size 2, τ=1 popcount, R=1, α alternating 1 and -1/2",
        ),
        BlockedExample(
            name="two-block-tau2-R2-n6",
            cnf=two,
            structured=flagged_structured,
            sample_points=(Fraction(0), Fraction(1), Fraction(2)),
            notes="two blocks, τ=2 (popcount and first-variable flag), R=2, written-out α",
        ),
    ]


def count_at_zero_matches_brute(cnf: CNF, structured: StructuredReplacement) -> bool:
    """H(0) from the blocked sum equals the brute-force model count.

    This is the T = 0 restoration on this family, checked through the blocked
    regrouping. It is the same identity Step 0 checks by full enumeration.
    """
    blocked = blocked_sum(cnf, structured, Fraction(0))
    brute = brute_count(cnf)
    return blocked.value == Fraction(brute) and blocked.value.denominator == 1
