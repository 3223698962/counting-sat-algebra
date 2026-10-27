"""
Simplified (r, τ, D)-design for a blocked sum.

This module does not build the matrix M0 and it does not implement Theorem 1.1.
The object below is the piece Step 2 needs in order to group variables, fix
shared measurements, and compare a blocked sum with the full sum H(t).

Parameters
----------
r
    Maximum block size. The variable set {0, ..., n-1} is partitioned into
    blocks B_0, ..., B_{b-1} with 1 ≤ |B_g| ≤ r.
τ
    Number of shared linear measurements. Measurement q has a rational weight
    vector w^{(q)} in Q^n. On a Boolean assignment x,

        μ_q(x) = Σ_i w^{(q)}_i x_i = Σ_g μ_q(x restricted to B_g).

    The second expression is the split a blocked sum uses: fix every block's
    contribution, add them to get the global measurement, then sum inside
    each block.
D
    Nonnegative integer carried for degree bookkeeping, the same role as the
    caps in Step 1. The blocked sum does not read D, does not interpolate,
    and does not assemble a matrix.

Clause restriction
------------------
A nonempty clause must take all of its variables from a single block. That is
what lets the summand factor once the measurements are fixed. A clause that
meets two blocks is rejected.

What is still missing for a paper-style design
----------------------------------------------
A design that could drive a faster counting algorithm would also have to keep
the measurement space small at large n and hand a low-degree matrix to Step 3.
Neither piece is in this module. For the popcount measurement (τ = 1 and every
weight equal to 1) the number of measurement tuples is Π_g (|B_g| + 1), which
grows exponentially with the number of blocks. ``blocked_sum.py`` records that
count; it does not claim 2^{o(n)}.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping, Sequence

from .cnf import CNF, Literal


def _as_index(value: object, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{what} must be an int, got {value!r}")
    return value


def _as_fraction(value: object, what: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise ValueError(f"{what} must be an int or Fraction, got {value!r}")
    return Fraction(value)


@dataclass(frozen=True)
class BlockDesign:
    """Partition into blocks of size at most r, plus τ shared linear measurements.

    Construct with ``make_design`` or ``popcount_design`` so the partition and
    the weight rows are checked. ``D`` is stored and otherwise unused.
    """

    n: int
    r: int
    tau: int
    D: int
    blocks: tuple[tuple[int, ...], ...]
    weights: tuple[tuple[Fraction, ...], ...]

    def owner_map(self) -> list[int]:
        """Block index of each variable, or -1 when the partition is incomplete."""
        owner = [-1] * self.n
        for group, block in enumerate(self.blocks):
            for var in block:
                if 0 <= var < self.n and owner[var] == -1:
                    owner[var] = group
                elif 0 <= var < self.n:
                    owner[var] = -2  # overlap marker; structural check reports it
        return owner

    def measure(self, assignment: Sequence[int]) -> tuple[Fraction, ...]:
        """Global measurement μ(x) for a full Boolean assignment."""
        if len(assignment) != self.n:
            raise ValueError(
                f"assignment has length {len(assignment)}, design has n={self.n}"
            )
        return tuple(
            sum(
                (weight[i] * assignment[i] for i in range(self.n)),
                start=Fraction(0),
            )
            for weight in self.weights
        )

    def local_measurement(
        self,
        block_index: int,
        assignment: Mapping[int, int],
    ) -> tuple[Fraction, ...]:
        """Contribution of one block. These add up to ``measure`` on a full assignment."""
        block = self.blocks[block_index]
        return tuple(
            sum((weight[var] * assignment[var] for var in block), start=Fraction(0))
            for weight in self.weights
        )


def structural_problems(design: BlockDesign) -> list[str]:
    """Partition, parameter, and weight-shape problems. Empty means the design is well formed."""
    problems: list[str] = []
    if isinstance(design.n, bool) or not isinstance(design.n, int) or design.n < 1:
        problems.append(f"n must be a positive int, got {design.n!r}")
        return problems
    if isinstance(design.r, bool) or not isinstance(design.r, int) or design.r < 1:
        problems.append(f"r must be a positive int, got {design.r!r}")
    if isinstance(design.tau, bool) or not isinstance(design.tau, int) or design.tau < 0:
        problems.append(f"tau must be a nonnegative int, got {design.tau!r}")
    if isinstance(design.D, bool) or not isinstance(design.D, int) or design.D < 0:
        problems.append(f"D must be a nonnegative int, got {design.D!r}")
    if not design.blocks:
        problems.append("need at least one block")

    seen: list[int] = []
    seen_set: set[int] = set()
    for group, block in enumerate(design.blocks):
        if len(block) < 1:
            problems.append(f"block {group} is empty")
        if len(block) > design.r:
            problems.append(
                f"block {group} has size {len(block)}, which is larger than r={design.r}"
            )
        for var in block:
            if isinstance(var, bool) or not isinstance(var, int):
                problems.append(f"block {group} contains non-int variable {var!r}")
                continue
            if not (0 <= var < design.n):
                problems.append(f"block {group} variable {var} is outside [0, {design.n})")
                continue
            if var in seen_set:
                problems.append(f"variable {var} is repeated in the partition")
            else:
                seen.append(var)
                seen_set.add(var)
    if len(seen) != design.n or seen_set != set(range(design.n)):
        missing = tuple(var for var in range(design.n) if var not in seen_set)
        if missing:
            problems.append(
                f"blocks do not partition {{0, ..., {design.n - 1}}}; missing {missing}"
            )

    if design.tau != len(design.weights):
        problems.append(
            f"tau={design.tau} does not match the number of weight vectors "
            f"({len(design.weights)})"
        )
    for row_index, row in enumerate(design.weights):
        if len(row) != design.n:
            problems.append(
                f"weight vector {row_index} has length {len(row)}, expected n={design.n}"
            )
            continue
        for entry in row:
            if not isinstance(entry, Fraction):
                problems.append(
                    f"weight vector {row_index} must contain Fraction entries, got {entry!r}"
                )
                break
    return problems


def clause_problems(design: BlockDesign, cnf: CNF) -> list[str]:
    """Reject clauses that do not sit inside a single block."""
    problems: list[str] = []
    if cnf.n != design.n:
        problems.append(f"formula has n={cnf.n}, design has n={design.n}")
        return problems
    owner = design.owner_map()
    for index, clause in enumerate(cnf.clauses):
        if not clause:
            continue
        groups: set[int] = set()
        for var, _neg in clause:
            if not (0 <= var < design.n) or owner[var] < 0:
                problems.append(f"clause {index} variable {var} is not in a unique block")
                groups = set()
                break
            groups.add(owner[var])
        if len(groups) > 1:
            joined = ", ".join(str(group) for group in sorted(groups))
            problems.append(
                f"clause {index} uses variables from blocks {joined}; "
                "this simplified design only sums clauses that sit inside one block"
            )
    return problems


def design_problems(design: BlockDesign, cnf: CNF | None = None) -> list[str]:
    """Structural problems, plus clause-support problems when a formula is given."""
    problems = structural_problems(design)
    if cnf is not None and not problems:
        problems.extend(clause_problems(design, cnf))
    return problems


def validate_design(design: BlockDesign, cnf: CNF | None = None) -> None:
    """Raise ``ValueError`` when the design, or its fit to ``cnf``, is illegal."""
    problems = design_problems(design, cnf)
    if problems:
        joined = "; ".join(problems)
        raise ValueError(f"illegal (r, τ, D)-design: {joined}")


def make_design(
    n: int,
    blocks: Sequence[Sequence[int]],
    weights: Sequence[Sequence[object]],
    D: int = 0,
    r: int | None = None,
) -> BlockDesign:
    """Build a design after checking the partition and the weight shape.

    ``r`` defaults to the size of the largest block. ``τ`` is ``len(weights)``.
    ``τ = 0`` is the empty measurement list: blocks share no linear measurement
    and the blocked sum is a single product of independent block sums.
    """
    n = _as_index(n, "n")
    D = _as_index(D, "D")
    if D < 0:
        raise ValueError(f"D must be a nonnegative int, got {D}")
    normalized_blocks: list[tuple[int, ...]] = []
    for group, block in enumerate(blocks):
        normalized_blocks.append(tuple(_as_index(var, f"block {group} variable") for var in block))
    if r is None:
        if not normalized_blocks:
            raise ValueError("need at least one block")
        r = max(len(block) for block in normalized_blocks)
    else:
        r = _as_index(r, "r")
    normalized_weights: list[tuple[Fraction, ...]] = []
    for row_index, row in enumerate(weights):
        normalized_weights.append(
            tuple(_as_fraction(entry, f"weight[{row_index}]") for entry in row)
        )
    design = BlockDesign(
        n=n,
        r=r,
        tau=len(normalized_weights),
        D=D,
        blocks=tuple(normalized_blocks),
        weights=tuple(normalized_weights),
    )
    validate_design(design)
    return design


def popcount_design(
    n: int,
    blocks: Sequence[Sequence[int]],
    D: int = 0,
    r: int | None = None,
) -> BlockDesign:
    """One shared measurement: the sum of all variables (every weight is 1).

    A block of size s then has s + 1 possible local measurements, the popcounts
    0, ..., s. Across b blocks the tuple count is Π_g (|B_g| + 1).
    """
    weights = [Fraction(1) for _ in range(n)]
    return make_design(n, blocks, [weights], D=D, r=r)


def clause_block_index(design: BlockDesign, clause: Sequence[Literal]) -> int | None:
    """Unique block of a nonempty clause.

    Returns ``None`` for an empty clause (the clause polynomial is identically
    0, so the whole product is 0). Raises ``ValueError`` when the clause meets
    more than one block.
    """
    if not clause:
        return None
    owner = design.owner_map()
    groups: set[int] = set()
    for var, _neg in clause:
        if not (0 <= var < design.n) or owner[var] < 0:
            raise ValueError(f"variable {var} is not in a unique block")
        groups.add(owner[var])
    if len(groups) != 1:
        raise ValueError(
            "clause meets more than one block; the simplified blocked sum "
            "requires each clause to sit inside a single block"
        )
    return next(iter(groups))
