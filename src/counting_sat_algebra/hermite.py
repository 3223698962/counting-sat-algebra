"""
Exact Hermite recovery of H(0) from Hasse coefficients.

H is univariate of degree at most D = m·k·R. The Step 0 path in
``replacements.interpolate_H0`` reads the value H(t) at the D+1 nodes
t = 1, …, D+1 and returns the Lagrange value at 0. This module reads several
Hasse coefficients at each node,

    H^{[k]}(t) = (1/k!) H^{(k)}(t),

so the same D+1 conditions determine H from fewer distinct nodes. The Newton
form of confluent interpolation uses those Hasse coefficients as the
divided differences on repeated nodes. Taylor coefficients are k! times the
Hasse coefficients; ``taylor_from_hasse`` converts them.

The jet at one node is the formal derivative of the clause product. For each
assignment x the linear forms are polynomials in T, the clause factor
Φ_e = 1 − ∏(1 − L_v) is evaluated in the truncated ring Q[u]/(u^{order+1})
with u = T − t, and the jets are summed. One pass over {0,1}^n returns every
requested coefficient at that t. Derivative orders are not separate passes
over the cube.

What this module does not do
----------------------------
The number of conditions stays D+1. A degree-D polynomial is not determined
by fewer conditions than that; the saving is only in the number of distinct
sample points. Each of those points still enumerates 2^n assignments. A full
jet of order D does Θ(D²) arithmetic per assignment, so one node that carries
every coefficient can do more rational arithmetic than the D+1 scalar
evaluations it replaces. Nothing here builds the matrix M0, implements
Theorem 1.1, or runs in 2^{o(n)} time.
"""

from __future__ import annotations

import itertools
import random
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .bounded_occurrence import bench_path, cnf_from_record, load_bench_document
from .cnf import CNF, brute_count, random_3cnf
from .replacements import (
    ReplacementFamily,
    degree_bound,
    eval_H_at,
    interpolate_H0,
    random_replacements,
)


def polynomial_hasse_jet(
    coeffs: Sequence[Fraction],
    t: Fraction,
    order: int,
) -> tuple[Fraction, ...]:
    """Hasse coefficients of a univariate polynomial at ``t``, indices 0..order.

    ``coeffs[j]`` is the coefficient of T^j. The result ``c`` satisfies
    p(t + u) = Σ_k c[k] u^k + O(u^{order+1}).
    """
    if order < 0:
        raise ValueError("order must be nonnegative")
    point = Fraction(t)
    acc = [Fraction(0)] * (order + 1)
    for raw in reversed(coeffs):
        updated = [Fraction(0)] * (order + 1)
        for k in range(order + 1):
            updated[k] = point * acc[k]
            if k:
                updated[k] += acc[k - 1]
        updated[0] += Fraction(raw)
        acc = updated
    return tuple(acc)


def taylor_from_hasse(hasse: Sequence[Fraction]) -> tuple[Fraction, ...]:
    """Ordinary derivatives: index k is k! times the Hasse coefficient."""
    factorial = 1
    out: list[Fraction] = []
    for k, coeff in enumerate(hasse):
        if k:
            factorial *= k
        out.append(Fraction(coeff) * factorial)
    return tuple(out)


def degree_from_hasse(hasse: Sequence[Fraction]) -> int:
    """Index of the last nonzero Hasse coefficient, or 0 when the jet is 0.

    Valid when the jet order is at least deg H. The leading Hasse coefficient
    of a polynomial is independent of the base point.
    """
    degree = 0
    for index, coeff in enumerate(hasse):
        if coeff != 0:
            degree = index
    return degree


def hasse_jet_value(
    hasse: Sequence[Fraction],
    base: Fraction,
    t: Fraction,
) -> Fraction:
    """Evaluate Σ hasse[k] (t − base)^k. Exact when the jet reaches deg H."""
    step = Fraction(t) - Fraction(base)
    acc = Fraction(0)
    for coeff in reversed(hasse):
        acc = acc * step + coeff
    return acc


def hermite_multiplicities(degree_cap: int, orders_per_point: int) -> tuple[int, ...]:
    """Multiplicities at t = 1, 2, … that sum to ``degree_cap + 1``.

    Each entry is the number of Hasse coefficients taken at that point
    (indices 0 .. m−1). Every entry equals ``orders_per_point`` except the
    last, which may be shorter. ``orders_per_point`` ≥ 1.
    """
    if degree_cap < 0:
        raise ValueError("degree cap must be nonnegative")
    if orders_per_point < 1:
        raise ValueError("orders_per_point must be at least 1")
    remaining = degree_cap + 1
    multiplicities: list[int] = []
    while remaining:
        take = orders_per_point if orders_per_point < remaining else remaining
        multiplicities.append(take)
        remaining -= take
    return tuple(multiplicities)


def hermite_schemes(degree_cap: int) -> tuple[tuple[str, int], ...]:
    """Distinct (scheme, orders_per_point) patterns for a degree cap.

    ``paired`` keeps at most two coefficients per point (value and first
    Hasse derivative). ``triple`` keeps at most three. ``single-point``
    puts all D+1 coefficients at t = 1. A pattern that repeats an earlier
    multiplicity tuple is omitted; when D+1 ≤ 2 the paired pattern is the
    single-point pattern.
    """
    if degree_cap < 0:
        raise ValueError("degree cap must be nonnegative")
    full = degree_cap + 1
    candidates = (
        ("paired", 2),
        ("triple", 3),
        ("single-point", full),
    )
    seen: list[tuple[int, ...]] = []
    schemes: list[tuple[str, int]] = []
    for name, orders in candidates:
        multiplicities = hermite_multiplicities(degree_cap, orders)
        if multiplicities in seen:
            continue
        seen.append(multiplicities)
        if orders >= full:
            schemes.append(("single-point", full))
        else:
            schemes.append((name, orders))
    return tuple(schemes)


def hermite_eval(
    samples: Sequence[tuple[Fraction, Sequence[Fraction]]],
    t: Fraction,
) -> Fraction:
    """Value at ``t`` of the Hermite interpolant through Hasse data.

    Each sample is ``(node, hasse)`` with ``hasse[k] = H^{[k]}(node)``.
    The length of ``hasse`` is the multiplicity. Nodes must be distinct;
    repetition is expressed by a longer coefficient tuple, not by a second
    sample at the same node.
    """
    if not samples:
        raise ValueError("need at least one Hermite node")
    nodes: list[Fraction] = []
    jets: list[tuple[Fraction, ...]] = []
    sample_of: list[int] = []
    seen: set[Fraction] = set()
    for index, (raw_node, raw_coeffs) in enumerate(samples):
        node = Fraction(raw_node)
        coeffs = tuple(Fraction(c) for c in raw_coeffs)
        if not coeffs:
            raise ValueError("each node needs a value")
        if node in seen:
            raise ValueError("sample points must be distinct; raise multiplicity instead")
        seen.add(node)
        nodes.extend([node] * len(coeffs))
        jets.append(coeffs)
        sample_of.extend([index] * len(coeffs))

    width = len(nodes)
    # Q[i][j] = f[z_i, z_{i-1}, ..., z_{i-j}], so the diagonal is the Newton row.
    triangle: list[list[Fraction]] = [
        [Fraction(0) for _ in range(i + 1)] for i in range(width)
    ]
    for i in range(width):
        triangle[i][0] = jets[sample_of[i]][0]
    for i in range(1, width):
        for j in range(1, i + 1):
            left_node = nodes[i - j]
            right_node = nodes[i]
            if left_node == right_node:
                if sample_of[i] != sample_of[i - j]:
                    raise ValueError("equal nodes are not one multiplicity group")
                jet = jets[sample_of[i]]
                if j >= len(jet):
                    raise ValueError("missing Hasse coefficient for a repeated node")
                triangle[i][j] = jet[j]
            else:
                span = right_node - left_node
                triangle[i][j] = (triangle[i][j - 1] - triangle[i - 1][j - 1]) / span

    acc = triangle[width - 1][width - 1]
    point = Fraction(t)
    for k in range(width - 2, -1, -1):
        acc = triangle[k][k] + (point - nodes[k]) * acc
    return acc


def _one_jet(order: int) -> list[Fraction]:
    jet = [Fraction(0)] * (order + 1)
    jet[0] = Fraction(1)
    return jet


def _sub_from_one(jet: Sequence[Fraction]) -> list[Fraction]:
    out = [-Fraction(c) for c in jet]
    out[0] += 1
    return out


def _jet_mul(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    if len(left) != len(right):
        raise ValueError("jet lengths differ")
    width = len(left)
    out = [Fraction(0)] * width
    for i, left_coeff in enumerate(left):
        if left_coeff == 0:
            continue
        room = width - i
        for j in range(room):
            right_coeff = right[j]
            if right_coeff == 0:
                continue
            out[i + j] += left_coeff * right_coeff
    return out


def _linear_hasse(
    family: ReplacementFamily,
    clause_index: int,
    occurrence: int,
    bits: Sequence[int],
    t: Fraction,
    order: int,
) -> tuple[Fraction, ...]:
    acc: list[Fraction] | None = None
    for variable, coeffs in enumerate(family.rows[clause_index][occurrence]):
        if not bits[variable]:
            continue
        if acc is None:
            acc = [Fraction(c) for c in coeffs]
            continue
        if len(coeffs) > len(acc):
            acc.extend(Fraction(0) for _ in range(len(coeffs) - len(acc)))
        for power, coeff in enumerate(coeffs):
            acc[power] += coeff
    if acc is None:
        acc = [Fraction(0)]
    return polynomial_hasse_jet(acc, t, order)


def eval_H_hasse_jet(
    cnf: CNF,
    family: ReplacementFamily,
    t: Fraction,
    order: int,
) -> tuple[Fraction, ...]:
    """Hasse coefficients H^{[0]}(t), …, H^{[order]}(t).

    One pass over {0,1}^n. For each assignment the clause product is
    differentiated in Q[u]/(u^{order+1}). The pass returns the whole jet;
    the orders are not separate enumerations. Arithmetic per assignment
    grows quadratically with ``order`` because jet products are truncated
    convolutions.
    """
    if order < 0:
        raise ValueError("order must be nonnegative")
    if family.n != cnf.n:
        raise ValueError("replacement family has a different n")
    point = Fraction(t)
    total = [Fraction(0)] * (order + 1)
    for bits in itertools.product((0, 1), repeat=cnf.n):
        clause_product = _one_jet(order)
        killed = False
        for clause_index, clause in enumerate(cnf.clauses):
            if not clause:
                killed = True
                break
            factor: Sequence[Fraction] = _one_jet(order)
            for occurrence in range(len(clause)):
                form = _linear_hasse(
                    family, clause_index, occurrence, bits, point, order
                )
                if family.occ_neg[clause_index][occurrence]:
                    piece: Sequence[Fraction] = form
                else:
                    piece = _sub_from_one(form)
                factor = _jet_mul(factor, piece)
            phi = _sub_from_one(factor)
            clause_product = _jet_mul(clause_product, phi)
        if killed:
            continue
        for index, coeff in enumerate(clause_product):
            total[index] += coeff
    return tuple(total)


def _check_family(cnf: CNF, family: ReplacementFamily) -> None:
    if family.n != cnf.n:
        raise ValueError("replacement family has a different n")
    if family.R < 0:
        raise ValueError("R must be nonnegative")
    if len(family.rows) != cnf.m or len(family.occ_neg) != cnf.m:
        raise ValueError("replacement family does not match the clauses")
    for clause_index, clause in enumerate(cnf.clauses):
        rows = family.rows[clause_index]
        if len(rows) != len(clause) or len(family.occ_neg[clause_index]) != len(clause):
            raise ValueError(f"clause {clause_index} has a mismatched occurrence row")
        for row in rows:
            if len(row) != cnf.n:
                raise ValueError("a coefficient row does not have length n")
            for coeffs in row:
                if len(coeffs) > family.R + 1:
                    raise ValueError("a coefficient polynomial exceeds degree R")


@dataclass(frozen=True)
class HermiteRecovery:
    """One exact Hermite recovery of H(0) against a degree cap.

    ``conditions`` equals D+1, the same count of pieces the Lagrange path
    uses. ``distinct_points`` is strictly smaller when some multiplicity is
    at least 2 and D ≥ 1. ``cube_passes`` is the standalone number of passes
    over {0,1}^n for this scheme: one per distinct point. A cache that reuses
    a jet across schemes does not change this count.
    """

    name: str
    scheme: str
    n: int
    m: int
    R: int
    k: int
    degree_cap: int
    brute: int
    h0: Fraction
    direct_at_zero: Fraction
    lagrange_h0: Fraction
    nodes: tuple[Fraction, ...]
    multiplicities: tuple[int, ...]
    orders_per_point: int
    conditions: int
    lagrange_points: int
    cube_passes: int
    assignments: int
    observed_degree: int
    max_derivative_order: int

    @property
    def distinct_points(self) -> int:
        return len(self.nodes)

    @property
    def fewer_than_lagrange(self) -> bool:
        return self.distinct_points < self.lagrange_points

    def sample_comparison(self) -> str:
        count = self.distinct_points
        noun = "point" if count == 1 else "points"
        mults = ", ".join(str(m) for m in self.multiplicities)
        return (
            f"{count} {noun} × multiplicities ({mults}) "
            f"= {self.conditions} conditions, against D+1 = {self.lagrange_points}"
        )


def recovery_problems(row: HermiteRecovery) -> list[str]:
    """Reasons a Step 4 row failed. Empty means the row met the check."""
    problems: list[str] = []
    if row.h0 != Fraction(row.brute):
        problems.append(f"H(0)={row.h0} != #F={row.brute}")
    if row.h0.denominator != 1:
        problems.append(f"H(0)={row.h0} is not an integer")
    if row.lagrange_h0 != row.h0:
        problems.append(f"Lagrange H(0)={row.lagrange_h0} != Hermite H(0)={row.h0}")
    if row.direct_at_zero != row.h0:
        problems.append(f"direct H(0)={row.direct_at_zero} != Hermite H(0)={row.h0}")
    if row.conditions != row.degree_cap + 1:
        problems.append(
            f"conditions {row.conditions} != D+1 = {row.degree_cap + 1}"
        )
    if row.lagrange_points != row.degree_cap + 1:
        problems.append("Lagrange point count is not D+1")
    if sum(row.multiplicities) != row.conditions:
        problems.append("multiplicities do not sum to the condition count")
    if not row.fewer_than_lagrange:
        problems.append(
            f"{row.distinct_points} distinct points are not fewer than "
            f"D+1 = {row.lagrange_points}"
        )
    if row.observed_degree > row.degree_cap:
        problems.append(
            f"observed degree {row.observed_degree} exceeds D = {row.degree_cap}"
        )
    if row.cube_passes != row.distinct_points:
        problems.append("cube pass count drifted from the distinct-point count")
    if row.assignments != row.distinct_points * (1 << row.n):
        problems.append("assignment count is not distinct points times 2^n")
    if row.max_derivative_order < 1:
        problems.append("no derivative was used")
    expected_nodes = tuple(Fraction(i) for i in range(1, row.distinct_points + 1))
    if row.nodes != expected_nodes:
        problems.append(f"sample points {row.nodes} are not 1..{row.distinct_points}")
    if row.orders_per_point < 2:
        problems.append("orders_per_point is below 2")
    return problems


def _jet_from_cache(
    cnf: CNF,
    family: ReplacementFamily,
    t: Fraction,
    order: int,
    cache: dict[Fraction, tuple[Fraction, ...]] | None,
) -> tuple[Fraction, ...]:
    if cache is not None:
        cached = cache.get(t)
        if cached is not None and len(cached) >= order + 1:
            return cached[: order + 1]
    fresh = eval_H_hasse_jet(cnf, family, t, order)
    if cache is not None:
        cached = cache.get(t)
        if cached is None or len(fresh) > len(cached):
            cache[t] = fresh
    return fresh


def recover_H0_hermite(
    cnf: CNF,
    family: ReplacementFamily,
    *,
    k: int,
    orders_per_point: int,
    name: str = "",
    scheme: str = "",
    lagrange_h0: Fraction | None = None,
    direct_at_zero: Fraction | None = None,
    observed_degree: int | None = None,
    jet_cache: dict[Fraction, tuple[Fraction, ...]] | None = None,
) -> HermiteRecovery:
    """Recover H(0) from Hasse data at t = 1, 2, … against the degree cap.

    The cap is ``degree_bound(cnf, family.R, k)``. The scheme takes
    ``orders_per_point`` Hasse coefficients at each node until D+1
    conditions are in hand. ``lagrange_h0``, when omitted, is computed by
    ``interpolate_H0`` on this same family and the same k.
    """
    if k < 0:
        raise ValueError("k must be nonnegative")
    _check_family(cnf, family)
    degree_cap = degree_bound(cnf, R=family.R, k=k)
    multiplicities = hermite_multiplicities(degree_cap, orders_per_point)
    samples: list[tuple[Fraction, tuple[Fraction, ...]]] = []
    for index, multiplicity in enumerate(multiplicities):
        node = Fraction(index + 1)
        jet = _jet_from_cache(cnf, family, node, multiplicity - 1, jet_cache)
        if len(jet) != multiplicity:
            raise RuntimeError("Hasse jet length drifted from the multiplicity")
        samples.append((node, jet))
    h0 = hermite_eval(samples, Fraction(0))

    if observed_degree is None:
        full = _jet_from_cache(cnf, family, Fraction(1), degree_cap, jet_cache)
        if len(full) != degree_cap + 1:
            raise RuntimeError("full Hasse jet does not have length D+1")
        observed_degree = degree_from_hasse(full)
    if direct_at_zero is None:
        direct_at_zero = eval_H_at(cnf, family, Fraction(0))
    if lagrange_h0 is None:
        lagrange_h0, lagrange_cap, _used = interpolate_H0(cnf, fam=family, k=k)
        if lagrange_cap != degree_cap:
            raise RuntimeError("Lagrange degree cap drifted")

    max_order = max(multiplicities) - 1
    points = len(multiplicities)
    return HermiteRecovery(
        name=name,
        scheme=scheme or f"coefficients-{orders_per_point}",
        n=cnf.n,
        m=cnf.m,
        R=family.R,
        k=k,
        degree_cap=degree_cap,
        brute=brute_count(cnf),
        h0=h0,
        direct_at_zero=direct_at_zero,
        lagrange_h0=lagrange_h0,
        nodes=tuple(node for node, _jet in samples),
        multiplicities=multiplicities,
        orders_per_point=orders_per_point,
        conditions=sum(multiplicities),
        lagrange_points=degree_cap + 1,
        cube_passes=points,
        assignments=points * (1 << cnf.n),
        observed_degree=observed_degree,
        max_derivative_order=max_order,
    )


@dataclass(frozen=True)
class Step4Case:
    """One formula and replacement family checked by Step 4."""

    name: str
    cnf: CNF
    family: ReplacementFamily
    k: int
    note: str


def _family(
    n: int,
    R: int,
    rows: list[list[list[list[Fraction]]]],
    occ_vars: list[list[int]],
    occ_neg: list[list[bool]],
) -> ReplacementFamily:
    return ReplacementFamily(
        n=n,
        R=R,
        rows=rows,
        occ_vars=occ_vars,
        occ_neg=occ_neg,
    )


def negative_unit_case() -> Step4Case:
    """One negative unit, ⟨P, X⟩ = X(1+T). Closed form H(T) = 1 − T, #F = 1."""
    cnf = CNF(1, [[(0, True)]])
    family = _family(
        1,
        1,
        [[[[Fraction(1), Fraction(1)]]]],
        [[0]],
        [[True]],
    )
    return Step4Case(
        name="negative-unit",
        cnf=cnf,
        family=family,
        k=1,
        note="H(T) = 1 - T. Width 1, so the cap is m·1·R = 1.",
    )


def unit_pair_case() -> Step4Case:
    """Two positive units, ⟨P, X⟩ = X_v + T(X_0 + X_1).

    Closed form H(T) = 1 + 6T + 6T^2, #F = 1.
    """
    cnf = CNF(2, [[(0, False)], [(1, False)]])

    def row(variable: int) -> list[list[Fraction]]:
        return [
            [Fraction(1 if i == variable else 0), Fraction(1)]
            for i in range(2)
        ]

    family = _family(
        2,
        1,
        [[row(0)], [row(1)]],
        [[0], [1]],
        [[False], [False]],
    )
    return Step4Case(
        name="unit-pair",
        cnf=cnf,
        family=family,
        k=1,
        note="H(T) = 1 + 6T + 6T^2. Width 1, so the cap is m·1·R = 2.",
    )


def contradiction_case() -> Step4Case:
    """(x) ∧ (¬x) with ⟨P, X⟩ = X(1+T) on both occurrences.

    Closed form H(T) = −T − T^2. H(0) = #F = 0, and H is not the zero
    polynomial, so the recovered value is not a trivial constant jet.
    """
    cnf = CNF(1, [[(0, False)], [(0, True)]])

    def occurrence() -> list[list[Fraction]]:
        return [[Fraction(1), Fraction(1)]]

    family = _family(
        1,
        1,
        [[occurrence()], [occurrence()]],
        [[0], [0]],
        [[False], [True]],
    )
    return Step4Case(
        name="contradiction-n1",
        cnf=cnf,
        family=family,
        k=1,
        note="H(T) = -T - T^2. #F = 0 and deg H = 2.",
    )


def _random_case(
    name: str,
    n: int,
    m: int,
    R: int,
    cnf_seed: int,
    replacement_seed: int,
) -> Step4Case:
    cnf = random_3cnf(n, m, rng=random.Random(cnf_seed))
    family = random_replacements(cnf, R=R, rng=random.Random(replacement_seed))
    return Step4Case(
        name=name,
        cnf=cnf,
        family=family,
        k=3,
        note=(
            f"random 3-CNF, formula seed {cnf_seed}, replacement seed {replacement_seed}, "
            f"cap m·3·R."
        ),
    )


def _bench_bounded_case() -> Step4Case:
    document = load_bench_document(bench_path())
    record = next(
        item for item in document["instances"] if item["name"] == "n6-m4-r1"
    )
    cnf = cnf_from_record(record)
    family = random_replacements(cnf, R=1, rng=random.Random(42))
    return Step4Case(
        name="bounded-n6-m4-r1",
        cnf=cnf,
        family=family,
        k=3,
        note=(
            "Step 1 bench n6-m4-r1, formula seed 42, replacement seed 42. "
            "Cap m·3·R = 12."
        ),
    )


def step4_cases() -> tuple[Step4Case, ...]:
    """Hand closed forms, Step 0 random formulas, and one Step 1 bench row.

    Every 3-CNF uses k = 3. The hand formulas are width 1 and use k = 1,
    which is the product m·(clause width)·R for those clauses.
    """
    one_clause = CNF(3, [[(0, False), (1, False), (2, False)]])
    one_clause_family = random_replacements(
        one_clause, R=1, rng=random.Random(0)
    )
    return (
        negative_unit_case(),
        unit_pair_case(),
        contradiction_case(),
        Step4Case(
            name="one-clause-n3",
            cnf=one_clause,
            family=one_clause_family,
            k=3,
            note="Positive 3-clause on n = 3, replacement seed 0, cap m·3·R = 3.",
        ),
        _random_case("random-n4-m4-r1", 4, 4, 1, 11, 111),
        _random_case("random-n4-m3-r2", 4, 3, 2, 13, 113),
        _random_case("random-n6-m6-r1", 6, 6, 1, 14, 114),
        _bench_bounded_case(),
    )


def analyze_step4_case(case: Step4Case) -> tuple[HermiteRecovery, ...]:
    """Run every distinct Hermite scheme on one case.

    One full jet at t = 1 measures the degree and feeds the single-point
    scheme. Lagrange's H(0) and the direct sum at 0 are computed once and
    shared. Reported cube passes stay equal to the scheme's own point count.
    """
    _check_family(case.cnf, case.family)
    degree_cap = degree_bound(case.cnf, R=case.family.R, k=case.k)
    cache: dict[Fraction, tuple[Fraction, ...]] = {
        Fraction(1): eval_H_hasse_jet(case.cnf, case.family, Fraction(1), degree_cap)
    }
    observed = degree_from_hasse(cache[Fraction(1)])
    direct = eval_H_at(case.cnf, case.family, Fraction(0))
    lagrange_h0, lagrange_cap, _family_used = interpolate_H0(
        case.cnf, fam=case.family, k=case.k
    )
    if lagrange_cap != degree_cap:
        raise RuntimeError("Lagrange degree cap drifted")
    rows: list[HermiteRecovery] = []
    for scheme, orders in hermite_schemes(degree_cap):
        rows.append(
            recover_H0_hermite(
                case.cnf,
                case.family,
                k=case.k,
                orders_per_point=orders,
                name=case.name,
                scheme=scheme,
                lagrange_h0=lagrange_h0,
                direct_at_zero=direct,
                observed_degree=observed,
                jet_cache=cache,
            )
        )
    return tuple(rows)
