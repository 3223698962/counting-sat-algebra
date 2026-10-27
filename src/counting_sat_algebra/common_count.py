"""Section 3's common counting routine for an explicitly placed CNF.

Lemma 3.2 streams integer measurement tuples and two scans of each local
Boolean cube. Lemma 3.3 recovers Taylor jets through fixed substitutions U
and restrictions E, then Hermite interpolation gives H(0). Neither the
global assignment cube nor an attained-measurement table is stored.

The operation budget is checked before any enumeration. This implements
the common routine on supplied designs; the rest of Theorem 1.1's pipeline
and its asymptotic parameter choices are separate obligations.

``separated_group_sum`` is Lemma 3.2 for an arbitrary coefficient ring.
The CNF counter uses it over ``Q[z]/(z^J)``. The Section 5 replacement sum
uses the same scan over ``A_ell[z]/(z^J)``.
"""

from __future__ import annotations

import itertools
import math
import operator
from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping, Sequence

from .hermite import _jet_mul, _one_jet, _sub_from_one, hermite_eval, polynomial_hasse_jet
from .linear_algebra import SparseRationalRowSpan
from .replacement_matrix import PlacedFormula, design_problems_at, groups_at, restoration_problems
from .sparse_placement import SparsePlacedFormula, check_sparse_design, sparse_groups_at

COMMON_VECTOR_BUDGET = 500_000
COMMON_PAIR_BUDGET = 20_000
COMMON_JET_BUDGET = 64


@dataclass(frozen=True)
class CommonCountPlan:
    degree_cap: int
    q: int
    J: int
    substitution_pairs: int
    measurement_counts: tuple[int, ...]
    measurement_tuples: tuple[int, ...]
    local_cube_sizes: tuple[tuple[int, ...], ...]
    boolean_vector_bound: int


@dataclass(frozen=True)
class CommonCountResult:
    status: str
    count: int | None
    reason: str
    plan: CommonCountPlan
    jets: tuple[tuple[Fraction, ...], ...] = ()
    boolean_vectors: int = 0
    measurement_tuples: int = 0
    factor_evaluations: int = 0

    @property
    def counted(self) -> bool:
        return self.status == "counted"


@dataclass(frozen=True)
class _Node:
    groups: tuple[tuple[int, ...], ...]
    rows: tuple[tuple[Mapping[int, Fraction], ...], ...]
    ranges: tuple[range, ...]
    counts: tuple[int, ...]


def _node(placed: PlacedFormula | SparsePlacedFormula, t: int) -> _Node:
    sparse = isinstance(placed, SparsePlacedFormula)
    groups = sparse_groups_at(placed, t) if sparse else groups_at(placed, t)
    matrices = []
    ranges = []
    for variables in groups:
        positions = {v: i for i, v in enumerate(variables)}
        if sparse:
            restricted = [{positions[v]: value for v, value in row.items() if v in positions}
                          for row in placed.measurement_rows]
        else:
            restricted = [{i: row[v] for i, v in enumerate(variables) if row[v]}
                          for row in placed.measurement_rows]
        # Keep original integer rows. Row operations are used only to select an
        # independent subset, so the absolute-sum bound and integrality survive.
        rows = SparseRationalRowSpan(restricted, len(variables)).independent_rows
        matrices.append(rows)
        for row in rows:
            lo = int(sum(value for value in row.values() if value < 0))
            hi = int(sum(value for value in row.values() if value > 0))
            ranges.append(range(lo, hi + 1))
    return _Node(groups, tuple(matrices), tuple(ranges), tuple(map(len, matrices)))


def _plan(placed: PlacedFormula | SparsePlacedFormula, nodes: Sequence[_Node]) -> CommonCountPlan:
    degree = placed.cap_occurrences
    q = degree // placed.auxiliary.Q
    n = placed.cnf.n
    if q >= n:
        pairs = 3 ** n
    else:
        pairs = term = 1
        for size in range(1, q + 1):
            term = term * (n - size + 1) * 2 // size
            pairs += term
    tuples = tuple(math.prod(len(values) for values in node.ranges) for node in nodes)
    local = tuple(tuple(2 ** len(group) for group in node.groups) for node in nodes)
    bound = pairs * sum(2 * count * sum(cubes) for count, cubes in zip(tuples, local))
    return CommonCountPlan(degree, q, q + 1, pairs, tuple(sum(node.counts) for node in nodes), tuples, local, bound)


def _substitutions(n: int, q: int):
    for size in range(min(n, q) + 1):
        for E in itertools.combinations(range(n), size):
            for u_size in range(size + 1):
                for U in itertools.combinations(E, u_size):
                    yield E, U, (-1) ** (size - u_size)


def _measure(rows, bits):
    return tuple(sum((value * bits[i] for i, value in row.items()), Fraction(0)) for row in rows)


def separated_group_sum(n, groups, rows, ranges, counts, fixed, evaluate, zero, one, add, multiply, stats):
    """Lemma 3.2. ``evaluate(group, bits)`` returns one element of the coefficient ring.

    ``add`` and ``multiply`` must return new values and leave their arguments
    unchanged. ``zero`` and ``one`` are reused, so they must not be mutated.
    The two scans of each local cube are sequential. No assignment table and
    no attained-measurement table is stored. ``stats[0]`` counts Boolean
    vectors and ``stats[1]`` counts measurement tuples.
    """
    total = zero
    for measurement in itertools.product(*ranges):
        stats[1] += 1
        representatives = []
        cursor = 0
        targets = []
        for group, matrix, count in zip(groups, rows, counts):
            target = measurement[cursor:cursor + count]
            cursor += count
            targets.append(target)
            representative = None
            for bits in itertools.product((0, 1), repeat=len(group)):
                stats[0] += 1
                if _measure(matrix, bits) == target:
                    representative = bits
                    break
            if representative is None:
                break
            representatives.append(representative)
        if len(representatives) != len(groups):
            continue
        baseline = [0] * n
        for group, representative in zip(groups, representatives):
            for variable, bit in zip(group, representative):
                baseline[variable] = bit
        product = one
        for g, (group, matrix, target) in enumerate(zip(groups, rows, targets)):
            subtotal = zero
            for local in itertools.product((0, 1), repeat=len(group)):
                stats[0] += 1
                if any(variable in fixed and bit == 0 for variable, bit in zip(group, local)):
                    continue
                if _measure(matrix, local) != target:
                    continue
                bits = list(baseline)
                for variable, bit in zip(group, local):
                    bits[variable] = bit
                subtotal = add(subtotal, evaluate(g, bits))
            product = multiply(product, subtotal)
        total = add(total, product)
    return total


def _add_jets(left, right):
    return tuple(left_coeff + right_coeff for left_coeff, right_coeff in zip(left, right))


def _mul_jets(left, right):
    return tuple(_jet_mul(left, right))


def _restricted_sum(placed, node, t, E, U, coeff_jets, q, stats):
    """Lemma 3.2 in Q[z]/(z**(q+1)), with X fixed to the indicator of U."""
    n = placed.cnf.n
    fixed = set(E)
    corrections = []
    for clause_jets in coeff_jets:
        clause = []
        for occurrence in clause_jets:
            jet = [Fraction(0)] * (q + 1)
            for v in U:
                for j in range(1, q + 1):
                    if v in occurrence:
                        jet[j] += occurrence[v][j]
            clause.append(tuple(jet))
        corrections.append(tuple(clause))
    assigned = [[] for _ in node.groups]
    for e, schedule in enumerate(placed.clause_schedules):
        assigned[schedule[t - 1]].append(e)

    def evaluate(g, bits):
        product = _one_jet(q)
        for e in assigned[g]:
            stats[2] += 1
            complement = _one_jet(q)
            for occ, (_v, neg) in enumerate(placed.cnf.clauses[e]):
                linear = list(corrections[e][occ])
                linear[0] = sum((coeff[0] * bits[v] for v, coeff in coeff_jets[e][occ].items()), Fraction(0))
                complement = _jet_mul(complement, linear if neg else _sub_from_one(linear))
            product = _jet_mul(product, _sub_from_one(complement))
            if not any(product):
                break
        return product

    zero = (Fraction(0),) * (q + 1)
    return separated_group_sum(
        n, node.groups, node.rows, node.ranges, node.counts, fixed, evaluate,
        zero, tuple(_one_jet(q)), _add_jets, _mul_jets, stats,
    )


def count_placed_formula(
    placed: PlacedFormula | SparsePlacedFormula, *, vector_budget: int = COMMON_VECTOR_BUDGET,
    pair_budget: int = COMMON_PAIR_BUDGET, jet_budget: int = COMMON_JET_BUDGET,
) -> CommonCountResult:
    """Recover the exact count using Sections 3.3–3.4 on a supplied placement.

    Budgets limit local Boolean vectors, (E,U) substitutions and jet length.
    Refusal returns no count and starts no enumeration. This routine checks
    restoration and the kernel condition, but makes no claim that the design
    parameters yield a subexponential bound on a family of general formulas.
    """
    vector_budget = operator.index(vector_budget)
    pair_budget = operator.index(pair_budget)
    jet_budget = operator.index(jet_budget)
    if min(vector_budget, pair_budget, jet_budget) < 0:
        raise ValueError("budgets must be nonnegative integers")
    if isinstance(placed, SparsePlacedFormula):
        design = check_sparse_design(placed)
        problems = list(design.problems)
        if design.budget_refused:
            empty = CommonCountPlan(placed.cap_occurrences, 0, 0, 0, (), (), (), 0)
            return CommonCountResult("refused", None, problems[-1], empty)
    else:
        problems = list(restoration_problems(placed))
        for t in range(1, placed.auxiliary.Q + 1):
            problems.extend(design_problems_at(placed, t))
    if problems:
        # Do not build measurement ranges from malformed/nonintegral rows.
        empty = CommonCountPlan(placed.cap_occurrences, 0, 0, 0, (), (), (), 0)
        return CommonCountResult("invalid-design", None, problems[0], empty)
    nodes = tuple(_node(placed, t) for t in range(1, placed.auxiliary.Q + 1))
    plan = _plan(placed, nodes)
    for value, limit, name in (
        (plan.boolean_vector_bound, vector_budget, "Boolean vector"),
        (plan.substitution_pairs, pair_budget, "substitution pair"),
        (plan.J, jet_budget, "jet length"),
    ):
        if value > limit:
            size = str(value) if value.bit_length() <= 256 else f"an integer of {value.bit_length()} bits"
            return CommonCountResult("refused", None, f"{name} bound {size} exceeds budget {limit}; no enumeration started", plan)
    stats = [0, 0, 0]
    jets = []
    for t, node in enumerate(nodes, start=1):
        source = placed.rows if isinstance(placed, SparsePlacedFormula) else placed.family.rows
        coeff_jets = []
        for clause in source:
            clause_jets = []
            for row in clause:
                items = row.items() if isinstance(placed, SparsePlacedFormula) else enumerate(row)
                clause_jets.append({v: polynomial_hasse_jet(poly, Fraction(t), plan.q)
                                    for v, poly in items if any(poly)})
            coeff_jets.append(clause_jets)
        jet = [Fraction(0)] * plan.J
        for E, U, sign in _substitutions(placed.cnf.n, plan.q):
            value = _restricted_sum(placed, node, t, E, U, coeff_jets, plan.q, stats)
            for j, coefficient in enumerate(value):
                jet[j] += sign * coefficient
        jets.append(tuple(jet))
    h0 = hermite_eval([(Fraction(t), jet) for t, jet in enumerate(jets, 1)], Fraction(0))
    if h0.denominator != 1 or h0 < 0 or h0 > 2 ** placed.cnf.n:
        raise ArithmeticError("restored count is not an integer in [0, 2**n]")
    if stats[0] > plan.boolean_vector_bound:
        raise ArithmeticError("local enumeration exceeded its planned vector bound")
    return CommonCountResult("counted", int(h0), "Section 3 shared-measurement jets and Hermite recovery completed", plan, tuple(jets), *stats)
