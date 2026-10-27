"""Factor the actual replaced polynomial by its variable-support components.

For disjoint support sets V_i, H(T) = 2**free * product_i H_i(T).
Each H_i(1+z) is summed over its small local Boolean cube, with a full
Taylor jet sufficient to recover H_i(0). This avoids the global measurement
box and global fixed substitutions when the replaced factors disconnect.
Large connected support sets still refuse; this is a finite contraction
method, not the general Appendix D algorithm or its complexity bound.
"""

from __future__ import annotations

import itertools
import operator
from dataclasses import dataclass
from fractions import Fraction

from .hermite import _jet_mul, _one_jet, _sub_from_one, hermite_eval, polynomial_hasse_jet
from .replacement_matrix import PlacedFormula, design_problems_at, restoration_problems
from .sparse_placement import SparsePlacedFormula, check_sparse_design


@dataclass(frozen=True)
class SupportCountPlan:
    components: int
    free_variables: int
    largest_component: int
    degree_cap: int
    max_component_degree: int
    boolean_vector_bound: int


@dataclass(frozen=True)
class SupportCountResult:
    status: str
    count: int | None
    reason: str
    plan: SupportCountPlan | None = None
    components_counted: int = 0
    local_evaluations: int = 0
    memo_hits: int = 0
    boolean_vectors: int = 0
    factor_evaluations: int = 0

    @property
    def counted(self) -> bool:
        return self.status == "counted"


def _component_jet(clauses, n, degree, stats):
    """Contract one small support component in Q[z]/(z**(degree+1))."""
    source = tuple(tuple((neg, tuple((v, polynomial_hasse_jet(poly, Fraction(1), degree))
                                    for v, poly in row)) for neg, row in clause)
                   for clause in clauses)
    total = [Fraction(0)] * (degree + 1)
    for bits in itertools.product((0, 1), repeat=n):
        stats['boolean_vectors'] += 1
        product = _one_jet(degree)
        for clause in source:
            stats['factor_evaluations'] += 1
            complement = _one_jet(degree)
            for neg, row in clause:
                linear = [Fraction(0)] * (degree + 1)
                for v, jet in row:
                    if bits[v]:
                        for j, coefficient in enumerate(jet):
                            linear[j] += coefficient
                complement = _jet_mul(complement, linear if neg else _sub_from_one(linear))
            product = _jet_mul(product, _sub_from_one(complement))
            if not any(product):
                break
        for j, coefficient in enumerate(product):
            total[j] += coefficient
    return tuple(total)


def _placed_design_error(placed):
    """Shared validation for evaluators of the actual replacement factors."""
    sparse = isinstance(placed, SparsePlacedFormula)
    if sparse:
        report = check_sparse_design(placed)
        if report.budget_refused:
            return "refused", report.problems[-1]
        problems = report.problems
    else:
        problems = list(restoration_problems(placed))
        for t in range(1, placed.auxiliary.Q + 1):
            problems.extend(design_problems_at(placed, t))
    if problems:
        return "invalid-design", problems[0]
    return None


@dataclass(frozen=True)
class _SupportComponent:
    variables: int
    clauses: tuple
    degree: int


def _support_components(placed):
    """Extract locally numbered factors after restoration has been validated."""
    n = placed.cnf.n
    sparse = isinstance(placed, SparsePlacedFormula)
    source = placed.rows if sparse else placed.family.rows
    parent = list(range(n))
    sizes = [1] * n
    active = [False] * n

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    factors, scopes, degrees = [], [], []
    for clause, rows in zip(placed.cnf.clauses, source):
        factor = []
        scope = set()
        degree = 0
        for (_original_v, neg), row in zip(clause, rows):
            items = row.items() if sparse else enumerate(row)
            entries = tuple(sorted((v, tuple(Fraction(value) for value in poly))
                                   for v, poly in items if any(poly)))
            scope.update(v for v, _poly in entries)
            degree += max((len(poly) - 1 for _v, poly in entries), default=0)
            factor.append((neg, entries))
        scope = tuple(sorted(scope))
        for v in scope:
            active[v] = True
            left, right = find(scope[0]), find(v)
            if left != right:
                if sizes[left] < sizes[right]:
                    left, right = right, left
                parent[right] = left
                sizes[left] += sizes[right]
        factors.append(tuple(factor))
        scopes.append(scope)
        degrees.append(degree)

    variables = {}
    for v, used in enumerate(active):
        if used:
            variables.setdefault(find(v), []).append(v)
    grouped = {root: [] for root in variables}
    component_degrees = {root: 0 for root in variables}
    for factor, scope, degree in zip(factors, scopes, degrees):
        root = find(scope[0])
        grouped[root].append(factor)
        component_degrees[root] += degree
    components = []
    for root, group in variables.items():
        owner = {v: i for i, v in enumerate(group)}
        clauses = tuple(tuple((neg, tuple((owner[v], poly) for v, poly in row))
                              for neg, row in factor) for factor in grouped[root])
        components.append(_SupportComponent(len(group), clauses, component_degrees[root]))
    return tuple(components), n - sum(active)


def count_placed_formula_by_support(
    placed: PlacedFormula | SparsePlacedFormula, *, max_variables: int = 1_000_000,
    max_component_variables: int = 12, max_component_degree: int = 64,
    vector_budget: int = 500_000, memo_entries: int = 4096,
) -> SupportCountResult:
    """Recover a count from the same placement using disjoint factor supports.

    Components are formed from every nonzero replacement coefficient in each
    clause, not from the original CNF. Thus cross-variable replacements cannot
    be incorrectly separated. The vector plan assumes no memoization; optional
    capped reuse can only reduce the actual work. Refusals return no partial
    count and start no local cube enumeration.
    """
    max_variables, max_component_variables, max_component_degree, vector_budget, memo_entries = (
        operator.index(value) for value in
        (max_variables, max_component_variables, max_component_degree, vector_budget, memo_entries)
    )
    if min(max_variables, max_component_variables, max_component_degree, vector_budget, memo_entries) < 0:
        raise ValueError("support budgets must be nonnegative integers")
    if placed.cnf.n > max_variables:
        return SupportCountResult("refused", None, f"declared variables exceed budget {max_variables}")
    error = _placed_design_error(placed)
    if error:
        return SupportCountResult(error[0], None, error[1])
    if any(not clause for clause in placed.cnf.clauses):
        return SupportCountResult("counted", 0, "an empty factor makes H(T) identically zero")

    components, free = _support_components(placed)
    plan = SupportCountPlan(
        len(components), free, max((c.variables for c in components), default=0),
        sum(c.degree for c in components), max((c.degree for c in components), default=0),
        sum(2 ** c.variables for c in components),
    )
    for value, limit, name in (
        (plan.largest_component, max_component_variables, "support component size"),
        (plan.max_component_degree, max_component_degree, "support component degree"),
        (plan.boolean_vector_bound, vector_budget, "local Boolean vector"),
    ):
        if value > limit:
            size = str(value) if value.bit_length() <= 256 else f"an integer of {value.bit_length()} bits"
            return SupportCountResult("refused", None,
                                      f"{name} bound {size} exceeds budget {limit}; no enumeration started", plan)

    stats = dict(components_counted=0, local_evaluations=0, memo_hits=0,
                 boolean_vectors=0, factor_evaluations=0)
    memo = {}
    total = 2 ** plan.free_variables
    for component in components:
        key = component.variables, component.clauses
        if key in memo:
            local_count = memo[key]
            stats['memo_hits'] += 1
        else:
            jet = _component_jet(component.clauses, component.variables, component.degree, stats)
            h0 = hermite_eval([(Fraction(1), jet)], Fraction(0))
            if h0.denominator != 1 or not 0 <= h0 <= 2 ** component.variables:
                raise ArithmeticError("recovered support count is not an integer in [0, 2**size]")
            local_count = int(h0)
            stats['local_evaluations'] += 1
            if len(memo) < memo_entries:
                memo[key] = local_count
        total *= local_count
        stats['components_counted'] += 1
        if not total:
            break
    if stats['boolean_vectors'] > plan.boolean_vector_bound:
        raise ArithmeticError("support contraction exceeded its planned vector bound")
    return SupportCountResult("counted", total,
                              "disjoint replacement-support sums and local Taylor recovery completed", plan, **stats)
