"""Exact CNF counting by independent components and sparse algebraic recovery.

Tautologies and duplicate clauses are removed without changing the models.
Declared unused variables contribute a factor of two each. Every nontrivial
local count is obtained from the shared-measurement routine, with optional
memoization of identical local formulas. No brute-force oracle is a backend.

This preprocessing can handle large disconnected formulas. A large connected
component is still subject to the algebraic evaluator's resource limits; the
reference's general-input asymptotic bound is a separate obligation.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass

from .cnf import CNF
from .common_count import COMMON_VECTOR_BUDGET, count_placed_formula
from .conditional_balance import choose_schedules
from .sparse_placement import (
    SparsePlacementBudgetExceeded, SparsePlacementBudgets, place_sparse_formula,
)
from .structured_auxiliary import StructuredAuxiliary, build_structured_auxiliary


@dataclass(frozen=True)
class ComponentCountResult:
    status: str
    count: int | None
    reason: str
    components: int = 0
    components_counted: int = 0
    free_variables: int = 0
    largest_component: int = 0
    normalized_clauses: int = 0
    local_evaluations: int = 0
    memo_hits: int = 0
    boolean_vectors: int = 0
    measurement_tuples: int = 0
    factor_evaluations: int = 0

    @property
    def counted(self) -> bool:
        return self.status == "counted"


def count_cnf_components(
    cnf: CNF, *, auxiliary: StructuredAuxiliary | None = None,
    max_variables: int = 1_000_000, max_component_variables: int = 20,
    vector_budget: int = COMMON_VECTOR_BUDGET, total_vector_budget: int = 5_000_000,
    memo_entries: int = 4096, placement_budgets: SparsePlacementBudgets = SparsePlacementBudgets(),
) -> ComponentCountResult:
    """Count components exactly, with an explicit refusal and no partial count.

    Memoization is optional and capped; identical renumbered formulas reuse
    only a completed exact count. Components are placed afresh, since their
    H(t) need not factor from a placement of the global formula.
    """
    max_variables, max_component_variables, vector_budget, total_vector_budget, memo_entries = (
        operator.index(value) for value in
        (max_variables, max_component_variables, vector_budget, total_vector_budget, memo_entries)
    )
    if min(max_variables, max_component_variables, vector_budget, total_vector_budget, memo_entries) < 0:
        raise ValueError("component budgets must be nonnegative integers")
    if cnf.n > max_variables:
        return ComponentCountResult("refused", None, f"declared variables exceed budget {max_variables}")

    normalized = set()
    for clause in cnf.clauses:
        literals = set(clause)
        if any((var, not neg) in literals for var, neg in literals):
            continue
        if not literals:
            return ComponentCountResult("counted", 0, "an empty clause makes the formula unsatisfiable")
        normalized.add(tuple(sorted(literals)))
    clauses = tuple(sorted(normalized))
    if not clauses:
        return ComponentCountResult("counted", 2 ** cnf.n, "all declared variables are free", free_variables=cnf.n)

    parent = list(range(cnf.n))
    sizes = [1] * cnf.n
    active = [False] * cnf.n

    def find(var):
        while parent[var] != var:
            parent[var] = parent[parent[var]]
            var = parent[var]
        return var

    for clause in clauses:
        anchor = clause[0][0]
        for var, _neg in clause:
            active[var] = True
            left, right = find(anchor), find(var)
            if left != right:
                if sizes[left] < sizes[right]:
                    left, right = right, left
                parent[right] = left
                sizes[left] += sizes[right]
    variables = {}
    for var, used in enumerate(active):
        if used:
            variables.setdefault(find(var), []).append(var)
    grouped_clauses = {root: [] for root in variables}
    for clause in clauses:
        grouped_clauses[find(clause[0][0])].append(clause)
    components = sorted(variables, key=lambda root: variables[root][0])
    free = cnf.n - sum(active)
    largest = max(map(len, variables.values()))
    stats = dict(components=len(components), components_counted=0, free_variables=free,
                 largest_component=largest, normalized_clauses=len(clauses), local_evaluations=0,
                 memo_hits=0, boolean_vectors=0, measurement_tuples=0, factor_evaluations=0)
    if largest > max_component_variables:
        return ComponentCountResult("refused", None,
                                    f"largest component has {largest} variables, above budget {max_component_variables}; no local enumeration started", **stats)
    if auxiliary is None:
        auxiliary = build_structured_auxiliary(2, 2, 1)
    if auxiliary.r != 2:
        raise ValueError("the local component placement requires r=2 for arbitrary small components")

    memo = {}
    total = 2 ** free
    for root in components:
        owner = {var: local for local, var in enumerate(variables[root])}
        local_clauses = tuple(tuple((owner[var], neg) for var, neg in clause) for clause in grouped_clauses[root])
        key = len(owner), local_clauses
        if key in memo:
            local_count = memo[key]
            stats['memo_hits'] += 1
        else:
            local = CNF(len(owner), [list(clause) for clause in local_clauses])
            alphabet = auxiliary.r ** auxiliary.Q
            if alphabet > placement_budgets.max_alphabet or (local.n + local.m) * alphabet > placement_budgets.max_trials:
                return ComponentCountResult("refused", None, "local schedule search exceeds its alphabet or trial budget", **stats)
            schedules = choose_schedules(local.n, tuple(tuple(var for var, _neg in clause) for clause in local_clauses),
                                         auxiliary.r, auxiliary.Q)
            try:
                placed = place_sparse_formula(local, auxiliary, schedules.variables, schedules.factors, budgets=placement_budgets)
            except SparsePlacementBudgetExceeded as exc:
                return ComponentCountResult("refused", None, str(exc), **stats)
            remaining = total_vector_budget - stats['boolean_vectors']
            result = count_placed_formula(placed, vector_budget=min(vector_budget, remaining))
            stats['local_evaluations'] += 1
            stats['boolean_vectors'] += result.boolean_vectors
            stats['measurement_tuples'] += result.measurement_tuples
            stats['factor_evaluations'] += result.factor_evaluations
            if not result.counted:
                return ComponentCountResult(result.status, None, result.reason, **stats)
            local_count = result.count
            if len(memo) < memo_entries:
                memo[key] = local_count
        stats['components_counted'] += 1
        total *= local_count
        if not total:
            return ComponentCountResult("counted", 0, "a completed local count is zero", **stats)
    return ComponentCountResult("counted", total, "independent component counts and declared free-bit factor completed", **stats)
