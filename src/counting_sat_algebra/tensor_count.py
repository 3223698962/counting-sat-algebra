"""Budgeted variable elimination of the actual replaced factors over Q[z].

Tables contain Taylor jets at T=1, through the full component degree bound.
Eliminating a Boolean variable preserves the sum-product polynomial. This
finite contraction can handle connected supports with small intermediate
tables. It is not Appendix D.4's recursive modular contraction or its bound.
"""

from __future__ import annotations

import heapq
import operator
from dataclasses import dataclass, fields
from fractions import Fraction

from .hermite import _jet_mul, _one_jet, _sub_from_one, hasse_jet_value, polynomial_hasse_jet
from .replacement_matrix import PlacedFormula
from .sparse_placement import SparsePlacedFormula
from .support_count import _placed_design_error, _support_components


@dataclass(frozen=True)
class ContractionBudgets:
    max_variables: int = 1_000_000
    max_component_degree: int = 64
    max_bucket_variables: int = 12
    max_boolean_vectors: int = 500_000
    max_jet_products: int = 2_000_000
    max_table_entries: int = 500_000
    max_table_coefficients: int = 2_000_000
    max_planning_steps: int = 2_000_000

    def __post_init__(self):
        for field in fields(self):
            value = operator.index(getattr(self, field.name))
            if value < 0:
                raise ValueError('contraction budgets must be nonnegative integers')
            object.__setattr__(self, field.name, value)


@dataclass(frozen=True)
class ContractionPlan:
    components: int
    free_variables: int
    largest_component: int
    max_component_degree: int
    max_bucket_variables: int
    boolean_vector_bound: int
    jet_product_bound: int
    peak_table_entries: int
    peak_table_coefficients: int
    planning_steps: int


@dataclass(frozen=True)
class ContractionResult:
    status: str
    count: int | None
    reason: str
    plan: ContractionPlan | None = None
    components_counted: int = 0
    boolean_vectors: int = 0
    factor_evaluations: int = 0
    jet_products: int = 0
    peak_table_entries: int = 0
    peak_table_coefficients: int = 0

    @property
    def counted(self) -> bool:
        return self.status == 'counted'


@dataclass(frozen=True)
class _Elimination:
    inputs: tuple[int, ...]
    output: int
    scope: tuple[int, ...]  # The eliminated variable is first.


@dataclass(frozen=True)
class _ComponentPlan:
    scopes: tuple[tuple[int, ...], ...]
    steps: tuple[_Elimination, ...]


@dataclass(frozen=True)
class _Tensor:
    scope: tuple[int, ...]
    values: tuple[tuple[Fraction, ...], ...]


class _BudgetExceeded(ValueError):
    pass


def _limit(value, limit, name):
    if value > limit:
        size = str(value) if value.bit_length() <= 256 else f'an integer of {value.bit_length()} bits'
        raise _BudgetExceeded(f'{name} bound {size} exceeds budget {limit}; no numerical tables created')


def _plan_components(components, free, budgets):
    """Plan all components before constructing any numerical factor table."""
    totals = dict(vectors=0, products=0, planning=0, width=0, entries=0, coefficients=0)

    def charge(key, amount, limit):
        totals[key] += amount
        _limit(totals[key], limit, key)

    def memory(entries, degree):
        totals['entries'] = max(totals['entries'], entries)
        totals['coefficients'] = max(totals['coefficients'], entries * (degree + 1))
        _limit(totals['entries'], budgets.max_table_entries, 'live tensor entries')
        _limit(totals['coefficients'], budgets.max_table_coefficients, 'live tensor coefficients')

    plans = []
    for component in components:
        _limit(component.degree, budgets.max_component_degree, 'component degree')
        adjacent = {v: set() for v in range(component.variables)}
        incident = {v: set() for v in adjacent}
        initial_scopes = []
        scopes, sizes = {}, {}
        live = 0
        for fid, clause in enumerate(component.clauses):
            scope = tuple(sorted({v for _neg, row in clause for v, _poly in row}))
            _limit(len(scope), budgets.max_bucket_variables, 'factor support')
            totals['width'] = max(totals['width'], len(scope))
            entries = 1 << len(scope)
            charge('vectors', entries, budgets.max_boolean_vectors)
            charge('products', entries * len(clause), budgets.max_jet_products)
            charge('planning', len(scope) ** 2 + len(scope) + 1, budgets.max_planning_steps)
            live += entries
            memory(live, component.degree)
            initial_scopes.append(scope)
            scopes[fid], sizes[fid] = scope, entries
            for v in scope:
                adjacent[v].update(u for u in scope if u != v)
                incident[v].add(fid)

        queue = [(len(neighbors), v) for v, neighbors in adjacent.items()]
        heapq.heapify(queue)
        steps = []
        for offset in range(component.variables):
            while True:
                degree, v = heapq.heappop(queue)
                charge('planning', 1, budgets.max_planning_steps)
                if v in adjacent and degree == len(adjacent[v]):
                    break
            neighbors = tuple(sorted(adjacent[v]))
            scope = (v, *neighbors)
            _limit(len(scope), budgets.max_bucket_variables, 'elimination bucket')
            totals['width'] = max(totals['width'], len(scope))
            inputs = tuple(sorted(incident[v]))
            output = len(component.clauses) + offset
            entries = 1 << len(neighbors)
            charge('vectors', 2 * entries, budgets.max_boolean_vectors)
            charge('products', 2 * entries * len(inputs), budgets.max_jet_products)
            charge('planning', len(neighbors) ** 2 + len(neighbors) + 1
                   + sum(len(scopes[fid]) for fid in inputs), budgets.max_planning_steps)
            # The output is allocated while the input tables are still live.
            memory(live + entries, component.degree)
            steps.append(_Elimination(inputs, output, scope))
            for fid in inputs:
                for u in scopes.pop(fid):
                    incident[u].remove(fid)
                live -= sizes.pop(fid)
            del incident[v]
            del adjacent[v]
            for u in neighbors:
                adjacent[u].discard(v)
                adjacent[u].update(w for w in neighbors if w != u)
                incident[u].add(output)
                heapq.heappush(queue, (len(adjacent[u]), u))
            scopes[output], sizes[output] = neighbors, entries
            live += entries
        plans.append(_ComponentPlan(tuple(initial_scopes), tuple(steps)))

    public = ContractionPlan(
        len(components), free, max((c.variables for c in components), default=0),
        max((c.degree for c in components), default=0), totals['width'],
        totals['vectors'], totals['products'], totals['entries'], totals['coefficients'], totals['planning'],
    )
    return public, tuple(plans)


def _factor_tensor(clause, scope, degree, stats):
    owner = {v: i for i, v in enumerate(scope)}
    source = tuple(tuple((owner[v], polynomial_hasse_jet(poly, Fraction(1), degree))
                         for v, poly in row) for _neg, row in clause)
    values = []
    for assignment in range(1 << len(scope)):
        stats['boolean_vectors'] += 1
        stats['factor_evaluations'] += 1
        complement = _one_jet(degree)
        for (neg, _row), jets in zip(clause, source):
            linear = [Fraction(0)] * (degree + 1)
            for position, jet in jets:
                if assignment & (1 << position):
                    for j, coefficient in enumerate(jet):
                        linear[j] += coefficient
            complement = _jet_mul(complement, linear if neg else _sub_from_one(linear))
            stats['jet_products'] += 1
        values.append(tuple(_sub_from_one(complement)))
    return _Tensor(scope, tuple(values))


def _contract_component(component, plan, stats):
    tables = {}
    live = 0

    def memory(entries):
        stats['peak_table_entries'] = max(stats['peak_table_entries'], entries)
        stats['peak_table_coefficients'] = max(stats['peak_table_coefficients'], entries * (component.degree + 1))

    for fid, (clause, scope) in enumerate(zip(component.clauses, plan.scopes)):
        tables[fid] = _factor_tensor(clause, scope, component.degree, stats)
        live += len(tables[fid].values)
        memory(live)
    for step in plan.steps:
        output_scope = step.scope[1:]
        zero = (Fraction(0),) * (component.degree + 1)
        values = [zero] * (1 << len(output_scope))
        memory(live + len(values))
        inputs = tuple((tables[fid], tuple(step.scope.index(v) for v in tables[fid].scope))
                       for fid in step.inputs)
        for assignment in range(1 << len(step.scope)):
            stats['boolean_vectors'] += 1
            product = _one_jet(component.degree)
            for tensor, positions in inputs:
                index = sum(((assignment >> position) & 1) << j for j, position in enumerate(positions))
                product = _jet_mul(product, tensor.values[index])
                stats['jet_products'] += 1
            index = assignment >> 1
            values[index] = tuple(left + right for left, right in zip(values[index], product))
        # Release references before removing input tables and allocating the
        # next output; the live-table accounting relies on these lifetimes.
        del inputs, tensor
        for fid in step.inputs:
            live -= len(tables.pop(fid).values)
        tables[step.output] = _Tensor(output_scope, tuple(values))
        live += len(values)
        del values
    if len(tables) != 1:
        raise ArithmeticError('connected support did not contract to one scalar tensor')
    scalar = next(iter(tables.values()))
    if scalar.scope or len(scalar.values) != 1:
        raise ArithmeticError('contraction left an unsummed variable')
    return scalar.values[0]


def count_placed_formula_by_contraction(
    placed: PlacedFormula | SparsePlacedFormula, *, budgets: ContractionBudgets = ContractionBudgets(),
) -> ContractionResult:
    """Count a validated placement by deterministic minimum-degree elimination.

    Full Taylor jets retain all coefficients needed for H_i(0). Budgets cover
    symbolic planning, intermediate support width, Boolean vectors, jet
    multiplications and live table slots. Slot counts exclude input polynomial
    storage, planner metadata, arithmetic temporaries and integer bit lengths.
    Every component is planned before any numerical table is constructed;
    a refused call returns no count. This backend does not use memoization.
    """
    if placed.cnf.n > budgets.max_variables:
        return ContractionResult('refused', None, f'declared variables exceed budget {budgets.max_variables}')
    error = _placed_design_error(placed)
    if error:
        return ContractionResult(error[0], None, error[1])
    if any(not clause for clause in placed.cnf.clauses):
        return ContractionResult('counted', 0, 'an empty factor makes H(T) identically zero')
    components, free = _support_components(placed)
    try:
        plan, local_plans = _plan_components(components, free, budgets)
    except _BudgetExceeded as exc:
        return ContractionResult('refused', None, str(exc))
    stats = dict(components_counted=0, boolean_vectors=0, factor_evaluations=0,
                 jet_products=0, peak_table_entries=0, peak_table_coefficients=0)
    count = 2 ** free
    for component, local_plan in zip(components, local_plans):
        jet = _contract_component(component, local_plan, stats)
        local_count = hasse_jet_value(jet, Fraction(1), Fraction(0))
        if local_count.denominator != 1 or not 0 <= local_count <= 2 ** component.variables:
            raise ArithmeticError('recovered contraction count is not an integer in [0, 2**size]')
        count *= int(local_count)
        stats['components_counted'] += 1
    if (stats['boolean_vectors'] != plan.boolean_vector_bound
            or stats['jet_products'] != plan.jet_product_bound
            or stats['peak_table_entries'] != plan.peak_table_entries
            or stats['peak_table_coefficients'] != plan.peak_table_coefficients):
        raise ArithmeticError('contraction disagreed with its symbolic resource plan')
    return ContractionResult('counted', count, 'variable elimination and complete local Taylor recovery finished', plan, **stats)
