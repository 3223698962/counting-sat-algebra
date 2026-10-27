"""Section 4.5 placement using only occupied coordinates and sparse rows.

No dense auxiliary matrix or occurrence-by-variable array is allocated.
Resource bounds are checked before measurement modes and replacement support
are enumerated. The schedule entry point checks every Lemma 4.7 numerical
hypothesis, then runs the exact conditional-expectation construction.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass, fields
from fractions import Fraction
from typing import Mapping, Sequence

from .cnf import CNF
from .conditional_balance import ConditionalSchedules, choose_schedules
from .linear_algebra import SparseEliminationBudgetExceeded, SparseRationalRowSpan
from .reference_parameters import ReferenceParameters, check_reference_parameters
from .replacements import eval_poly
from .scheduler import scopes_from_cnf
from .structured_auxiliary import Label, StructuredAuxiliary, build_structured_auxiliary

SparsePolynomialRow = Mapping[int, tuple[Fraction, ...]]
SparseMeasurementRow = Mapping[int, Fraction]


class SparsePlacementBudgetExceeded(ValueError):
    """A planned sparse construction exceeded a caller's resource limit."""


@dataclass(frozen=True)
class SparsePlacementBudgets:
    max_groups: int = 256
    max_alphabet: int = 32
    max_trials: int = 300_000
    max_modes: int = 4096
    max_mode_coordinates: int = 2_000_000
    max_pullback_candidates: int = 2_000_000
    max_coefficient_entries: int = 2_000_000
    max_measurement_entries: int = 2_000_000
    max_elimination_entries: int = 2_000_000

    def __post_init__(self) -> None:
        for field in fields(self):
            value = operator.index(getattr(self, field.name))
            if value < 0:
                raise ValueError("sparse budgets must be nonnegative integers")
            object.__setattr__(self, field.name, value)


@dataclass(frozen=True)
class SparsePlacementStatistics:
    occupied_coordinates: int
    mode_coordinate_queries: int
    pullback_candidates: int
    polynomial_entries: int
    measurement_entries: int
    dense_polynomial_slots: int


@dataclass(frozen=True)
class SparsePlacedFormula:
    cnf: CNF
    auxiliary: StructuredAuxiliary
    variable_schedules: tuple[tuple[int, ...], ...]
    clause_schedules: tuple[tuple[int, ...], ...]
    rows: tuple[tuple[SparsePolynomialRow, ...], ...]
    measurement_rows: tuple[SparseMeasurementRow, ...]
    copies: int
    occurrence_mass: int
    max_width: int
    statistics: SparsePlacementStatistics

    @property
    def cap_occurrences(self) -> int:
        return self.occurrence_mass * (self.auxiliary.R - 1)

    @property
    def cap_step1(self) -> int:
        return self.cnf.m * self.max_width * self.auxiliary.R


@dataclass(frozen=True)
class SparseDesignReport:
    problems: tuple[str, ...]
    group_ranks: tuple[tuple[int, ...], ...] = ()
    max_group_size: int = 0
    elimination_entries: int = 0
    budget_refused: bool = False

    @property
    def ok(self) -> bool:
        return not self.problems and not self.budget_refused


@dataclass(frozen=True)
class SparseScheduleOutcome:
    status: str
    reason: str
    parameters: ReferenceParameters
    placed: SparsePlacedFormula | None = None
    schedules: ConditionalSchedules | None = None
    design: SparseDesignReport | None = None

    @property
    def accepted(self) -> bool:
        return (self.status == "placed" and self.parameters.admitted
                and self.placed is not None and self.design is not None and self.design.ok)


def _size(value: int) -> str:
    return str(value) if value.bit_length() <= 256 else f"an integer of {value.bit_length()} bits"


def _limit(value: int, budget: int, name: str) -> None:
    if value > budget:
        raise SparsePlacementBudgetExceeded(f"{name} bound {_size(value)} exceeds budget {budget}; construction not started")


def _schedule(schedule: Sequence[int], auxiliary: StructuredAuxiliary) -> tuple[int, ...]:
    parsed = tuple(operator.index(symbol) for symbol in schedule)
    if len(parsed) != auxiliary.Q or any(symbol < 0 or symbol >= auxiliary.r for symbol in parsed):
        raise ValueError("schedules must have length Q and symbols in {0, ..., r-1}")
    return parsed


def _add_poly(row: dict[int, tuple[Fraction, ...]], var: int, poly: tuple[Fraction, ...]) -> None:
    previous = row.get(var, ())
    combined = [Fraction(0)] * max(len(previous), len(poly))
    for i, value in enumerate(previous):
        combined[i] += value
    for i, value in enumerate(poly):
        combined[i] += value
    while combined and not combined[-1]:
        combined.pop()
    if combined:
        row[var] = tuple(combined)
    else:
        row.pop(var, None)


def place_sparse_formula(
    cnf: CNF, auxiliary: StructuredAuxiliary,
    variable_schedules: Sequence[Sequence[int]], clause_schedules: Sequence[Sequence[int]], *,
    budgets: SparsePlacementBudgets = SparsePlacementBudgets(),
) -> SparsePlacedFormula:
    """Pull back sparse P and C rows along distinct occurrence coordinates.

    This accepts supplied schedules without asserting the balance hypotheses.
    ``schedule_sparse_placement`` is the gate for admitted reference schedules.
    """
    if len(variable_schedules) != cnf.n or len(clause_schedules) != cnf.m:
        raise ValueError("expected one schedule per declared variable and per clause")
    variables = tuple(_schedule(schedule, auxiliary) for schedule in variable_schedules)
    clauses = tuple(_schedule(schedule, auxiliary) for schedule in clause_schedules)
    _limit(auxiliary.r, budgets.max_groups, "group count")
    modes_count = auxiliary.measurement_mode_count
    _limit(modes_count, budgets.max_modes, "measurement mode")

    counts: dict[Label, int] = {}
    placements: list[dict[Label, int]] = []
    occurrences = []
    for e, clause in enumerate(cnf.clauses):
        for occ, (var, _neg) in enumerate(clause):
            label = variables[var], clauses[e]
            copy = counts.get(label, 0)
            counts[label] = copy + 1
            if copy == len(placements):
                placements.append({})
            placements[copy][label] = var
            occurrences.append((e, occ, var, label, copy))
    mass = len(occurrences)
    mode_queries = mass * modes_count
    pullback_bound = sum(min(len(placements[copy]), auxiliary.replacement_support_bound(label))
                         for _e, _occ, _var, label, copy in occurrences)
    _limit(mode_queries, budgets.max_mode_coordinates, "mode-coordinate query")
    _limit(mode_queries, budgets.max_measurement_entries, "measurement entry")
    _limit(pullback_bound, budgets.max_pullback_candidates, "pullback candidate")
    _limit(pullback_bound * auxiliary.R, budgets.max_coefficient_entries, "polynomial coefficient")

    modes = tuple(auxiliary.measurement_modes(max_modes=budgets.max_modes))
    measurement_rows = []
    for placement in placements:
        for mode in modes:
            pulled: dict[int, Fraction] = {}
            for label, var in placement.items():
                value = auxiliary.measurement_entry(mode, label)
                if value:
                    pulled[var] = pulled.get(var, Fraction(0)) + value
            pulled = {var: value for var, value in pulled.items() if value}
            if pulled:
                measurement_rows.append(pulled)

    rows = [[{} for _ in clause] for clause in cnf.clauses]
    support_cache = {}
    candidates = 0
    for e, occ, _var, label, copy in occurrences:
        placement = placements[copy]
        support_bound = auxiliary.replacement_support_bound(label)
        if support_bound <= len(placement):
            if label not in support_cache:
                support_cache[label] = tuple(auxiliary.replacement_support(label, max_terms=support_bound))
            for column, poly in support_cache[label]:
                candidates += 1
                var = placement.get(column)
                if var is not None:
                    _add_poly(rows[e][occ], var, poly)
        else:
            for column, var in placement.items():
                candidates += 1
                poly = auxiliary.replacement_entry(label, column)
                if any(poly):
                    _add_poly(rows[e][occ], var, poly)

    polynomial_entries = sum(len(poly) for clause in rows for row in clause for poly in row.values())
    measurement_entries = sum(map(len, measurement_rows))
    if candidates > pullback_bound or polynomial_entries > pullback_bound * auxiliary.R or measurement_entries > mode_queries:
        raise ArithmeticError("sparse placement exceeded its planned allocation bound")
    statistics = SparsePlacementStatistics(mass, mode_queries, candidates, polynomial_entries,
                                          measurement_entries, mass * cnf.n * auxiliary.R)
    return SparsePlacedFormula(cnf, auxiliary, variables, clauses, tuple(tuple(clause) for clause in rows),
                               tuple(measurement_rows), len(placements), mass,
                               max((len(clause) for clause in cnf.clauses), default=0), statistics)


def sparse_groups_at(placed: SparsePlacedFormula, t: int) -> tuple[tuple[int, ...], ...]:
    t = operator.index(t)
    if not 1 <= t <= placed.auxiliary.Q:
        raise ValueError("t must lie inside {1, ..., Q}")
    groups = [[] for _ in range(placed.auxiliary.r)]
    for var, schedule in enumerate(placed.variable_schedules):
        groups[schedule[t - 1]].append(var)
    return tuple(tuple(group) for group in groups)


def check_sparse_design(
    placed: SparsePlacedFormula, *, tau: Fraction | None = None, max_entries: int = 2_000_000,
) -> SparseDesignReport:
    """Check restoration, denominators, partitions and rational kernel conditions.

    Restricted matrices and targets use global sparse coordinates. Independent
    original integer rows provide the actual per-group measurement count.
    """
    max_entries = operator.index(max_entries)
    if max_entries < 0:
        raise ValueError("max_entries must be nonnegative")
    if tau is not None:
        tau = Fraction(tau)
        if not 0 < tau <= Fraction(1, 2):
            raise ValueError("need 0 < tau <= 1/2")
    n, auxiliary = placed.cnf.n, placed.auxiliary
    problems = []
    try:
        if len(placed.variable_schedules) != n or len(placed.clause_schedules) != placed.cnf.m:
            raise ValueError("schedule counts do not match the formula")
        for schedule in placed.variable_schedules + placed.clause_schedules:
            _schedule(schedule, auxiliary)
    except (ValueError, TypeError) as exc:
        return SparseDesignReport((str(exc),))
    if placed.occurrence_mass != sum(map(len, placed.cnf.clauses)):
        return SparseDesignReport(("occurrence mass does not match the formula",))
    if len(placed.rows) != placed.cnf.m:
        return SparseDesignReport(("replacement clause count does not match the formula",))
    denominator = math.factorial(auxiliary.R - 1) * math.factorial(auxiliary.Q)
    coefficient_bound = denominator * (1 + (1 + 2 ** auxiliary.Q) * auxiliary.R * (4 * auxiliary.Q) ** (auxiliary.R - 1))
    for e, clause in enumerate(placed.cnf.clauses):
        if len(placed.rows[e]) != len(clause):
            problems.append(f"clause {e}: replacement occurrence count does not match")
            continue
        for occ, (var, _neg) in enumerate(clause):
            row = placed.rows[e][occ]
            if any(v < 0 or v >= n or not poly for v, poly in row.items()):
                problems.append(f"clause {e} occurrence {occ}: invalid sparse polynomial coordinate")
                continue
            if row.get(var, (Fraction(0),))[0] != 1 or any(poly[0] for v, poly in row.items() if v != var):
                problems.append(f"clause {e} occurrence {occ}: T=0 does not restore X_{var}")
            if any(len(poly) > auxiliary.R for poly in row.values()):
                problems.append(f"clause {e} occurrence {occ}: degree exceeds R-1")
            scaled = [denominator * value for poly in row.values() for value in poly]
            if any(value.denominator != 1 for value in scaled):
                problems.append(f"clause {e} occurrence {occ}: coefficient denominator exceeds the common denominator")
            if sum(map(abs, scaled), Fraction(0)) > coefficient_bound:
                problems.append(f"clause {e} occurrence {occ}: coefficient norm exceeds Lemma 4.6's bound")
    for i, row in enumerate(placed.measurement_rows):
        if any(v < 0 or v >= n for v in row):
            problems.append(f"measurement {i}: coordinate outside the declared variables")
        if any(value.denominator != 1 for value in row.values()):
            problems.append(f"measurement {i}: coefficients must be integral")
        if sum(map(abs, row.values()), Fraction(0)) > auxiliary.dimension:
            problems.append(f"measurement {i}: absolute row sum exceeds D")
    if problems:
        return SparseDesignReport(tuple(problems))

    ranks = []
    max_group = 0
    peak_entries = 0
    try:
        for t in range(1, auxiliary.Q + 1):
            groups = sparse_groups_at(placed, t)
            owner = [0] * n
            for g, variables in enumerate(groups):
                max_group = max(max_group, len(variables))
                if len(variables) * auxiliary.r > 2 * n:
                    problems.append(f"t={t}, group {g}: group size exceeds 2n/r")
                for var in variables:
                    owner[var] = g
            restricted = [[] for _ in groups]
            for row in placed.measurement_rows:
                pieces = {}
                for var, value in row.items():
                    pieces.setdefault(owner[var], {})[var] = value
                for g, piece in pieces.items():
                    restricted[g].append(piece)
            spans = []
            entries = 0
            for matrix in restricted:
                span = SparseRationalRowSpan(matrix, n, max_entries=max_entries - entries)
                spans.append(span)
                entries += span.stored_entries
            peak_entries = max(peak_entries, entries)
            node_ranks = tuple(span.rank for span in spans)
            ranks.append(node_ranks)
            if tau is not None and sum(node_ranks) > 2 * tau * n:
                problems.append(f"t={t}: independent measurement count exceeds 2*tau*n")
            for e, clause in enumerate(placed.cnf.clauses):
                assigned = placed.clause_schedules[e][t - 1]
                for occ, _literal in enumerate(clause):
                    targets = {}
                    for var, poly in placed.rows[e][occ].items():
                        group = owner[var]
                        if group == assigned:
                            continue
                        value = eval_poly(list(poly), Fraction(t))
                        if value:
                            targets.setdefault(group, {})[var] = value
                    for group, target in targets.items():
                        if not spans[group].contains(target):
                            problems.append(f"t={t}, clause {e}, occurrence {occ}: outside group {group} is not in the measurement span")
    except SparseEliminationBudgetExceeded as exc:
        return SparseDesignReport(tuple(problems) + (str(exc),), tuple(ranks), max_group, peak_entries, True)
    return SparseDesignReport(tuple(problems), tuple(ranks), max_group, peak_entries)


def schedule_sparse_placement(
    cnf: CNF, *, r: int, Q: int, R: int, tau: Fraction = Fraction(1, 2),
    A0: int | None = None, budgets: SparsePlacementBudgets = SparsePlacementBudgets(),
) -> SparseScheduleOutcome:
    """Admit every Lemma 4.7 hypothesis, then choose and verify sparse schedules."""
    scopes = scopes_from_cnf(cnf)
    parameters = check_reference_parameters(cnf.n, scopes, r=r, Q=Q, R=R, tau=tau, A0=A0)
    if not parameters.admitted:
        return SparseScheduleOutcome("refused", "; ".join(parameters.reasons), parameters)
    schedules = None
    try:
        alphabet = r ** Q
        _limit(r, budgets.max_groups, "group count")
        _limit(alphabet, budgets.max_alphabet, "schedule alphabet")
        _limit((cnf.n + cnf.m) * alphabet, budgets.max_trials, "candidate trial")
        auxiliary = build_structured_auxiliary(r, Q, R, parameters.interval[0])
        _limit(auxiliary.measurement_mode_count, budgets.max_modes, "measurement mode")
        schedules = choose_schedules(cnf.n, scopes, r, Q)
        placed = place_sparse_formula(cnf, auxiliary, schedules.variables, schedules.factors, budgets=budgets)
    except SparsePlacementBudgetExceeded as exc:
        return SparseScheduleOutcome("refused", str(exc), parameters, schedules=schedules)
    design = check_sparse_design(placed, tau=tau, max_entries=budgets.max_elimination_entries)
    if design.budget_refused:
        return SparseScheduleOutcome("refused", design.problems[-1], parameters, placed, schedules, design)
    if not design.ok or placed.copies > Fraction(2 * parameters.balance.M, auxiliary.dimension):
        return SparseScheduleOutcome("invalid-design", "sparse placement failed its exact design or multiplicity check", parameters, placed, schedules, design)
    return SparseScheduleOutcome("placed", "all Lemma 4.7 hypotheses, sparse kernel and restoration checks passed", parameters, placed, schedules, design)
