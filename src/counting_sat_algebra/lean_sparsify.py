"""Stream the disjoint counting frontier described in SqrtLogSpace.

SparsePriority selects a smallest width and then a largest literal core.
DisjointBranch replaces *all* clauses of that width containing the core.
SparseCountProgram discards completed leaves and keeps one sibling per depth.
The common, tensor and native leaves are existing finite algebraic evaluators;
the hash leaf assembles the scheduled sparse code, runs the violation-mass
greedy row selection and counts the separate final accepted query. Their
general input time and space bounds have not been transferred from Lean to
Python.
"""

from __future__ import annotations

import operator
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, fields, replace
from itertools import combinations
from math import comb

from .cnf import CNF, Literal
from .component_count import count_cnf_components
from .lean_encoding import clause_bit_length, formula_bit_length, literal_code
from .lean_greedy import GreedyBudgets, count_cnf_greedy
from .lean_parameters import choose_sparse_parameter, occurrence_cap
from .lean_query import count_cnf_ternary
from .lean_ring import LeanAlgebraBudgets
from .sparse_placement import SparsePlacementBudgets
from .tensor_count import ContractionBudgets

ClauseSet = frozenset[Literal]
Formula = tuple[ClauseSet, ...]
STEP23_SCOPE = (
    'Lean-guided canonical syntax and streaming disjoint sparsification with finite algebraic '
    'leaves, including a native ternary squarefree query backend. The SqrtLogSpace code '
    'search, packing, streamed Taylor leaf algorithm and complete Python bit-time and '
    'space analysis remain open. No Theorem 1.1, no 2^{o(n)} claim.'
)


class StreamingRefusal(ValueError):
    """An explicit limit prevents completion; no partial total is published."""


@dataclass(frozen=True)
class StreamingBudgets:
    max_variables: int = 1_000_000
    max_width: int = 8
    max_clauses: int = 100_000
    max_formula_bits: int = 16_000_000
    max_frontier_bits: int = 16_000_000
    max_nodes: int = 100_000
    max_depth: int = 10_000
    max_search_candidates: int = 2_000_000
    max_parameter_exponent: int = 128

    def __post_init__(self) -> None:
        for field in fields(self):
            value = operator.index(getattr(self, field.name))
            if value < 0:
                raise ValueError('streaming budgets must be nonnegative integers')
            object.__setattr__(self, field.name, value)


@dataclass(frozen=True)
class PriorityChoice:
    width: int
    core: ClauseSet
    degree: int
    threshold: int


@dataclass(frozen=True)
class StreamingCountResult:
    status: str
    count: int | None
    reason: str
    n: int = 0
    width: int = 0
    M: int = 0
    parameter_test_passed: bool = False
    parameter_candidates: int = 0
    input_bits: int = 0
    normalized_bits: int = 0
    normalized_clauses: int = 0
    nodes: int = 0
    planning_nodes: int = 0
    branches: int = 0
    leaves: int = 0
    leaves_counted: int = 0
    depth: int = 0
    max_pending_nodes: int = 0
    peak_frontier_bits: int = 0
    search_candidates: int = 0
    occurrence_bound: int = 0
    max_leaf_occurrence: int = 0
    accumulator_bits: int = 0
    count_word_bits: int = 0
    boolean_vectors: int = 0
    factor_evaluations: int = 0
    local_evaluations: int = 0
    memo_hits: int = 0
    jet_products: int = 0
    peak_table_entries: int = 0
    peak_table_coefficients: int = 0
    native_queries: int = 0
    native_source_polynomials: int = 0
    monomial_products: int = 0
    scalar_operations: int = 0
    peak_polynomial_terms: int = 0
    peak_live_coefficients: int = 0
    peak_algebra_bits: int = 0
    max_ring_precision: int = 0
    hash_rounds: int = 0
    hash_candidates: int = 0
    hash_marked_queries: int = 0
    hash_accepted_queries: int = 0
    hash_history_rows: int = 0
    hash_initial_mass: int = 0
    hash_final_mass: int = 0
    same_declared_n: bool = True
    completed_leaves_retained: int = 0
    backend: str = 'lean'
    leaf_backend: str = 'common'
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False
    recursive_contraction: bool = False

    @property
    def counted(self) -> bool:
        return self.status == 'counted'


def _clause_key(clause: ClauseSet) -> tuple[int, ...]:
    return tuple(sorted(literal_code(literal) for literal in clause))


def _normalize(clauses: Sequence[ClauseSet]) -> Formula:
    unique = sorted(set(clauses), key=lambda clause: (len(clause), _clause_key(clause)))
    minimal: list[ClauseSet] = []
    by_size: dict[int, set[ClauseSet]] = {}
    for clause in unique:
        if any((variable, not negated) in clause for variable, negated in clause):
            continue
        absorbed = False
        for size, group in by_size.items():
            if size >= len(clause):
                continue
            if len(group) <= comb(len(clause), size):
                absorbed = any(smaller <= clause for smaller in group)
            else:
                absorbed = any(frozenset(subset) in group for subset in combinations(clause, size))
            if absorbed:
                break
        if not absorbed:
            minimal.append(clause)
            by_size.setdefault(len(clause), set()).add(clause)
    return tuple(sorted(minimal, key=_clause_key))


def normalize_lean_cnf(cnf: CNF) -> CNF:
    """Remove tautologies, repetitions and absorbed clauses while retaining n."""
    formula_bit_length(cnf, canonical=False)  # Validate even mutable CNF inputs.
    return _to_cnf(cnf.n, _normalize(tuple(frozenset(clause) for clause in cnf.clauses)))


def _to_cnf(n: int, formula: Formula) -> CNF:
    return CNF(n, [sorted(clause, key=literal_code) for clause in formula])


def _search_size(formula: Formula) -> int:
    return sum((1 << len(clause)) - 2 for clause in formula if len(clause) > 1)


def select_priority(
    clauses: Sequence[ClauseSet], M: int, *, max_candidates: int = 1_000_000,
) -> PriorityChoice | None:
    """SparsePriority.select without tuples of parameter-dependent size.

    Equal width/core-size candidates retain their first canonical occurrence,
    matching the ordered-clause/powerset enumeration used by the Lean model.
    """
    M = operator.index(M)
    limit = operator.index(max_candidates)
    if M < 2 or limit < 0:
        raise ValueError('M must be >= 2 and the candidate budget nonnegative')
    formula = tuple(sorted(set(clauses), key=_clause_key))
    if _search_size(formula) > limit:
        raise StreamingRefusal('core candidate search exceeds its remaining budget')
    degrees: dict[tuple[int, ClauseSet], int] = {}
    for clause in formula:
        literals = sorted(clause, key=literal_code)
        for size in range(1, len(literals)):
            for subset in combinations(literals, size):
                key = len(clause), frozenset(subset)
                degrees[key] = degrees.get(key, 0) + 1
    best: PriorityChoice | None = None
    for (width, core), degree in degrees.items():
        gap = width - len(core)
        # A passing threshold cannot exceed the known finite degree.
        threshold = 1
        for _ in range(gap):
            threshold *= M
            if threshold > degree:
                break
        if threshold > degree:
            continue
        if best is None or (width, -len(core)) < (best.width, -len(best.core)):
            best = PriorityChoice(width, core, degree, threshold)
    return best


def split_priority(formula: Sequence[ClauseSet], choice: PriorityChoice) -> tuple[Formula, Formula]:
    """DisjointBranch.childA/childB, including every selected clause and petal."""
    if not choice.core or len(choice.core) >= choice.width:
        raise ValueError('a branch needs a nonempty proper core')
    selected = frozenset(
        clause for clause in formula if len(clause) == choice.width and choice.core <= clause
    )
    if len(selected) != choice.degree or len(selected) < choice.threshold:
        raise ValueError('branch degree does not describe the current formula')
    others = tuple(clause for clause in formula if clause not in selected)
    petals = tuple(clause - choice.core for clause in selected)
    units = tuple(frozenset(((variable, not negated),)) for variable, negated in choice.core)
    return _normalize((*others, choice.core)), _normalize((*others, *petals, *units))


def _body_bits(formula: Formula) -> int:
    return 1 + sum(clause_bit_length(tuple(clause)) for clause in formula)


def _occurrences(formula: Formula) -> int:
    occurrences: dict[int, int] = {}
    for clause in formula:
        for variable in {literal[0] for literal in clause}:
            occurrences[variable] = occurrences.get(variable, 0) + 1
    return max(occurrences.values(), default=0)


@dataclass
class _Traversal:
    nodes: int = 0
    branches: int = 0
    leaves: int = 0
    depth: int = 0
    max_pending: int = 0
    peak_bits: int = 0
    candidates: int = 0
    max_occurrence: int = 0


def _leaves(
    root: Formula, n: int, M: int, bound: int, limits: StreamingBudgets, stats: _Traversal,
) -> Iterator[Formula]:
    body = _body_bits(root)
    pending = [(root, 0, body)]
    frontier_bits = 1 + 2 * body + 2
    if frontier_bits > limits.max_frontier_bits:
        raise StreamingRefusal('initial serialized frontier exceeds its bit budget')
    stats.max_pending = 1
    stats.peak_bits = frontier_bits
    while pending:
        if stats.nodes >= limits.max_nodes:
            raise StreamingRefusal('sparsification node budget exceeded')
        formula, depth, body = pending.pop()
        frontier_bits -= 2 * body + 2
        stats.nodes += 1
        stats.depth = max(stats.depth, depth)
        candidate_size = _search_size(formula)
        choice = select_priority(
            formula, M, max_candidates=limits.max_search_candidates // 2 - stats.candidates,
        )
        stats.candidates += candidate_size
        if choice is None:
            maximum = _occurrences(formula)
            if maximum > bound:
                raise RuntimeError('terminal occurrence bound was violated')
            stats.max_occurrence = max(stats.max_occurrence, maximum)
            stats.leaves += 1
            yield formula
            del formula
            continue
        if depth >= limits.max_depth:
            raise StreamingRefusal('sparsification depth budget exceeded')
        core_bits = clause_bit_length(tuple(choice.core))
        units_bits = sum(
            clause_bit_length(((variable, not negated),)) for variable, negated in choice.core
        )
        if n + 1 + body + max(core_bits, units_bits) > limits.max_formula_bits:
            raise StreamingRefusal('child formula encoding upper bound exceeds its bit budget')
        upper_frontier = frontier_bits + 4 * body + 2 * core_bits + 2 * units_bits + 4
        if upper_frontier > limits.max_frontier_bits:
            raise StreamingRefusal('child frontier encoding upper bound exceeds its bit budget')
        left, right = split_priority(formula, choice)
        left_bits, right_bits = _body_bits(left), _body_bits(right)
        pending.append((right, depth + 1, right_bits))
        pending.append((left, depth + 1, left_bits))
        frontier_bits += 2 * (left_bits + right_bits) + 4
        stats.branches += 1
        stats.max_pending = max(stats.max_pending, len(pending))
        stats.peak_bits = max(stats.peak_bits, frontier_bits)
        if len(pending) > depth + 2:
            raise RuntimeError('depth-first frontier retained too many siblings')
        del formula, left, right, choice


def count_cnf_streaming(
    cnf: CNF, *, k: int | None = None, M: int | None = None,
    sigma_numerator: int = 1, sigma_denominator: int = 1,
    budgets: StreamingBudgets | None = None, leaf_backend: str = 'common',
    max_component_variables: int = 20, vector_budget: int = 500_000,
    total_vector_budget: int = 5_000_000, memo_entries: int = 0,
    placement_budgets: SparsePlacementBudgets = SparsePlacementBudgets(),
    tensor_budgets: ContractionBudgets = ContractionBudgets(),
    native_budgets: LeanAlgebraBudgets = LeanAlgebraBudgets(),
    greedy_budgets: GreedyBudgets = GreedyBudgets(),
) -> StreamingCountResult:
    """Preflight the frontier, then replay it with finite algebraic leaf counts.

    Structural refusals happen before any leaf evaluator runs. A later leaf
    refusal records work already done but publishes neither count nor answer.
    Completed formulas, counts and branch histories are not stored.
    """
    limits = budgets or StreamingBudgets()
    if leaf_backend not in ('common', 'tensor', 'native', 'hash'):
        raise ValueError("leaf_backend must be 'common', 'tensor', 'native' or 'hash'")
    max_component_variables, vector_budget, total_vector_budget, memo_entries = (
        operator.index(value) for value in
        (max_component_variables, vector_budget, total_vector_budget, memo_entries)
    )
    if min(max_component_variables, vector_budget, total_vector_budget, memo_entries) < 0:
        raise ValueError('leaf budgets must be nonnegative integers')
    receipt = dict(n=cnf.n, leaf_backend=leaf_backend)
    plan = _Traversal()
    execution = _Traversal()
    work = dict(
        leaves_counted=0, boolean_vectors=0, factor_evaluations=0, local_evaluations=0,
        memo_hits=0, jet_products=0, peak_table_entries=0, peak_table_coefficients=0,
        accumulator_bits=1,
        native_queries=0, native_source_polynomials=0, monomial_products=0, scalar_operations=0,
        peak_polynomial_terms=0, peak_live_coefficients=0, peak_algebra_bits=0,
        max_ring_precision=0, hash_rounds=0, hash_candidates=0, hash_marked_queries=0,
        hash_accepted_queries=0, hash_history_rows=0, hash_initial_mass=0, hash_final_mass=0,
    )

    def finish(status: str, count: int | None, reason: str) -> StreamingCountResult:
        return StreamingCountResult(
            status, count, reason, **receipt, **work, nodes=execution.nodes,
            planning_nodes=plan.nodes, branches=plan.branches, leaves=plan.leaves,
            depth=plan.depth, max_pending_nodes=max(plan.max_pending, execution.max_pending),
            peak_frontier_bits=max(plan.peak_bits, execution.peak_bits),
            search_candidates=plan.candidates + execution.candidates,
            max_leaf_occurrence=plan.max_occurrence,
            count_word_bits=0 if count is None else max(1, count.bit_length()),
        )

    try:
        if cnf.n > limits.max_variables:
            raise StreamingRefusal('declared n exceeds the variable budget')
        if cnf.m > limits.max_clauses:
            raise StreamingRefusal('raw clause count exceeds its budget')
        size = formula_bit_length(cnf, canonical=False)
        receipt['input_bits'] = size
        if size > limits.max_formula_bits:
            raise StreamingRefusal('raw formula encoding exceeds its bit budget')
        root = _normalize(tuple(frozenset(clause) for clause in cnf.clauses))
        actual_width = max(map(len, root), default=0)
        width = max(2, actual_width) if k is None else operator.index(k)
        if width < 2 or actual_width > width:
            raise ValueError('k must be >= 2 and bound every normalized clause')
        receipt.update(width=width, normalized_clauses=len(root), normalized_bits=cnf.n + 1 + _body_bits(root))
        if width > limits.max_width:
            raise StreamingRefusal('clause width exceeds the core search budget')
        if M is None:
            try:
                parameter = choose_sparse_parameter(
                    width, sigma_numerator, sigma_denominator,
                    max_exponent=limits.max_parameter_exponent,
                )
            except ValueError as exc:
                if 'exponent budget' not in str(exc):
                    raise
                raise StreamingRefusal(str(exc)) from exc
            M = parameter.M
            receipt.update(parameter_test_passed=True, parameter_candidates=parameter.candidates_tested)
        else:
            M = operator.index(M)
            if M < 2:
                raise ValueError('M must be >= 2')
            if M.bit_length() - 1 > limits.max_parameter_exponent:
                raise StreamingRefusal('explicit M exceeds the parameter exponent budget')
        bound = occurrence_cap(width, M)
        receipt.update(M=M, occurrence_bound=bound)
        for _leaf in _leaves(root, cnf.n, M, bound, limits, plan):
            del _leaf
        total = 0
        stream = _leaves(root, cnf.n, M, bound, limits, execution)
        try:
            for leaf in stream:
                if leaf_backend == 'native':
                    result = count_cnf_ternary(_to_cnf(cnf.n, leaf), budgets=replace(
                        native_budgets,
                        max_monomial_products=native_budgets.max_monomial_products - work['monomial_products'],
                        max_scalar_operations=native_budgets.max_scalar_operations - work['scalar_operations'],
                    ))
                    work['native_queries'] += 1
                    work['native_source_polynomials'] += result.source_polynomials
                    work['local_evaluations'] += int(result.counted)
                    for name in ('monomial_products', 'scalar_operations'):
                        work[name] += getattr(result, name)
                    for name in ('peak_polynomial_terms', 'peak_live_coefficients', 'peak_algebra_bits'):
                        work[name] = max(work[name], getattr(result, name))
                    work['max_ring_precision'] = max(work['max_ring_precision'], result.precision)
                elif leaf_backend == 'hash':
                    result = count_cnf_greedy(
                        _to_cnf(cnf.n, leaf), depth=greedy_budgets.depth,
                        budgets=replace(
                            greedy_budgets,
                            max_marked_queries=(
                                greedy_budgets.max_marked_queries - work['hash_marked_queries']
                            ),
                        ),
                        native_budgets=replace(
                            native_budgets,
                            max_monomial_products=(
                                native_budgets.max_monomial_products - work['monomial_products']
                            ),
                            max_scalar_operations=(
                                native_budgets.max_scalar_operations - work['scalar_operations']
                            ),
                        ),
                    )
                    work['native_queries'] += 1
                    work['local_evaluations'] += int(result.counted)
                    work['hash_rounds'] += result.rounds
                    work['hash_candidates'] += result.candidates_scanned
                    work['hash_marked_queries'] += result.marked_queries
                    work['hash_accepted_queries'] += result.accepted_queries
                    work['hash_history_rows'] += result.history_rows
                    work['hash_initial_mass'] += result.initial_mass
                    work['hash_final_mass'] += result.final_mass
                    work['native_source_polynomials'] += result.source_polynomials
                    for name in ('monomial_products', 'scalar_operations'):
                        work[name] += getattr(result, name)
                    for name in (
                        'peak_polynomial_terms', 'peak_live_coefficients', 'peak_algebra_bits'
                    ):
                        work[name] = max(work[name], getattr(result, name))
                    work['max_ring_precision'] = max(work['max_ring_precision'], result.precision)
                else:
                    remaining = total_vector_budget - work['boolean_vectors']
                    result = count_cnf_components(
                        _to_cnf(cnf.n, leaf), max_variables=limits.max_variables,
                        max_component_variables=max_component_variables,
                        vector_budget=min(vector_budget, remaining), total_vector_budget=remaining,
                        memo_entries=memo_entries, placement_budgets=placement_budgets,
                        backend=leaf_backend, tensor_budgets=replace(
                            tensor_budgets, max_boolean_vectors=min(tensor_budgets.max_boolean_vectors, remaining),
                        ),
                    )
                    for name in ('boolean_vectors', 'factor_evaluations', 'local_evaluations', 'memo_hits', 'jet_products'):
                        work[name] += getattr(result, name)
                    for name in ('peak_table_entries', 'peak_table_coefficients'):
                        work[name] = max(work[name], getattr(result, name))
                if not result.counted:
                    reason = result.reason.removesuffix('; no partial count returned')
                    raise StreamingRefusal(f'algebraic leaf refused: {reason}')
                value = result.count
                if value is None or value < 0 or value.bit_length() > cnf.n + 1 or value > (1 << cnf.n):
                    raise RuntimeError('an algebraic leaf returned an invalid exact count')
                total += value
                if total > (1 << cnf.n):
                    raise RuntimeError('the disjoint count accumulator exceeds 2^n')
                work['accumulator_bits'] = max(work['accumulator_bits'], max(1, total.bit_length()))
                work['leaves_counted'] += 1
                del leaf, result
        finally:
            stream.close()
        return finish('counted', total, 'streamed disjoint leaf counts completed on the original declared n')
    except StreamingRefusal as exc:
        return finish('refused', None, f'{exc}; no partial count returned')
