"""Finite Section 5.3–5.5 replacement sum in the ring A_ell.

Each sparse equation becomes the factor in (5.4):

    B_i = sum_{j < ell} e_j(L_i) h(L_i)^j V^j

in A_ell = Q[V]/(V^ell), with the linear polynomials e_j from Lemma 5.1 and
h(Y) = Y^2 - Y. Neutral factors equal to 1, each with scope {v}, make every
declared variable occur. A finite design meeting Definition 3.1 supplies one
replacement row per factor–variable incidence. Proposition 3.4 then sums the
product by the Section 3 measurement scan, with coefficients in
A_ell[z]/(z^J). The least nonnegative residue of S_ell(H(0)) is N(P) when
2^ell > 2^n.

The rows used here are an explicit finite design: P(0) = delta_v and the
degree is at most R. They are not the Section 4 auxiliary array and they are
not produced by Lemma 4.8. On every finite branch checked by this module the
Lemma 4.8 preliminary guard fails, so tau, epsilon and D are not formed.
Section 6 equation selection, Appendix D, and Theorem 2.2's bit-complexity
bound are not run. ``theorem_1_1_complete`` stays false.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass, fields
from fractions import Fraction
from typing import Mapping, Sequence

from .common_count import _substitutions, separated_group_sum
from .hermite import _jet_mul, _one_jet, hermite_eval, polynomial_hasse_jet
from .linear_algebra import SparseRationalRowSpan
from .reference_parameters import lemma_4_8_proposal
from .sparse_systems import (
    BudgetRefusal, ParityPolynomials, Section5Budgets, SparsePolynomial, SparseSystem,
    construct_parity_polynomials,
)

STEP18_SCOPE = (
    'Finite Section 5.3–5.5 replacement sum in A_ell. Equation factors use '
    'Lemma 5.1 and the total h-power truncation. A checked finite design '
    'meets Definition 3.1, and the Section 3 routine sums the product in '
    'A_ell[z]/(z^J). Lemma 4.8 does not form D on these branches. Section 6, '
    'Appendix D, and Theorem 2.2 are not run. This does not implement '
    'Theorem 1.1 and it is not a 2^{o(n)} algorithm.'
)


@dataclass(frozen=True)
class AlgebraBudgets:
    """Caps for the replacement sum. Each refusal happens before the corresponding allocation."""

    max_ell: int = 24
    max_coefficient_bits: int = 200_000
    max_boolean_vectors: int = 100_000
    max_substitution_pairs: int = 20_000
    max_jet_length: int = 64
    max_algebra_slots: int = 100_000
    max_accumulator_bits: int = 200_000
    max_branches: int = 64
    max_equations: int = 32

    def __post_init__(self):
        for field in fields(self):
            value = operator.index(getattr(self, field.name))
            if value < 0:
                raise ValueError('algebra budgets must be nonnegative integers')
            object.__setattr__(self, field.name, value)


@dataclass(frozen=True)
class PolynomialRow:
    """One replacement row P(T) in Q[T]^{1 x n}. ``terms[s]`` is the coefficient of T^s."""

    terms: tuple[dict[int, Fraction], ...]

    @property
    def degree(self) -> int:
        return len(self.terms) - 1


@dataclass(frozen=True)
class EquationFactor:
    """A Section 5.3 factor. Neutral factors are the constant 1 with a singleton scope."""

    kind: str
    scope: tuple[int, ...]
    monomials: tuple[frozenset[int], ...]
    rows: tuple[PolynomialRow, ...]


@dataclass(frozen=True)
class AlgebraDesign:
    """Definition 3.1 data for one finite sum. Points are indexed from 0 and evaluated at 1..Q."""

    n: int
    r: int
    Q: int
    R: int
    D: int
    groups: tuple[tuple[tuple[int, ...], ...], ...]
    assignments: tuple[tuple[int, ...], ...]
    measurements: tuple[tuple[tuple[dict[int, int], ...], ...], ...]
    factors: tuple[EquationFactor, ...]

    @property
    def equation_count(self) -> int:
        return sum(factor.kind == 'equation' for factor in self.factors)

    @property
    def neutral_count(self) -> int:
        return sum(factor.kind == 'neutral' for factor in self.factors)


@dataclass(frozen=True)
class IncidenceCut:
    """The high-incidence variables of Section 5.5, after discarding 0 = 0 and adding neutral scopes."""

    n: int
    incidence_total: int
    log_n: int
    sqrt_log_n: int
    occurrence_cap: int
    heavy: tuple[int, ...]
    equation_scopes: tuple[tuple[int, ...], ...]
    neutral_scopes: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class AlgebraPlan:
    input_degree_bound: int
    h_degree_bound: int
    q: int
    J: int
    substitution_pairs: int
    boolean_vector_bound: int
    algebra_slots: int
    accumulator_bit_bound: int = 0


@dataclass(frozen=True)
class ReplacementSumResult:
    status: str
    count: int | None
    reason: str
    ell: int
    boolean_vectors: int = 0
    factor_evaluations: int = 0
    vector_bound: int = 0
    branches: int = 0
    branches_enumerated: int = 0
    branch_counts: tuple[int, ...] = ()
    heavy_variables: tuple[int, ...] = ()
    equations_used: int = 0
    neutral_factors: int = 0
    input_degree_bound: int = 0
    h_degree_bound: int = 0
    q: int = 0
    jet_length: int = 0
    design_Q: int = 0
    replacement_degree: int = 0
    substitution_pairs: int = 0
    algebra_slots: int = 0
    h0_coefficients: tuple[int, ...] = ()
    unreduced_sum: int | None = None
    sample_at_one: str | None = None
    slope: str = '0'
    lemma_48_guard_passed: bool | None = None
    lemma_48_dimension_formed: bool = False
    lemma_48_r: int | None = None
    lemma_48_Q: int | None = None
    lemma_48_R: int | None = None
    jets_allocated: bool = False
    coefficients_constructed: bool = False
    enumerated: bool = False
    occurrence_cap: int = 0
    note: str = STEP18_SCOPE

    @property
    def counted(self) -> bool:
        return self.status == 'counted'

    @property
    def refused(self) -> bool:
        return self.status == 'refused'


def input_degree_bound(k: int, m_test: int, ell: int) -> int:
    """Lemma 5.2: ``d = k (m_test + 2(ell - 1))`` for the product in A_ell."""
    degree = operator.index(k)
    equations = operator.index(m_test)
    width = operator.index(ell)
    if degree < 0 or equations < 0 or width < 1:
        raise ValueError('degree bound expects k >= 0, m_test >= 0 and ell >= 1')
    return degree * (equations + 2 * (width - 1))


def h_degree_bound(k: int, m_test: int, ell: int, R: int) -> int:
    """Equation (5.5): ``d_H = k R (m_test + 2(ell - 1))``."""
    replacement = operator.index(R)
    if replacement < 0:
        raise ValueError('R must be nonnegative')
    return replacement * input_degree_bound(k, m_test, ell)


def incidence_cut(system: SparseSystem) -> IncidenceCut:
    """Section 5.5's occurrence cap, including one neutral scope per declared variable.

    ``b0 = ceil((M / n) * ceil(sqrt(lambda)))`` with ``lambda = ceil(log2 n)``.
    A variable is heavy when it belongs to more than ``b0`` scopes. Discard
    equations ``0 = 0`` before forming the scopes. Constant ``1 = 0`` is left
    for the caller and is not a scope.
    """
    variables = system.n
    equation_scopes = []
    for equation in system.equations:
        if equation.is_zero or equation.is_one:
            continue
        equation_scopes.append(tuple(sorted(equation.support)))
    neutral_scopes = tuple((variable,) for variable in range(variables))
    scopes = tuple(equation_scopes) + neutral_scopes
    total = sum(len(scope) for scope in scopes)
    log_n = _ceil_log2(variables) if variables >= 1 else 0
    root = _ceil_sqrt(log_n)
    if variables == 0:
        cap = 0
        heavy: tuple[int, ...] = ()
    else:
        cap = (total * root + variables - 1) // variables if root else 0
        counts = [0] * variables
        for scope in scopes:
            for variable in scope:
                counts[variable] += 1
        heavy = tuple(index for index, count in enumerate(counts) if count > cap)
    return IncidenceCut(
        variables, total, log_n, root, cap, heavy,
        tuple(equation_scopes), neutral_scopes,
    )


def restrict_system(system: SparseSystem, fixed: Mapping[int, int]) -> SparseSystem:
    """Substitute Boolean bits and relabel the free variables by ``[n']``.

    Duplicate monomials cancel modulo 2. Neither degree nor the number of
    monomials increases. The original allowance is not rechecked against the
    smaller ``n'``: ``(log2 n')^c`` can shrink even though the proof's lists
    do not grow.
    """
    assignment = {}
    for raw_variable, raw_bit in fixed.items():
        variable = operator.index(raw_variable)
        bit = operator.index(raw_bit)
        if variable < 0 or variable >= system.n:
            raise ValueError(f'fixed variable {variable} is outside 0..{system.n - 1}')
        if bit not in (0, 1):
            raise ValueError('fixed bits must be 0 or 1')
        assignment[variable] = bit
    free = [variable for variable in range(system.n) if variable not in assignment]
    relabel = {old: new for new, old in enumerate(free)}
    equations = []
    for equation in system.equations:
        present: set[frozenset[int]] = set()
        for monomial in equation.monomials:
            include = True
            kept = []
            for variable in monomial:
                if variable in assignment:
                    if assignment[variable] == 0:
                        include = False
                        break
                else:
                    kept.append(relabel[variable])
            if not include:
                continue
            reduced = frozenset(kept)
            if reduced in present:
                present.remove(reduced)
            else:
                present.add(reduced)
        ordered = tuple(sorted(present, key=lambda item: (len(item), tuple(sorted(item)))))
        equations.append(SparsePolynomial(len(free), ordered))
    return SparseSystem(len(free), system.k, system.c, tuple(equations))


def finite_admitted_design(
    system: SparseSystem,
    *,
    R: int,
    Q: int,
    slope: Fraction | None = None,
) -> AlgebraDesign:
    """One block contains every variable. The other block is empty. Every factor uses the full block.

    ``P_{f,v}(T) = (1 + slope * T) e_v`` when ``slope`` is nonzero, and ``e_v``
    when it is zero. The default slope is ``1/Q``, so ``(R-1)! Q!`` clears the
    denominator. This meets Definition 3.1 for these scopes. It is not a
    Lemma 4.8 design.
    """
    replacement = operator.index(R)
    points = operator.index(Q)
    if replacement < 1 or points <= replacement:
        raise ValueError('Definition 3.1 needs integers 1 <= R < Q')
    if system.n < 1:
        raise ValueError('the finite design needs at least one variable')
    chosen = Fraction(1, points) if slope is None else Fraction(slope)
    factors = _active_factors(system, replacement, chosen)
    if factors is None:
        raise ValueError('a constant equation 1 = 0 has no replacement design')
    if not any(factor.kind == 'equation' for factor in factors):
        raise ValueError('an empty equation list has no replacement design')
    groups = tuple(
        (tuple(range(system.n)), tuple())
        for _ in range(points)
    )
    assignments = tuple(tuple(0 for _ in range(points)) for _ in factors)
    measurements = tuple(tuple(() for _ in range(2)) for _ in range(points))
    return AlgebraDesign(
        system.n, 2, points, replacement, 1, groups, assignments, measurements, tuple(factors),
    )


def grouped_unit_design(
    system: SparseSystem,
    groups: Sequence[Sequence[int]],
    equation_groups: Sequence[int],
    *,
    R: int,
    Q: int,
) -> AlgebraDesign:
    """Repeat one partition at every point and use the unit rows ``P = e_v``.

    Each neutral variable is assigned to the unique block that contains it.
    ``equation_groups`` assigns the non-constant equations, in order. The
    kernel check, not this constructor, decides whether the assignment is
    admitted.
    """
    replacement = operator.index(R)
    points = operator.index(Q)
    if replacement < 1 or points <= replacement:
        raise ValueError('Definition 3.1 needs integers 1 <= R < Q')
    partition = tuple(tuple(block) for block in groups)
    if len(partition) < 2:
        raise ValueError('Definition 3.1 asks for r >= 2, and empty blocks are allowed')
    equations = [equation for equation in system.equations if not equation.is_zero]
    if any(equation.is_one for equation in equations):
        raise ValueError('a constant equation 1 = 0 has no replacement design')
    if len(equation_groups) != len(equations):
        raise ValueError('one group index is required per non-constant equation')
    owner = {}
    for index, block in enumerate(partition):
        for variable in block:
            if variable in owner:
                raise ValueError(f'variable {variable} lies in two blocks')
            owner[variable] = index
    if set(owner) != set(range(system.n)):
        raise ValueError('the blocks must partition every declared variable')
    factors = []
    for equation, group in zip(equations, equation_groups):
        group_index = operator.index(group)
        if group_index < 0 or group_index >= len(partition):
            raise ValueError('equation group index is outside the partition')
        scope = tuple(sorted(equation.support))
        factors.append(_equation_factor(scope, equation.monomials, system.n, replacement, Fraction(0)))
    for variable in range(system.n):
        factors.append(_neutral_factor(variable, system.n, replacement, Fraction(0)))
    assignments = []
    for equation, group in zip(equations, equation_groups):
        assignments.append(tuple(operator.index(group) for _ in range(points)))
    for variable in range(system.n):
        assignments.append(tuple(owner[variable] for _ in range(points)))
    repeated = tuple(partition for _ in range(points))
    measurements = tuple(
        tuple(tuple() for _ in partition)
        for _ in range(points)
    )
    return AlgebraDesign(
        system.n, len(partition), points, replacement, 1,
        repeated, tuple(assignments), measurements, tuple(factors),
    )


def design_problems(design: AlgebraDesign) -> tuple[str, ...]:
    """Exact Definition 3.1 checks for a supplied finite design. No jet is allocated."""
    problems = []
    if design.r < 2:
        problems.append(f'r={design.r} is below 2')
    if not 1 <= design.R < design.Q:
        problems.append(f'the design needs 1 <= R < Q; got R={design.R}, Q={design.Q}')
    if design.D < 1:
        problems.append('D must be an integer >= 1')
    if len(design.groups) != design.Q or len(design.measurements) != design.Q:
        problems.append('each of the Q points needs a partition and measurement matrices')
    if len(design.factors) != len(design.assignments):
        problems.append('each factor needs an assignment at every point')
    denominator = math.factorial(max(design.R - 1, 0)) * math.factorial(max(design.Q, 0))
    for point, (partition, matrices) in enumerate(zip(design.groups, design.measurements), start=1):
        if len(partition) != design.r or len(matrices) != design.r:
            problems.append(f'point {point} does not have r blocks')
            continue
        seen: set[int] = set()
        for group_index, (block, rows) in enumerate(zip(partition, matrices)):
            if len(block) * design.r > 2 * design.n:
                problems.append(
                    f'point {point} block {group_index} has {len(block)} variables; '
                    f'Definition 3.1 allows at most 2n/r'
                )
            for variable in block:
                if variable in seen or variable < 0 or variable >= design.n:
                    problems.append(f'point {point} is not a partition of 0..{design.n - 1}')
                    break
                seen.add(variable)
            for row in rows:
                absolute = 0
                for local, value in row.items():
                    coefficient = Fraction(value)
                    if coefficient.denominator != 1:
                        problems.append('measurement entry is not an integer')
                    if local < 0 or local >= len(block):
                        problems.append('measurement column is outside its block')
                    absolute += abs(int(coefficient))
                if absolute > design.D:
                    problems.append(f'measurement row has absolute sum {absolute} above D={design.D}')
        if seen != set(range(design.n)) and design.n:
            problems.append(f'point {point} misses a declared variable')
    for index, (factor, assignment) in enumerate(zip(design.factors, design.assignments)):
        if len(assignment) != design.Q:
            problems.append(f'factor {index} is not assigned at every point')
        if any(group < 0 or group >= design.r for group in assignment):
            problems.append(f'factor {index} is assigned outside 0..r-1')
        if factor.kind not in ('equation', 'neutral'):
            problems.append(f'factor {index} has unknown kind {factor.kind}')
        if len(factor.rows) != len(factor.scope):
            problems.append(f'factor {index} does not have one row per scope variable')
        if factor.kind == 'neutral' and (len(factor.scope) != 1 or factor.monomials):
            problems.append(f'neutral factor {index} must be the constant 1 on a singleton scope')
        scope = set(factor.scope)
        for monomial in factor.monomials:
            if not monomial <= scope:
                problems.append(f'factor {index} has a monomial outside its scope')
        for variable, row in zip(factor.scope, factor.rows):
            if row.degree > design.R or not row.terms:
                problems.append(f'replacement degree exceeds R={design.R}')
                continue
            constant = row.terms[0]
            if constant != {variable: Fraction(1)}:
                problems.append(f'factor {index} does not restore variable {variable} at T=0')
            for term in row.terms:
                for coordinate, coefficient in term.items():
                    if coordinate < 0 or coordinate >= design.n:
                        problems.append('replacement coordinate is outside 0..n-1')
                    if (denominator * Fraction(coefficient)).denominator != 1:
                        problems.append('replacement coefficient is not integral after clearing (R-1)! Q!')
    problems.extend(_kernel_problems(design))
    return tuple(problems)


def count_replacement_sum(
    system: SparseSystem,
    budgets: AlgebraBudgets | None = None,
    *,
    ell: int | None = None,
    slope: Fraction | None = None,
    Q: int | None = None,
    R: int = 1,
) -> ReplacementSumResult:
    """Section 5.5's branch split, then the Section 3 sum in A_ell on each branch.

    An empty list is ``2^n`` and a constant equation ``1 = 0`` is 0, both
    before any design or scan. Heavy variables are fixed one branch at a time.
    Every branch is planned before the first scan. A budget failure returns
    ``count is None`` and scans nothing. Lemma 4.8 is consulted only to record
    that its preliminary guard does not form D; the rows actually summed are
    the finite design.
    """
    limits = budgets or AlgebraBudgets()
    width = system.n + 1 if ell is None else operator.index(ell)
    if width < 1:
        raise ValueError('ell must be a positive integer')
    if width <= system.n:
        raise ValueError('ell must be at least n+1 so the residue equals N(P)')
    if any(equation.is_one for equation in system.equations):
        return _counted(
            0, width, reason='constant equation 1 = 0 has no Boolean solutions',
        )
    active = [equation for equation in system.equations if not equation.is_zero]
    if not active:
        return _counted(
            1 << system.n, width, reason='empty equation list has 2^n common zeros',
        )
    cut = incidence_cut(system)
    branch_count = 1 << len(cut.heavy)
    if branch_count > limits.max_branches:
        return _refused(
            width,
            reason=(
                f'branch count 2^{len(cut.heavy)} exceeds budget {limits.max_branches}; '
                'no branch enumerated'
            ),
            heavy=cut.heavy, branches=branch_count, occurrence_cap=cut.occurrence_cap,
        )
    if len(active) > limits.max_equations:
        return _refused(
            width,
            reason=(
                f'm_test={len(active)} exceeds budget {limits.max_equations}; '
                'replacement rows not allocated'
            ),
            heavy=cut.heavy, branches=branch_count, occurrence_cap=cut.occurrence_cap,
        )
    replacement = operator.index(R)
    prepared = []
    for mask in range(branch_count):
        fixed = {variable: (mask >> index) & 1 for index, variable in enumerate(cut.heavy)}
        restricted = restrict_system(system, fixed)
        if any(equation.is_one for equation in restricted.equations):
            prepared.append(('zero', restricted, None, None, None))
            continue
        equations = [equation for equation in restricted.equations if not equation.is_zero]
        if not equations:
            prepared.append(('power', restricted, None, None, None))
            continue
        degree = h_degree_bound(restricted.k, len(equations), width, replacement)
        points = max(replacement + 1, degree + 1) if Q is None else operator.index(Q)
        try:
            design = finite_admitted_design(restricted, R=replacement, Q=points, slope=slope)
        except ValueError as exc:
            return _invalid(width, str(exc), heavy=cut.heavy, branches=branch_count)
        problems = design_problems(design)
        if problems:
            raise RuntimeError(f'finite admitted design failed Definition 3.1: {problems[0]}')
        plan = _plan(design, restricted.k, len(equations), width)
        proposal = lemma_4_8_proposal(restricted.n, [factor.scope for factor in design.factors])
        prepared.append(('algebra', restricted, design, plan, proposal))
    for kind, _restricted, _design, plan, _proposal in prepared:
        if kind != 'algebra':
            continue
        refusal = _budget_refusal(plan, limits, width, cut, branch_count)
        if refusal is not None:
            return refusal
    algebra = [item for item in prepared if item[0] == 'algebra']
    parity = None
    if algebra:
        if width > limits.max_ell:
            return _refused(
                width,
                reason=f'ell={width} exceeds budget {limits.max_ell}; parity coefficients not constructed',
                heavy=cut.heavy, branches=branch_count, occurrence_cap=cut.occurrence_cap,
                **_plan_fields(algebra[0][3]),
            )
        try:
            parity = construct_parity_polynomials(width, _parity_budgets(limits))
        except BudgetRefusal as exc:
            return _refused(
                width, reason=str(exc), heavy=cut.heavy, branches=branch_count,
                occurrence_cap=cut.occurrence_cap, **_plan_fields(algebra[0][3]),
            )
        for _kind, _restricted, design, plan, _proposal in algebra:
            bound = _accumulator_bit_bound(design, parity, plan)
            if bound > limits.max_accumulator_bits:
                return _refused(
                    width,
                    reason=(
                        f'accumulator bit bound {bound} exceeds budget {limits.max_accumulator_bits}; '
                        'jets not allocated and no assignments scanned'
                    ),
                    heavy=cut.heavy, branches=branch_count, occurrence_cap=cut.occurrence_cap,
                    coefficients_constructed=True, coefficient_bits=parity.coefficient_bits,
                    **_plan_fields(plan),
                )
    total = 0
    vectors = 0
    evaluations = 0
    enumerated = 0
    branch_counts = []
    jets_allocated = False
    sample = None
    coefficients: tuple[int, ...] = ()
    unreduced = None
    guard_flags: list[bool] = []
    dimension_formed = False
    proposal_r = proposal_Q = proposal_R = None
    used = 0
    neutrals = 0
    degree_bound = 0
    replaced_degree = 0
    query = 0
    jet_length = 0
    design_points = 0
    replacement_degree = 0
    pairs = 0
    slots = 0
    vector_bound = 0
    chosen_slope = '0'
    for kind, restricted, design, plan, proposal in prepared:
        if kind == 'zero':
            branch_counts.append(0)
            continue
        if kind == 'power':
            branch_counts.append(1 << restricted.n)
            total += 1 << restricted.n
            continue
        assert parity is not None and design is not None and plan is not None and proposal is not None
        outcome = _enumerate_design(design, parity, plan)
        jets_allocated = True
        enumerated += 1
        vectors += outcome[0]
        evaluations += outcome[1]
        branch_counts.append(outcome[2])
        total += outcome[2]
        if sample is None:
            sample = outcome[3]
            coefficients = outcome[4]
            unreduced = outcome[5]
            used = design.equation_count
            neutrals = design.neutral_count
            degree_bound = plan.input_degree_bound
            replaced_degree = plan.h_degree_bound
            query = plan.q
            jet_length = plan.J
            design_points = design.Q
            replacement_degree = design.R
            pairs = plan.substitution_pairs
            slots = plan.algebra_slots
            chosen_slope = _slope_text(design)
            proposal_r, proposal_Q, proposal_R = proposal.r, proposal.Q, proposal.R
        guard_flags.append(proposal.r >= 2 and proposal.Q > proposal.R)
        dimension_formed = dimension_formed or bool(proposal.dimension_formed)
        vector_bound += plan.boolean_vector_bound
    if total < 0 or total > 1 << system.n:
        raise ArithmeticError('branch total is outside 0..2^n')
    modulus = 1 << width
    if total >= modulus:
        raise ArithmeticError('branch total is not determined by its residue modulo 2^ell')
    return ReplacementSumResult(
        status='counted',
        count=total,
        reason='least nonnegative residue of S_ell(H(0)) summed over Section 5.5 branches',
        ell=width,
        boolean_vectors=vectors,
        factor_evaluations=evaluations,
        vector_bound=vector_bound,
        branches=branch_count,
        branches_enumerated=enumerated,
        branch_counts=tuple(branch_counts),
        heavy_variables=cut.heavy,
        equations_used=used,
        neutral_factors=neutrals,
        input_degree_bound=degree_bound,
        h_degree_bound=replaced_degree,
        q=query,
        jet_length=jet_length,
        design_Q=design_points,
        replacement_degree=replacement_degree,
        substitution_pairs=pairs,
        algebra_slots=slots,
        h0_coefficients=coefficients,
        unreduced_sum=unreduced,
        sample_at_one=sample,
        slope=chosen_slope,
        lemma_48_guard_passed=all(guard_flags) if guard_flags else None,
        lemma_48_dimension_formed=dimension_formed,
        lemma_48_r=proposal_r,
        lemma_48_Q=proposal_Q,
        lemma_48_R=proposal_R,
        jets_allocated=jets_allocated,
        coefficients_constructed=parity is not None,
        enumerated=enumerated > 0,
        occurrence_cap=cut.occurrence_cap,
        note=STEP18_SCOPE,
    )


def sum_supplied_design(
    system: SparseSystem,
    design: AlgebraDesign,
    budgets: AlgebraBudgets | None = None,
    *,
    ell: int | None = None,
) -> ReplacementSumResult:
    """Proposition 3.4 in A_ell for one supplied design. No Section 5.5 branch split.

    The design's equation factors must be exactly the system's non-constant
    equations, in order, followed by one neutral factor for each variable.
    An invalid design returns no count and starts no scan.
    """
    limits = budgets or AlgebraBudgets()
    width = system.n + 1 if ell is None else operator.index(ell)
    if width < 1:
        raise ValueError('ell must be a positive integer')
    if width <= system.n:
        raise ValueError('ell must be at least n+1 so the residue equals N(P)')
    if design.n != system.n:
        return _invalid(width, 'design dimension does not equal n')
    equations = [equation for equation in system.equations if not equation.is_zero]
    if any(equation.is_one for equation in equations):
        return _counted(0, width, reason='constant equation 1 = 0 has no Boolean solutions')
    if not equations:
        return _counted(1 << system.n, width, reason='empty equation list has 2^n common zeros')
    expected = _active_factors(system, design.R, Fraction(0))
    if expected is None:
        return _counted(0, width, reason='constant equation 1 = 0 has no Boolean solutions')
    if len(design.factors) != len(expected):
        return _invalid(width, 'supplied design does not have one factor per equation and variable')
    for got, want in zip(design.factors, expected):
        if got.kind != want.kind or got.scope != want.scope or got.monomials != want.monomials:
            return _invalid(width, 'supplied factors do not match the sparse system')
    problems = design_problems(design)
    if problems:
        return _invalid(width, problems[0])
    plan = _plan(design, system.k, len(equations), width)
    empty_cut = IncidenceCut(system.n, 0, 0, 0, 0, (), (), ())
    refusal = _budget_refusal(plan, limits, width, empty_cut, 1)
    if refusal is not None:
        return refusal
    if width > limits.max_ell:
        return _refused(
            width,
            reason=f'ell={width} exceeds budget {limits.max_ell}; parity coefficients not constructed',
            **_plan_fields(plan),
        )
    try:
        parity = construct_parity_polynomials(width, _parity_budgets(limits))
    except BudgetRefusal as exc:
        return _refused(width, reason=str(exc), **_plan_fields(plan))
    bound = _accumulator_bit_bound(design, parity, plan)
    if bound > limits.max_accumulator_bits:
        return _refused(
            width,
            reason=(
                f'accumulator bit bound {bound} exceeds budget {limits.max_accumulator_bits}; '
                'jets not allocated and no assignments scanned'
            ),
            coefficients_constructed=True,
            **_plan_fields(plan),
        )
    vectors, evaluations, count, sample, coefficients, unreduced = _enumerate_design(design, parity, plan)
    return ReplacementSumResult(
        status='counted',
        count=count,
        reason='least nonnegative residue of S_ell(H(0)) for the supplied design',
        ell=width,
        boolean_vectors=vectors,
        factor_evaluations=evaluations,
        vector_bound=plan.boolean_vector_bound,
        branches=1,
        branches_enumerated=1,
        branch_counts=(count,),
        equations_used=design.equation_count,
        neutral_factors=design.neutral_count,
        input_degree_bound=plan.input_degree_bound,
        h_degree_bound=plan.h_degree_bound,
        q=plan.q,
        jet_length=plan.J,
        design_Q=design.Q,
        replacement_degree=design.R,
        substitution_pairs=plan.substitution_pairs,
        algebra_slots=plan.algebra_slots,
        h0_coefficients=coefficients,
        unreduced_sum=unreduced,
        sample_at_one=sample,
        slope=_slope_text(design),
        jets_allocated=True,
        coefficients_constructed=True,
        enumerated=True,
        note=STEP18_SCOPE,
    )


def _enumerate_design(design: AlgebraDesign, parity: ParityPolynomials, plan: AlgebraPlan):
    """Allocate the jets and run Proposition 3.4. Callers must budget-check first."""
    nodes = tuple(_node(design, point) for point in range(design.Q))
    prepared = tuple(_prepare_point(design, point, plan.q) for point in range(design.Q))
    assigned = []
    for point in range(design.Q):
        groups = [[] for _ in range(design.r)]
        for index, assignment in enumerate(design.assignments):
            groups[assignment[point]].append(index)
        assigned.append(tuple(tuple(group) for group in groups))
    jets = []
    stats = [0, 0, 0]
    zero = _ajet_zero(parity.ell, plan.J)
    one = _ajet_one(parity.ell, plan.J)
    for point, node in enumerate(nodes):
        jet = zero
        point_assigned = assigned[point]
        point_inputs = prepared[point]
        # U is fixed before the measurement scan, as in (3.13).
        for E, U, sign in _substitutions(design.n, plan.q):
            def evaluate_U(group, bits, _E=E, _U=U, _assigned=point_assigned, _inputs=point_inputs):
                del _E
                product = one
                for index in _assigned[group]:
                    stats[2] += 1
                    factor = design.factors[index]
                    if factor.kind == 'neutral':
                        continue
                    values = {}
                    for local, variable in enumerate(factor.scope):
                        values[variable] = _input_jet(_inputs[index][local], bits, _U, plan.J)
                    product = _ajet_mul(product, _equation_jet(factor.monomials, values, parity, plan.J))
                return product

            value = separated_group_sum(
                design.n, node[0], node[1], node[2], node[3], set(E), evaluate_U,
                zero, one, _ajet_add, _ajet_mul, stats,
            )
            jet = _ajet_add(jet, _ajet_scale(value, sign))
        jets.append(jet)
    recovered = []
    for slot in range(parity.ell):
        samples = []
        for point, jet in enumerate(jets, start=1):
            samples.append((Fraction(point), tuple(coefficient[slot] for coefficient in jet)))
        recovered.append(hermite_eval(samples, Fraction(0)))
    if any(coefficient.denominator != 1 for coefficient in recovered):
        raise ArithmeticError('S_ell coefficients at T=0 are not integers')
    integers = tuple(int(coefficient) for coefficient in recovered)
    unreduced = sum(integers)
    modulus = 1 << parity.ell
    residue = unreduced % modulus
    if residue > 1 << design.n:
        raise ArithmeticError('restored residue is outside 0..2^n')
    sample = _ajet_sum(jets[0][0]) if jets else Fraction(0)
    if stats[0] > plan.boolean_vector_bound:
        raise ArithmeticError('local enumeration exceeded its planned vector bound')
    return stats[0], stats[2], residue, str(sample), integers, unreduced


def _plan(design: AlgebraDesign, k: int, m_test: int, ell: int) -> AlgebraPlan:
    degree = input_degree_bound(k, m_test, ell)
    replaced = design.R * degree
    q = replaced // design.Q
    length = q + 1
    pairs = _substitution_count(design.n, q)
    bound = 0
    for point in range(design.Q):
        partition = design.groups[point]
        matrices = design.measurements[point]
        tuples = 1
        cubes = 0
        for block, rows in zip(partition, matrices):
            span = SparseRationalRowSpan(tuple(dict(row) for row in rows), len(block))
            for row in span.independent_rows:
                lo = int(sum(value for value in row.values() if value < 0))
                hi = int(sum(value for value in row.values() if value > 0))
                tuples *= hi - lo + 1
            cubes += 1 << len(block)
        bound += pairs * tuples * 2 * cubes
    return AlgebraPlan(
        degree, replaced, q, length, pairs, bound, design.Q * length * ell,
    )


def _budget_refusal(plan, limits, ell, cut, branch_count):
    checks = (
        (plan.boolean_vector_bound, limits.max_boolean_vectors, 'Boolean vector',
         'no enumeration started'),
        (plan.substitution_pairs, limits.max_substitution_pairs, 'substitution pair',
         'no enumeration started'),
        (plan.J, limits.max_jet_length, 'jet length', 'jets not allocated'),
        (plan.algebra_slots, limits.max_algebra_slots, 'algebra slot', 'jets not allocated'),
    )
    for value, limit, name, stop in checks:
        if value > limit:
            return _refused(
                ell,
                reason=f'{name} bound {value} exceeds budget {limit}; {stop}',
                heavy=cut.heavy, branches=branch_count, occurrence_cap=cut.occurrence_cap,
                **_plan_fields(plan),
            )
    return None


def _node(design: AlgebraDesign, point: int):
    groups = design.groups[point]
    matrices = []
    ranges = []
    for block, rows in zip(groups, design.measurements[point]):
        independent = SparseRationalRowSpan(tuple(dict(row) for row in rows), len(block)).independent_rows
        matrices.append(independent)
        for row in independent:
            lo = int(sum(value for value in row.values() if value < 0))
            hi = int(sum(value for value in row.values() if value > 0))
            ranges.append(range(int(lo), int(hi) + 1))
    return groups, tuple(matrices), tuple(ranges), tuple(map(len, matrices))


def _prepare_point(design: AlgebraDesign, point: int, order: int):
    base = Fraction(point + 1)
    prepared = []
    for factor in design.factors:
        factor_inputs = []
        for row in factor.rows:
            coordinates = []
            width = row.degree + 1
            for variable in range(design.n):
                coeffs = [Fraction(term.get(variable, 0)) for term in row.terms]
                if any(coeffs):
                    coordinates.append((variable, polynomial_hasse_jet(coeffs, base, order)))
            del width
            factor_inputs.append(tuple(coordinates))
        prepared.append(tuple(factor_inputs))
    return tuple(prepared)


def _input_jet(coordinates, bits, U, length: int):
    jet = [Fraction(0)] * length
    fixed = set(U)
    for variable, hasse in coordinates:
        if bits[variable]:
            jet[0] += hasse[0]
        if variable in fixed:
            for power in range(1, length):
                jet[power] += hasse[power]
    return jet


def _equation_jet(monomials, inputs, parity: ParityPolynomials, length: int):
    lift = [Fraction(0)] * length
    for monomial in monomials:
        term = _one_jet(length - 1)
        for variable in monomial:
            term = _jet_mul(term, inputs[variable])
        for index, coefficient in enumerate(term):
            lift[index] += coefficient
    height_argument = list(lift)
    height_argument[0] -= 1
    height = _jet_mul(lift, height_argument)
    result = [_alg_zero(parity.ell) for _ in range(length)]
    height_power = _one_jet(length - 1)
    for power, (constant, slope) in enumerate(parity.coefficients):
        linear = [slope * coefficient for coefficient in lift]
        linear[0] += constant
        term = _jet_mul(linear, height_power)
        for index, coefficient in enumerate(term):
            if coefficient:
                slot = list(result[index])
                slot[power] += coefficient
                result[index] = tuple(slot)
        if power + 1 < parity.ell:
            height_power = _jet_mul(height_power, height)
    return tuple(result)


def _kernel_problems(design: AlgebraDesign) -> list[str]:
    problems = []
    for point, partition in enumerate(design.groups):
        base = Fraction(point + 1)
        for factor_index, (factor, assignment) in enumerate(zip(design.factors, design.assignments)):
            owner = assignment[point]
            for group_index, block in enumerate(partition):
                if group_index == owner or not block:
                    continue
                rows = design.measurements[point][group_index]
                span = SparseRationalRowSpan(tuple(dict(row) for row in rows), len(block))
                position = {variable: local for local, variable in enumerate(block)}
                for variable, row in zip(factor.scope, factor.rows):
                    target = {}
                    for degree, term in enumerate(row.terms):
                        scale = base ** degree
                        for coordinate, coefficient in term.items():
                            local = position.get(coordinate)
                            if local is None:
                                continue
                            target[local] = target.get(local, Fraction(0)) + Fraction(coefficient) * scale
                    target = {local: value for local, value in target.items() if value}
                    if not span.contains(target):
                        problems.append(
                            f'point {point + 1} factor {factor_index} variable {variable} '
                            f'fails the kernel condition on block {group_index}'
                        )
    return problems


def _active_factors(system: SparseSystem, R: int, slope: Fraction):
    factors = []
    for equation in system.equations:
        if equation.is_zero:
            continue
        if equation.is_one:
            return None
        scope = tuple(sorted(equation.support))
        factors.append(_equation_factor(scope, equation.monomials, system.n, R, slope))
    for variable in range(system.n):
        factors.append(_neutral_factor(variable, system.n, R, slope))
    return factors


def _equation_factor(scope, monomials, n: int, R: int, slope: Fraction) -> EquationFactor:
    rows = tuple(_unit_row(variable, n, R, slope) for variable in scope)
    return EquationFactor('equation', scope, tuple(monomials), rows)


def _neutral_factor(variable: int, n: int, R: int, slope: Fraction) -> EquationFactor:
    return EquationFactor('neutral', (variable,), (), (_unit_row(variable, n, R, slope),))


def _unit_row(variable: int, n: int, R: int, slope: Fraction) -> PolynomialRow:
    del n
    terms = [{variable: Fraction(1)}]
    if slope != 0:
        if R < 1:
            raise ValueError('a nonzero slope needs R >= 1')
        terms.append({variable: Fraction(slope)})
    return PolynomialRow(tuple(terms))


def _accumulator_bit_bound(design: AlgebraDesign, parity: ParityPolynomials, plan: AlgebraPlan) -> int:
    """A planning bound. It does not allocate a jet or scan an assignment."""
    coefficient_bits = 1
    for factor in design.factors:
        for row in factor.rows:
            for degree, term in enumerate(row.terms):
                for coefficient in term.values():
                    size = abs(coefficient.numerator) * (design.Q ** degree)
                    bits = size.bit_length() + coefficient.denominator.bit_length()
                    if bits > coefficient_bits:
                        coefficient_bits = bits
    equations = max(design.equation_count, 1)
    return (
        (coefficient_bits + parity.coefficient_bits + design.Q.bit_length() + 8)
        * plan.J * parity.ell * equations * 4
        + design.n + 64
    )


def _parity_budgets(limits: AlgebraBudgets) -> Section5Budgets:
    return Section5Budgets(
        max_ell=limits.max_ell,
        max_coefficient_bits=limits.max_coefficient_bits,
    )


def _plan_fields(plan: AlgebraPlan) -> dict:
    return dict(
        input_degree_bound=plan.input_degree_bound,
        h_degree_bound=plan.h_degree_bound,
        q=plan.q,
        jet_length=plan.J,
        substitution_pairs=plan.substitution_pairs,
        algebra_slots=plan.algebra_slots,
        vector_bound=plan.boolean_vector_bound,
    )


def _slope_text(design: AlgebraDesign) -> str:
    for factor in design.factors:
        for row in factor.rows:
            if row.degree >= 1:
                for coefficient in row.terms[1].values():
                    return str(Fraction(coefficient))
    return '0'


def _substitution_count(n: int, q: int) -> int:
    if q >= n:
        return 3 ** n
    total = term = 1
    for size in range(1, q + 1):
        term = term * (n - size + 1) * 2 // size
        total += term
    return total


def _ceil_log2(n: int) -> int:
    if n <= 1:
        return 0
    return (n - 1).bit_length()


def _ceil_sqrt(n: int) -> int:
    root = math.isqrt(n)
    if root * root < n:
        root += 1
    return root


def _counted(count: int, ell: int, *, reason: str) -> ReplacementSumResult:
    return ReplacementSumResult(
        status='counted', count=count, reason=reason, ell=ell,
        enumerated=False, jets_allocated=False, coefficients_constructed=False,
    )


def _invalid(ell: int, reason: str, **extra) -> ReplacementSumResult:
    return ReplacementSumResult(
        status='invalid-design', count=None, reason=reason, ell=ell,
        boolean_vectors=0, enumerated=False, jets_allocated=False,
        coefficients_constructed=False, **extra,
    )


def _refused(ell: int, *, reason: str, heavy=(), branches=0, occurrence_cap=0,
             coefficients_constructed=False, coefficient_bits=0, **plan) -> ReplacementSumResult:
    del coefficient_bits
    return ReplacementSumResult(
        status='refused', count=None, reason=reason, ell=ell,
        boolean_vectors=0, branches=branches, branches_enumerated=0,
        heavy_variables=tuple(heavy), occurrence_cap=occurrence_cap,
        jets_allocated=False, coefficients_constructed=coefficients_constructed,
        enumerated=False, **plan,
    )


def _alg_zero(ell: int):
    return (Fraction(0),) * ell


def _alg_one(ell: int):
    return (Fraction(1),) + (Fraction(0),) * (ell - 1)


def _alg_add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def _alg_mul(left, right):
    ell = len(left)
    out = [Fraction(0)] * ell
    for i, left_value in enumerate(left):
        if not left_value:
            continue
        for j, right_value in enumerate(right):
            if i + j >= ell or not right_value:
                continue
            out[i + j] += left_value * right_value
    return tuple(out)


def _ajet_zero(ell: int, length: int):
    return (_alg_zero(ell),) * length


def _ajet_one(ell: int, length: int):
    return (_alg_one(ell),) + (_alg_zero(ell),) * (length - 1)


def _ajet_add(left, right):
    return tuple(_alg_add(a, b) for a, b in zip(left, right))


def _ajet_scale(jet, scalar: int):
    if scalar == 0:
        return _ajet_zero(len(jet[0]), len(jet))
    if scalar == 1:
        return jet
    if scalar == -1:
        return tuple(tuple(-coefficient for coefficient in element) for element in jet)
    return tuple(tuple(scalar * coefficient for coefficient in element) for element in jet)


def _ajet_mul(left, right):
    length = len(left)
    ell = len(left[0])
    out = [_alg_zero(ell) for _ in range(length)]
    for i, left_value in enumerate(left):
        if all(coefficient == 0 for coefficient in left_value):
            continue
        for j, right_value in enumerate(right):
            if i + j >= length or all(coefficient == 0 for coefficient in right_value):
                continue
            out[i + j] = _alg_add(out[i + j], _alg_mul(left_value, right_value))
    return tuple(out)


def _ajet_sum(element) -> Fraction:
    return sum(element, Fraction(0))
