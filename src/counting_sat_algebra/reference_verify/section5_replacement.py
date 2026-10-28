"""Step 18 finite checks: the Section 5.3–5.5 replacement sum in A_ell.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from .protocol import _row
from ..sparse_systems import (
    brute_common_zeros,
    construct_parity_polynomials,
    integer_lift,
    sparse_system,
)
from ..algebra_count import (
    AlgebraBudgets,
    count_replacement_sum,
    h_degree_bound,
    incidence_cut,
    input_degree_bound,
)


def _boolean_coefficient_sum(system, ell: int) -> tuple[int, ...]:
    """Cube sum of the truncated factors in Z[V]/(V^ell). Not the Section 3 routine."""
    parity = construct_parity_polynomials(ell)
    active = [equation for equation in system.equations if not equation.is_zero]
    total = [0] * ell
    for mask in range(1 << system.n):
        bits = [(mask >> variable) & 1 for variable in range(system.n)]
        accumulator = [0] * ell
        accumulator[0] = 1
        for equation in active:
            lift = integer_lift(equation, bits)
            height = lift * (lift - 1)
            column = []
            height_power = 1
            for constant, slope in parity.coefficients:
                column.append((constant + slope * lift) * height_power)
                height_power *= height
            combined = [0] * ell
            for left, left_value in enumerate(accumulator):
                if left_value == 0:
                    continue
                for right, right_value in enumerate(column):
                    if left + right < ell and right_value:
                        combined[left + right] += left_value * right_value
            accumulator = combined
        for slot, value in enumerate(accumulator):
            total[slot] += value
    return tuple(total)


def run_section5_replacement_suite() -> list[dict]:
    """Finite Sections 5.3–5.5. Not Section 6, Appendix D, or Theorem 2.2."""
    print('=== (T) Section 5.3–5.5 replacement sum in A_ell ===', flush=True)
    rows = []
    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)
    cut = incidence_cut(system)
    degree = input_degree_bound(system.k, 2, 5)
    replaced = h_degree_bound(system.k, 2, 5, 1)
    bounds_ok = (
        degree == 20
        and replaced == 20
        and cut.incidence_total == 8
        and cut.occurrence_cap == 4
        and cut.heavy == ()
        and cut.neutral_scopes == ((0,), (1,), (2,), (3,))
        and cut.equation_scopes == ((0, 1, 2), (3,))
        and cut.log_n == 2
        and cut.sqrt_log_n == 2
    )
    rows.append(
        _row(
            'replacement-degree-and-incidence',
            system.n,
            system.m_test,
            bounds_ok,
            'CHECKED',
            input_degree_bound=degree,
            h_degree_bound=replaced,
            incidence_total=cut.incidence_total,
            occurrence_cap=cut.occurrence_cap,
            heavy=[],
            neutral_factors=4,
            section_6_implemented=False,
            appendix_d_implemented=False,
            count_checked=False,
            note='Equation (5.5) and the Section 5.5 occurrence cap. No sum is run in this row.',
        )
    )
    print(
        f'  [T] degree d={degree}, d_H={replaced}, b0={cut.occurrence_cap}: {"CHECKED" if bounds_ok else "FAIL"}',
        flush=True,
    )

    counted = count_replacement_sum(system)
    oracle = _boolean_coefficient_sum(system, 5)
    brute = brute_common_zeros(system)
    counted_ok = (
        counted.counted
        and counted.count == brute == 4
        and counted.h0_coefficients == oracle
        and counted.unreduced_sum == sum(oracle) == 2308
        and counted.unreduced_sum % 32 == 4
        and counted.boolean_vectors == 399
        and counted.vector_bound == 714
        and counted.q == 0
        and counted.jet_length == 1
        and counted.design_Q == 21
        and counted.replacement_degree == 1
        and counted.slope == '1/21'
        and counted.input_degree_bound == 20
        and counted.h_degree_bound == 20
        and counted.neutral_factors == 4
        and counted.equations_used == 2
        and counted.branches == 1
        and counted.branch_counts == (4,)
        and counted.sample_at_one != '4'
        and counted.enumerated
        and counted.jets_allocated
        and counted.lemma_48_guard_passed is False
        and counted.lemma_48_dimension_formed is False
        and counted.lemma_48_r == 1
        and counted.lemma_48_Q == 0
        and counted.lemma_48_R == 10
    )
    rows.append(
        _row(
            'replacement-sum-n4',
            system.n,
            system.m_test,
            counted_ok,
            'COUNTED',
            brute=brute,
            count=counted.count,
            ell=counted.ell,
            h0_coefficients=list(counted.h0_coefficients),
            unreduced_sum=counted.unreduced_sum,
            boolean_vectors=counted.boolean_vectors,
            vector_bound=counted.vector_bound,
            factor_evaluations=counted.factor_evaluations,
            q=counted.q,
            jet_length=counted.jet_length,
            design_Q=counted.design_Q,
            replacement_degree=counted.replacement_degree,
            slope=counted.slope,
            input_degree_bound=counted.input_degree_bound,
            h_degree_bound=counted.h_degree_bound,
            sample_at_one=counted.sample_at_one,
            neutral_factors=counted.neutral_factors,
            note='S_ell(H(0)) residue 4 matches the cube oracle and the F2 scan. The row is not Lemma 4.8.',
        )
    )
    print(
        f'  [T] n=4 replacement sum: N(P)={counted.count}; vectors={counted.boolean_vectors}; '
        f'{"COUNTED" if counted_ok else "FAIL"}',
        flush=True,
    )

    branched_system = sparse_system(4, [[[0]], [[0]], [[0]], [[0]], [[0]], [[0], [1]]], k=1, c=1)
    branched = count_replacement_sum(branched_system)
    branched_brute = brute_common_zeros(branched_system)
    branched_ok = (
        branched.counted
        and branched.count == branched_brute == 4
        and branched.heavy_variables == (0,)
        and branched.occurrence_cap == 6
        and branched.branches == 2
        and branched.branch_counts == (4, 0)
        and branched.boolean_vectors == 110
        and branched.branches_enumerated == 1
        and branched.enumerated
        and branched.jets_allocated
    )
    rows.append(
        _row(
            'high-incidence-branches-n4',
            branched_system.n,
            branched_system.m_test,
            branched_ok,
            'COUNTED',
            brute=branched_brute,
            count=branched.count,
            branch_counts=list(branched.branch_counts),
            heavy=list(branched.heavy_variables),
            occurrence_cap=branched.occurrence_cap,
            boolean_vectors=branched.boolean_vectors,
            branches_enumerated=branched.branches_enumerated,
            note='Fixing the heavy variable splits the cube. The branches add to N(P).',
        )
    )
    print(
        f'  [T] high-incidence branches: N(P)={branched.count}; counts={list(branched.branch_counts)}; '
        f'{"COUNTED" if branched_ok else "FAIL"}',
        flush=True,
    )

    guard_ok = (
        branched.lemma_48_guard_passed is False
        and branched.lemma_48_dimension_formed is False
        and branched.lemma_48_r == 1
        and branched.lemma_48_Q == 0
        and branched.lemma_48_R == 8
        and branched.counted
    )
    rows.append(
        _row(
            'lemma-4.8-guard-on-branch',
            3,
            1,
            guard_ok,
            'CHECKED',
            r=branched.lemma_48_r,
            Q=branched.lemma_48_Q,
            R=branched.lemma_48_R,
            guard_passed=False,
            dimension_formed=False,
            count_checked=False,
            note="The restricted branch has n'=3. Tau and D are not formed. The finite design still summed.",
        )
    )
    print(
        f'  [T] Lemma 4.8 on the branch: r={branched.lemma_48_r}, Q={branched.lemma_48_Q}, '
        f'R={branched.lemma_48_R}; {"CHECKED" if guard_ok else "FAIL"}',
        flush=True,
    )

    vector_refusal = count_replacement_sum(system, AlgebraBudgets(max_boolean_vectors=1))
    vector_ok = (
        vector_refusal.refused
        and vector_refusal.count is None
        and vector_refusal.boolean_vectors == 0
        and vector_refusal.jets_allocated is False
        and vector_refusal.coefficients_constructed is False
        and vector_refusal.enumerated is False
        and vector_refusal.vector_bound == 714
        and 'no enumeration started' in vector_refusal.reason
    )
    rows.append(
        _row(
            'replacement-vector-budget',
            system.n,
            system.m_test,
            vector_ok,
            'REFUSED',
            count=None,
            boolean_vectors=0,
            vector_bound=vector_refusal.vector_bound,
            jets_allocated=False,
            coefficients_constructed=False,
            boolean_vector_budget=1,
            reason=vector_refusal.reason,
        )
    )
    print(
        f'  [T] vector budget: {"REFUSED before enumeration" if vector_ok else "FAIL"}', flush=True
    )

    branch_refusal = count_replacement_sum(branched_system, AlgebraBudgets(max_branches=1))
    branch_ok = (
        branch_refusal.refused
        and branch_refusal.count is None
        and branch_refusal.boolean_vectors == 0
        and branch_refusal.branches == 2
        and branch_refusal.branches_enumerated == 0
        and branch_refusal.enumerated is False
        and branch_refusal.jets_allocated is False
        and branch_refusal.coefficients_constructed is False
        and 'no branch enumerated' in branch_refusal.reason
    )
    rows.append(
        _row(
            'replacement-branch-budget',
            branched_system.n,
            branched_system.m_test,
            branch_ok,
            'REFUSED',
            count=None,
            boolean_vectors=0,
            branches=2,
            branches_enumerated=0,
            jets_allocated=False,
            branch_budget=1,
            reason=branch_refusal.reason,
        )
    )
    print(
        f'  [T] branch budget: {"REFUSED before enumeration" if branch_ok else "FAIL"}', flush=True
    )

    bit_refusal = count_replacement_sum(system, AlgebraBudgets(max_coefficient_bits=1))
    bit_ok = (
        bit_refusal.refused
        and bit_refusal.count is None
        and bit_refusal.boolean_vectors == 0
        and bit_refusal.coefficients_constructed is False
        and bit_refusal.jets_allocated is False
        and bit_refusal.enumerated is False
        and 'not constructed' in bit_refusal.reason
    )
    rows.append(
        _row(
            'replacement-coefficient-budget',
            system.n,
            system.m_test,
            bit_ok,
            'REFUSED',
            count=None,
            boolean_vectors=0,
            jets_allocated=False,
            coefficients_constructed=False,
            coefficient_bit_budget=1,
            reason=bit_refusal.reason,
        )
    )
    print(f'  [T] coefficient bits: {"REFUSED before the jets" if bit_ok else "FAIL"}', flush=True)
    return rows
