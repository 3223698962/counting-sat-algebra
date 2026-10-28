"""Step 17 finite checks: Section 5.1–5.2 sparse systems. Not Theorem 2.2.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from .protocol import _row
from ..sparse_systems import (
    Section5Budgets,
    beigel_tarui_complement,
    brute_common_zeros,
    build_truncated_product,
    construct_parity_polynomials,
    count_common_zeros,
    evaluate_truncated_product,
    integer_lift,
    monomial_allowance,
    sparse_polynomial,
    sparse_system,
    truncated_indicator_residue,
)


def _height(y: int) -> int:
    return y * (y - 1)


def _paper_ell2_truncated(y1: int, y2: int) -> int:
    """The three-term product displayed after Lemma 5.2."""
    return (
        (1 - y1) * (1 - y2)
        + (2 * y1 - 1) * _height(y1) * (1 - y2)
        + (1 - y1) * (2 * y2 - 1) * _height(y2)
    )


def _paper_ell2_omitted(y1: int, y2: int) -> int:
    return (2 * y1 - 1) * (2 * y2 - 1) * _height(y1) * _height(y2)


def run_section5_suite() -> list[dict]:
    """Finite Section 5.1–5.2 checks. Not Theorem 2.2's running time."""
    print('=== (S) Section 5 sparse systems, Lemma 5.1, truncated product ===', flush=True)
    rows = []
    max_bits = 0
    identity_ok = True
    samples = range(-6, 13)
    for ell in range(1, 5):
        parity = construct_parity_polynomials(ell)
        again = construct_parity_polynomials(ell)
        raw = beigel_tarui_complement(ell)
        max_bits = max(max_bits, parity.coefficient_bits)
        identity_ok = (
            identity_ok
            and parity.coefficients == again.coefficients
            and raw[0] == 1
            and all(coefficient == 0 for coefficient in raw[1:ell])
            and all(parity.indicator_agrees(y) for y in samples)
        )
    rows.append(
        _row(
            'lemma-5.1-integers',
            0,
            0,
            identity_ok,
            'CHECKED',
            ell_values=[1, 2, 3, 4],
            integer_samples=len(samples),
            max_coefficient_bits=max_bits,
            count_checked=False,
            note='Even y is congruent to 1 and odd y to 0 modulo 2^ell. Finite check only.',
        )
    )
    print(f'  [S] Lemma 5.1 on ell=1..4: {"CHECKED" if identity_ok else "FAIL"}', flush=True)

    parity = construct_parity_polynomials(2)
    raw = beigel_tarui_complement(2)
    product = build_truncated_product(2, 2)
    ell2_ok = (
        parity.coefficients == ((1, -1), (-1, 2))
        and raw == (1, 0, -3, 2)
        and [list(term.powers) for term in product.terms] == [[0, 0], [0, 1], [1, 0]]
        and product.shape_ok()
        and all(term.linear_factors <= 2 and term.quadratic_factors <= 1 for term in product.terms)
        and (1, 1) not in {term.powers for term in product.terms}
    )
    for y1 in range(-4, 7):
        for y2 in range(-4, 7):
            got = evaluate_truncated_product(product, parity, (y1, y2))
            full = parity.evaluate(y1) * parity.evaluate(y2)
            if (
                got != _paper_ell2_truncated(y1, y2)
                or got + _paper_ell2_omitted(y1, y2) != full
                or _paper_ell2_omitted(y1, y2) % 4
                or truncated_indicator_residue(parity, (y1, y2)) != got % 4
            ):
                ell2_ok = False
    rows.append(
        _row(
            'lemma-5.1-ell-2-product',
            0,
            2,
            ell2_ok,
            'CHECKED',
            coefficients=[[1, -1], [-1, 2]],
            p_coefficients=[1, 0, -3, 2],
            retained_powers=[[0, 0], [0, 1], [1, 0]],
            omitted_powers=[1, 1],
            linear_factor_bound=2,
            quadratic_factor_bound=1,
            count_checked=False,
            note='The omitted h(y1)h(y2) term is divisible by 4 at integer inputs.',
        )
    )
    print(f'  [S] ell=2 displayed product: {"CHECKED" if ell2_ok else "FAIL"}', flush=True)

    cancelled = sparse_polynomial(4, [[0], [0], [1]], degree_bound=2, monomial_bound=3)
    constant = sparse_polynomial(2, [[]], degree_bound=1, monomial_bound=1)
    zero = sparse_polynomial(4, [], degree_bound=2, monomial_bound=2)
    collapsed = sparse_polynomial(4, [[0, 0, 1]], degree_bound=2, monomial_bound=1)
    lift_system = sparse_system(4, [[[0, 1], [2]]], k=2, c=1)
    pointwise_lift = True
    for mask in range(16):
        assignment = [(mask >> variable) & 1 for variable in range(4)]
        lift = integer_lift(lift_system.equations[0], assignment)
        x0, x1, x2 = assignment[:3]
        pointwise_lift = pointwise_lift and (lift % 2 == 0) == ((x0 * x1 + x2) % 2 == 0)
    rejected_long = False
    try:
        sparse_system(4, [[[0], [1], [2]]], k=1, c=1)
    except ValueError:
        rejected_long = True
    rejected_degree = False
    try:
        sparse_system(4, [[[0, 1, 2]]], k=2, c=1)
    except ValueError:
        rejected_degree = True
    normalization_ok = (
        cancelled.monomials == (frozenset({1}),)
        and constant.is_one
        and integer_lift(constant, [0, 1]) == 1
        and zero.is_zero
        and integer_lift(zero, [1, 1, 1, 1]) == 0
        and collapsed.monomials == (frozenset({0, 1}),)
        and collapsed.degree == 2
        and monomial_allowance(1, 1) == 0
        and monomial_allowance(3, 1) == 1
        and monomial_allowance(3, 2) == 2
        and monomial_allowance(4, 1) == 2
        and monomial_allowance(1024, 3) == 1000
        and lift_system.equations[0].degree == 2
        and pointwise_lift
        and rejected_long
        and rejected_degree
        and lift_system.meets_theorem_size is False
    )
    rows.append(
        _row(
            'sparse-f2-lift',
            4,
            1,
            normalization_ok,
            'CHECKED',
            monomial_allowance_n4_c1=2,
            monomial_allowance_n1024_c3=1000,
            duplicates_cancelled=True,
            empty_monomial_is_one=True,
            count_checked=False,
            note='(log2 n)^c is an exact integer comparison. n=4 is below the theorem size.',
        )
    )
    print(
        f'  [S] sparse normalization and lift parity: {"CHECKED" if normalization_ok else "FAIL"}',
        flush=True,
    )

    system = sparse_system(4, [[[0, 1], [2]], [[3], []]], k=2, c=1)
    counted = count_common_zeros(system)
    brute = brute_common_zeros(system)
    counted_parity = construct_parity_polynomials(5)
    counted_product = build_truncated_product(2, 5)
    pointwise = counted_product.shape_ok()
    for mask in range(16):
        assignment = [(mask >> variable) & 1 for variable in range(4)]
        lifts = [integer_lift(equation, assignment) for equation in system.equations]
        residue = evaluate_truncated_product(counted_product, counted_parity, lifts) % 32
        indicator = 0 if any(lift % 2 for lift in lifts) else 1
        pointwise = pointwise and residue == indicator
        pointwise = pointwise and truncated_indicator_residue(counted_parity, lifts) == indicator
    counted_ok = (
        counted.counted
        and counted.count == brute == 4
        and counted.boolean_vectors == 16
        and counted.terms == 15
        and counted.ell == 5
        and counted.coefficient_bits == 7
        and counted.input_degree_bound == 20
        and counted.linear_factor_bound == 2
        and counted.quadratic_factor_bound == 4
        and counted.term_table_allocated
        and counted.coefficients_constructed
        and counted.enumerated
        and pointwise
        and counted.count is not None
    )
    rows.append(
        _row(
            'tiny-common-zeros-n4',
            system.n,
            system.m_test,
            counted_ok,
            'COUNTED',
            brute=brute,
            count=counted.count,
            ell=counted.ell,
            terms=counted.terms,
            boolean_vectors=counted.boolean_vectors,
            coefficient_bits=counted.coefficient_bits,
            input_degree_bound=counted.input_degree_bound,
            linear_factor_bound=counted.linear_factor_bound,
            quadratic_factor_bound=counted.quadratic_factor_bound,
            pointwise_indicators=pointwise,
            term_table_allocated=True,
            note='Exhaustive residue on 16 points. Not Section 5.3, Section 6, or Theorem 1.1.',
        )
    )
    print(
        f'  [S] n=4 common zeros: N(P)={counted.count}; {"COUNTED" if counted_ok else "FAIL"}',
        flush=True,
    )

    empty = count_common_zeros(sparse_system(5, [], k=2, c=1))
    empty_ok = (
        empty.counted
        and empty.count == brute_common_zeros(sparse_system(5, [], k=2, c=1)) == 32
        and empty.boolean_vectors == 0
        and empty.term_table_allocated is False
        and empty.coefficients_constructed is False
        and empty.enumerated is False
    )
    rows.append(
        _row(
            'empty-system-n5',
            5,
            0,
            empty_ok,
            'COUNTED',
            brute=32,
            count=empty.count,
            boolean_vectors=0,
            term_table_allocated=False,
            note='An empty equation list is 2^n. No cube was scanned.',
        )
    )
    print(
        f'  [S] empty system: N(P)={empty.count}; {"COUNTED" if empty_ok else "FAIL"}', flush=True
    )

    constant_system = sparse_system(4, [[[]]], k=2, c=1)
    constant = count_common_zeros(constant_system)
    constant_ok = (
        constant_system.equations[0].is_one
        and constant.counted
        and constant.count == 0
        and constant.boolean_vectors == 0
        and constant.term_table_allocated is False
        and constant.enumerated is False
        and brute_common_zeros(constant_system) == 0
    )
    rows.append(
        _row(
            'constant-one-equation',
            4,
            1,
            constant_ok,
            'COUNTED',
            brute=0,
            count=constant.count,
            boolean_vectors=0,
            term_table_allocated=False,
            note='The constant equation 1 = 0 has no solutions. No partial search.',
        )
    )
    print(
        f'  [S] constant 1 = 0: N(P)={constant.count}; {"COUNTED" if constant_ok else "FAIL"}',
        flush=True,
    )

    term_refusal = count_common_zeros(system, Section5Budgets(max_terms=2))
    term_ok = (
        term_refusal.refused
        and term_refusal.count is None
        and term_refusal.boolean_vectors == 0
        and term_refusal.term_table_allocated is False
        and term_refusal.coefficients_constructed is False
        and 'term table not allocated' in term_refusal.reason
    )
    rows.append(
        _row(
            'term-table-budget',
            system.n,
            system.m_test,
            term_ok,
            'REFUSED',
            count=None,
            boolean_vectors=0,
            term_table_allocated=False,
            coefficients_constructed=False,
            term_budget=2,
            reason=term_refusal.reason,
        )
    )
    print(f'  [S] term budget: {"REFUSED before allocation" if term_ok else "FAIL"}', flush=True)

    bit_refusal = count_common_zeros(system, Section5Budgets(max_coefficient_bits=1))
    bit_ok = (
        bit_refusal.refused
        and bit_refusal.count is None
        and bit_refusal.boolean_vectors == 0
        and bit_refusal.term_table_allocated is False
        and bit_refusal.coefficients_constructed is False
        and 'not constructed' in bit_refusal.reason
    )
    rows.append(
        _row(
            'coefficient-bit-budget',
            system.n,
            system.m_test,
            bit_ok,
            'REFUSED',
            count=None,
            boolean_vectors=0,
            term_table_allocated=False,
            coefficients_constructed=False,
            coefficient_bit_budget=1,
            reason=bit_refusal.reason,
        )
    )
    print(
        f'  [S] coefficient bits: {"REFUSED before the term table" if bit_ok else "FAIL"}',
        flush=True,
    )

    cube_refusal = count_common_zeros(system, Section5Budgets(max_boolean_vectors=1))
    cube_ok = (
        cube_refusal.refused
        and cube_refusal.count is None
        and cube_refusal.boolean_vectors == 0
        and cube_refusal.term_table_allocated is False
        and cube_refusal.coefficients_constructed is False
        and cube_refusal.enumerated is False
        and 'no assignments scanned' in cube_refusal.reason
    )
    rows.append(
        _row(
            'cube-budget',
            system.n,
            system.m_test,
            cube_ok,
            'REFUSED',
            count=None,
            boolean_vectors=0,
            term_table_allocated=False,
            coefficients_constructed=False,
            boolean_vector_budget=1,
            reason=cube_refusal.reason,
        )
    )
    print(f'  [S] cube budget: {"REFUSED before the scan" if cube_ok else "FAIL"}', flush=True)
    return rows
