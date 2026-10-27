"""Step 3: tiny M0 = V0 f_I(W0), restoration, degree, and a legal design.

The side length is r**(2Q) and the slice stops at 64. These tests do not claim
Theorem 1.1, a subexponential runtime, or the matrix at R ≥ 32, 6R ≤ Q.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.design import popcount_design, validate_design
from counting_sat_algebra.replacement_matrix import (
    FULL_M0_LOG2_DIMENSION_AT_R2,
    FULL_M0_MIN_Q,
    FULL_M0_MIN_R,
    SLICE_MAX_DIMENSION,
    auxiliary_problems,
    binomial_exceptional_count,
    build_auxiliary_slice,
    default_interval_start,
    degree_cap_report,
    design_problems_at,
    groups_at,
    place_formula,
    restoration_problems,
    scalar_f_coefficients,
    step3_examples,
    vector_in_row_span,
)
from counting_sat_algebra.replacements import eval_H_at, eval_poly


def _example(name: str):
    found = [example for example in step3_examples() if example.name == name]
    assert len(found) == 1
    return found[0]


def _closed_form_remark(t: Fraction) -> Fraction:
    """Four assignments of (X0 - (T/2) X1)(X1 - (T/2) X0)."""
    total = Fraction(0)
    for x0 in (0, 1):
        for x1 in (0, 1):
            total += (x0 - (t / 2) * x1) * (x1 - (t / 2) * x0)
    return total


def test_scalar_interpolation_hits_alternating_signs() -> None:
    cases = [(0, 1), (0, 2), (0, 3), (1, 1), (1, 2)]
    for start, length in cases:
        coeffs = list(scalar_f_coefficients(start, length))
        for offset in range(length):
            point = start + offset
            assert eval_poly(coeffs, Fraction(point)) == Fraction((-1) ** point)


def test_vector_in_row_span_exact() -> None:
    assert vector_in_row_span([], [Fraction(0), Fraction(0)])
    assert not vector_in_row_span([], [Fraction(1)])
    assert vector_in_row_span([[Fraction(0), Fraction(1)]], [Fraction(0), Fraction(2)])
    assert not vector_in_row_span([[Fraction(1), Fraction(1)]], [Fraction(2), Fraction(3)])
    assert vector_in_row_span(
        [[Fraction(1), Fraction(1)], [Fraction(1), Fraction(1)]],
        [Fraction(2), Fraction(2)],
    )


def test_swap_matrix_is_the_4_by_4_exchange() -> None:
    # Order (a, d): (0,0), (0,1), (1,0), (1,1). f is the constant 1, so M0 = V0
    # exchanges the two off-diagonal coordinates and fixes the diagonal.
    sl = build_auxiliary_slice(2, 1, 1)
    assert (sl.r, sl.Q, sl.R, sl.dimension, sl.A0, sl.interval) == (2, 1, 1, 4, 0, (0,))
    assert sl.f_coefficients == (Fraction(1),)
    assert sl.M0 == (
        (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
    )
    assert sl.C_rows == ((Fraction(0), Fraction(1), Fraction(-1), Fraction(0)),)
    assert auxiliary_problems(sl) == []
    for u in range(4):
        for y in range(4):
            poly = sl.rows[u][y]
            assert poly == ((Fraction(1),) if u == y else (Fraction(0),))


def test_remark_hand_row_and_closed_form() -> None:
    # r = 2, Q = 2, R = 2. R = Q, outside the hypothesis R < Q.
    # M0 row of (00, 11) is +ζ_(10,01) + ζ_(01,10) - ζ_(11,00).
    # The replacement row is ζ_(00,11) - (T/2) ζ_(10,01) - T ζ_(01,10).
    sl = build_auxiliary_slice(2, 2, 2)
    assert sl.R == sl.Q == 2
    assert sl.dimension == 16
    assert sl.interval == (0, 1)
    assert sl.f_coefficients == (Fraction(1), Fraction(-2))
    assert len(sl.C_rows) == 1
    assert auxiliary_problems(sl) == []

    u = sl.index_of((0, 0), (1, 1))
    z1 = sl.index_of((1, 0), (0, 1))
    z2 = sl.index_of((0, 1), (1, 0))
    z12 = sl.index_of((1, 1), (0, 0))
    expected_m = [Fraction(0) for _ in range(sl.dimension)]
    expected_m[z1] = Fraction(1)
    expected_m[z2] = Fraction(1)
    expected_m[z12] = Fraction(-1)
    assert list(sl.M0[u]) == expected_m

    for y in range(sl.dimension):
        poly = sl.rows[u][y]
        if y == u:
            assert poly == (Fraction(1),)
        elif y == z1:
            assert poly == (Fraction(0), Fraction(-1, 2))
        elif y == z2:
            assert poly == (Fraction(0), Fraction(-1))
        else:
            assert poly == (Fraction(0),)

    measurement = [Fraction(0) for _ in range(sl.dimension)]
    measurement[u] = Fraction(1)
    measurement[z1] = Fraction(-1)
    measurement[z2] = Fraction(-1)
    measurement[z12] = Fraction(1)
    assert list(sl.C_rows[0]) == measurement

    # At t = 1 and t = 2 the same row equals the measurement plus assigned-group terms.
    # ⟨P(1), ζ⟩ = c(ζ) + (1/2) ζ1 - ζ12, ⟨P(2), ζ⟩ = c(ζ) - ζ2 - ζ12.
    def coeff(index: int, t: Fraction) -> Fraction:
        return eval_poly(list(sl.rows[u][index]), t)

    for t, extra in (
        (Fraction(1), {z1: Fraction(1, 2), z12: Fraction(-1)}),
        (Fraction(2), {z2: Fraction(-1), z12: Fraction(-1)}),
    ):
        for y in range(sl.dimension):
            base = measurement[y]
            got = coeff(y, t)
            assert got == base + extra.get(y, Fraction(0))

    example = _example("remark-4-4")
    placed = example.placed
    assert restoration_problems(placed) == []
    assert placed.family.rows[0][0][0] == [Fraction(1)]
    assert placed.family.rows[0][0][1] == [Fraction(0), Fraction(-1, 2)]
    assert placed.family.rows[1][0][0] == [Fraction(0), Fraction(-1, 2)]
    assert placed.family.rows[1][0][1] == [Fraction(1)]
    assert placed.measurement_rows == ((Fraction(1), Fraction(-1)),)
    assert groups_at(placed, 1) == ((0,), (1,))
    assert groups_at(placed, 2) == ((0, 1), ())
    assert design_problems_at(placed, 1) == []
    assert design_problems_at(placed, 2) == []

    for t in (Fraction(0), Fraction(1), Fraction(2), Fraction(1, 2)):
        closed = Fraction(1) - 2 * t + (t * t) / 4
        assert _closed_form_remark(t) == closed
        assert eval_H_at(placed.cnf, placed.family, t) == closed
    report = degree_cap_report(placed)
    assert report.ok
    assert report.brute == 1
    assert report.observed_degree == 2
    assert report.cap_occurrences == 2
    assert report.cap_step1 == 4
    assert report.observed_degree == report.cap_occurrences


def test_examples_match_brute_force_degree_and_design() -> None:
    expected = {
        "swap-4": {
            "D": 4,
            "rows": 1,
            "degree": 0,
            "cap_m": 0,
            "cap_step1": 2,
            "brute": 1,
            "values": {Fraction(0): Fraction(1), Fraction(1): Fraction(1), Fraction(2): Fraction(1)},
        },
        "remark-4-4": {
            "D": 16,
            "rows": 1,
            "degree": 2,
            "cap_m": 2,
            "cap_step1": 4,
            "brute": 1,
            "values": {
                Fraction(0): Fraction(1),
                Fraction(1): Fraction(-3, 4),
                Fraction(2): Fraction(-2),
                Fraction(1, 2): Fraction(1, 16),
            },
        },
        "cross-r3": {
            "D": 9,
            "rows": 3,
            "degree": 0,
            "cap_m": 0,
            "cap_step1": 3,
            "brute": 7,
            "values": {Fraction(0): Fraction(7), Fraction(1): Fraction(7), Fraction(2): Fraction(7)},
        },
        "r2-q3": {
            "D": 64,
            "rows": 10,
            "degree": 2,
            "cap_m": 3,
            "cap_step1": 6,
            "brute": 7,
            "values": {Fraction(0): Fraction(7), Fraction(1): Fraction(8), Fraction(2): Fraction(29, 3)},
        },
    }
    examples = step3_examples()
    assert [example.name for example in examples] == list(expected)
    assert any(example.auxiliary.R < example.auxiliary.Q for example in examples)
    assert any(example.auxiliary.R == example.auxiliary.Q for example in examples)
    for example in examples:
        spec = expected[example.name]
        sl = example.auxiliary
        placed = example.placed
        assert sl.dimension == spec["D"]
        assert len(sl.C_rows) == spec["rows"]
        assert len(sl.C_rows) == binomial_exceptional_count(sl.r, sl.Q, sl.interval)
        assert auxiliary_problems(sl) == []
        assert restoration_problems(placed) == []
        assert brute_count(placed.cnf) == spec["brute"]
        for t in range(1, sl.Q + 1):
            assert design_problems_at(placed, t) == []
            sizes = [len(group) for group in groups_at(placed, t)]
            assert all(size * sl.r <= 2 * placed.cnf.n for size in sizes)
        report = degree_cap_report(placed)
        assert report.ok
        assert report.observed_degree == spec["degree"]
        assert report.cap_occurrences == spec["cap_m"]
        assert report.cap_step1 == spec["cap_step1"]
        assert report.direct_at_zero == Fraction(spec["brute"])
        for t in example.sample_points:
            assert eval_H_at(placed.cnf, placed.family, t) == spec["values"][t]
        assignment = tuple(1 for _ in range(placed.cnf.n))
        for clause_index, clause in enumerate(placed.cnf.clauses):
            for occ, (var, _neg) in enumerate(clause):
                value = placed.family.linear_form(clause_index, occ, assignment, Fraction(0))
                assert value == Fraction(assignment[var])


def test_cross_r3_meets_three_groups_and_step2_rejects_it() -> None:
    example = _example("cross-r3")
    placed = example.placed
    assert placed.cnf.n == 3
    assert groups_at(placed, 1) == ((0,), (1,), (2,))
    assert placed.clause_schedules == ((0,),)
    # The clause uses variables 1 and 2, which are outside its assigned group.
    outside = []
    for var, _neg in placed.cnf.clauses[0]:
        if var not in groups_at(placed, 1)[placed.clause_schedules[0][0]]:
            outside.append(var)
    assert outside == [1, 2]
    assert design_problems_at(placed, 1) == []
    assert brute_count(placed.cnf) == 7
    for u, label in enumerate(example.auxiliary.labels):
        swapped = example.auxiliary.index_of(label[1], label[0])
        row = example.auxiliary.M0[u]
        assert row[swapped] == 1
        assert sum(1 for entry in row if entry != 0) == 1

    blocked = popcount_design(3, [(0,), (1,), (2,)], D=0, r=1)
    with pytest.raises(ValueError, match="one block"):
        validate_design(blocked, placed.cnf)


def test_full_interval_makes_M0_the_identity() -> None:
    # Default start is -1, which is rejected. The explicit interval {0, 1, 2}
    # fixes every basis direction, so there are no measurement rows and M0 = I.
    assert default_interval_start(2, 2, 3) == -1
    with pytest.raises(ValueError, match="not contained"):
        build_auxiliary_slice(2, 2, 3)
    sl = build_auxiliary_slice(2, 2, 3, A0=0)
    assert sl.C_rows == ()
    assert sl.interval == (0, 1, 2)
    assert auxiliary_problems(sl) == []
    for u in range(sl.dimension):
        for y in range(sl.dimension):
            assert sl.M0[u][y] == (Fraction(1) if u == y else Fraction(0))
    u = sl.index_of((0, 0), (1, 0))
    # One mismatch, at position s = 1, so the diagonal entry is 1 - T.
    assert sl.rows[u][u] == (Fraction(1), Fraction(-1))
    assert all(sl.rows[u][y] == (Fraction(0),) for y in range(sl.dimension) if y != u)


def test_side_length_above_the_cap_names_the_obstacle() -> None:
    assert FULL_M0_MIN_R == 32
    assert FULL_M0_MIN_Q == 192
    assert FULL_M0_LOG2_DIMENSION_AT_R2 == 384
    assert SLICE_MAX_DIMENSION == 64
    assert 2**FULL_M0_LOG2_DIMENSION_AT_R2 > SLICE_MAX_DIMENSION
    with pytest.raises(ValueError, match=r"2\*\*384"):
        build_auxiliary_slice(2, 4, 1)
    with pytest.raises(ValueError, match="cannot sit inside"):
        build_auxiliary_slice(2, 1, 3)
    with pytest.raises(ValueError, match="r must"):
        build_auxiliary_slice(1, 1, 1)


def test_repeated_label_uses_a_second_copy() -> None:
    cnf = CNF(1, [[(0, False)], [(0, False)]])
    sl = build_auxiliary_slice(2, 1, 1)
    placed = place_formula(cnf, sl, variable_schedules=((0,),), clause_schedules=((0,), (0,)))
    assert placed.copies == 2
    assert len(placed.measurement_rows) == 2
    assert restoration_problems(placed) == []
    assert design_problems_at(placed, 1) == []
    report = degree_cap_report(placed)
    assert report.ok
    assert report.brute == 1
    assert report.observed_degree == 0


def test_bad_schedules_are_rejected() -> None:
    sl = build_auxiliary_slice(2, 1, 1)
    cnf = CNF(1, [[(0, False)]])
    with pytest.raises(ValueError, match="length"):
        place_formula(cnf, sl, variable_schedules=((0, 1),), clause_schedules=((0,),))
    with pytest.raises(ValueError, match="outside"):
        place_formula(cnf, sl, variable_schedules=((2,),), clause_schedules=((0,),))
    with pytest.raises(ValueError, match="True"):
        build_auxiliary_slice(True, 1, 1)  # type: ignore[arg-type]
