"""Step 4: Hermite / Hasse recovery of H(0) on tiny formulas.

The condition count stays D+1. Distinct sample points are fewer than the
Lagrange path. Each point still enumerates 2^n assignments. These tests do
not claim Theorem 1.1, a subexponential runtime, or 2^{o(n)}.
"""

from __future__ import annotations

import random
from fractions import Fraction

import pytest

from counting_sat_algebra.bounded_occurrence import equispaced_degree
from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.hermite import (
    analyze_step4_case,
    contradiction_case,
    degree_from_hasse,
    eval_H_hasse_jet,
    hasse_jet_value,
    hermite_eval,
    hermite_multiplicities,
    hermite_schemes,
    negative_unit_case,
    polynomial_hasse_jet,
    recovery_problems,
    step4_cases,
    taylor_from_hasse,
    unit_pair_case,
)
from counting_sat_algebra.replacements import (
    degree_bound,
    eval_H_at,
    eval_poly,
    interpolate_H0,
    lagrange_eval,
    lagrange_eval_at_zero,
    random_replacements,
)


def _closed_unit_pair(t: Fraction) -> Fraction:
    return Fraction(1) + 6 * t + 6 * t * t


def _closed_negative_unit(t: Fraction) -> Fraction:
    return Fraction(1) - t


def _closed_contradiction(t: Fraction) -> Fraction:
    return -t - t * t


def test_polynomial_hasse_jet_of_the_unit_pair_quadratic() -> None:
    coeffs = [Fraction(1), Fraction(6), Fraction(6)]
    jet = polynomial_hasse_jet(coeffs, Fraction(1), 3)
    # 1 + 6(1+u) + 6(1+u)^2 = 13 + 18u + 6u^2.
    assert jet == (Fraction(13), Fraction(18), Fraction(6), Fraction(0))
    assert taylor_from_hasse(jet[:3]) == (
        Fraction(13),
        Fraction(18),
        Fraction(12),
    )
    assert degree_from_hasse(jet) == 2
    assert hasse_jet_value(jet, Fraction(1), Fraction(0)) == Fraction(1)
    assert polynomial_hasse_jet(coeffs, Fraction(1), 4)[:4] == jet


def test_hermite_eval_reproduces_a_cubic_from_mixed_nodes() -> None:
    # (1/2) T^3 - T. Four conditions: multiplicity 2 at t = 1 and at t = 2.
    coeffs = [Fraction(0), Fraction(-1), Fraction(0), Fraction(1, 2)]
    samples = [
        (Fraction(1), polynomial_hasse_jet(coeffs, Fraction(1), 1)),
        (Fraction(2), polynomial_hasse_jet(coeffs, Fraction(2), 1)),
    ]
    assert hermite_eval(samples, Fraction(0)) == Fraction(0)
    assert hermite_eval(samples, Fraction(2)) == eval_poly(coeffs, Fraction(2))
    assert hermite_eval(samples, Fraction(5)) == eval_poly(coeffs, Fraction(5))

    single = [(Fraction(0), polynomial_hasse_jet(coeffs, Fraction(0), 3))]
    assert hermite_eval(single, Fraction(0)) == Fraction(0)
    assert hermite_eval(single, Fraction(4)) == eval_poly(coeffs, Fraction(4))


def test_value_only_nodes_match_lagrange() -> None:
    xs = [Fraction(t) for t in (1, 2, 4, Fraction(1, 2))]
    ys = [Fraction(t * t - 3) for t in xs]
    samples = [(x, (y,)) for x, y in zip(xs, ys)]
    assert hermite_eval(samples, Fraction(0)) == lagrange_eval_at_zero(xs, ys)
    assert hermite_eval(samples, Fraction(7)) == lagrange_eval(xs, ys, Fraction(7))


def test_multiplicities_pack_d_plus_one_conditions_into_fewer_points() -> None:
    assert hermite_multiplicities(3, 2) == (2, 2)
    assert hermite_multiplicities(3, 3) == (3, 1)
    assert hermite_multiplicities(3, 4) == (4,)
    assert hermite_multiplicities(12, 2) == (2, 2, 2, 2, 2, 2, 1)
    assert hermite_multiplicities(12, 3) == (3, 3, 3, 3, 1)
    assert hermite_multiplicities(0, 2) == (1,)
    assert hermite_schemes(1) == (("single-point", 2),)
    assert hermite_schemes(2) == (("paired", 2), ("single-point", 3))
    assert hermite_schemes(3) == (
        ("paired", 2),
        ("triple", 3),
        ("single-point", 4),
    )
    for degree_cap, orders in ((1, 2), (5, 2), (5, 3), (8, 4)):
        multiplicities = hermite_multiplicities(degree_cap, orders)
        assert sum(multiplicities) == degree_cap + 1
        assert len(multiplicities) < degree_cap + 1
        assert max(multiplicities) >= 2


def test_rejects_bad_hermite_input() -> None:
    with pytest.raises(ValueError):
        hermite_multiplicities(1, 0)
    with pytest.raises(ValueError):
        hermite_multiplicities(-1, 2)
    with pytest.raises(ValueError):
        polynomial_hasse_jet([Fraction(1)], Fraction(0), -1)
    with pytest.raises(ValueError):
        hermite_eval([], Fraction(0))
    with pytest.raises(ValueError):
        hermite_eval(
            [(Fraction(1), (Fraction(1),)), (Fraction(1), (Fraction(2),))],
            Fraction(0),
        )
    with pytest.raises(ValueError):
        eval_H_hasse_jet(CNF(2, []), negative_unit_case().family, Fraction(1), 0)


def test_negative_unit_closed_form_jet() -> None:
    case = negative_unit_case()
    assert brute_count(case.cnf) == 1
    for point in (Fraction(0), Fraction(1), Fraction(2), Fraction(1, 2)):
        assert eval_H_at(case.cnf, case.family, point) == _closed_negative_unit(point)
    jet = eval_H_hasse_jet(case.cnf, case.family, Fraction(1), 2)
    # 1 - (1+u) = -u.
    assert jet == (Fraction(0), Fraction(-1), Fraction(0))
    assert taylor_from_hasse(jet) == (Fraction(0), Fraction(-1), Fraction(0))
    assert degree_from_hasse(jet) == 1


def test_unit_pair_closed_form_and_centered_difference() -> None:
    case = unit_pair_case()
    assert brute_count(case.cnf) == 1
    for point in (Fraction(0), Fraction(1), Fraction(2), Fraction(1, 2)):
        assert eval_H_at(case.cnf, case.family, point) == _closed_unit_pair(point)
    jet = eval_H_hasse_jet(case.cnf, case.family, Fraction(1), 3)
    assert jet == (Fraction(13), Fraction(18), Fraction(6), Fraction(0))
    assert taylor_from_hasse(jet[:3]) == (
        Fraction(13),
        Fraction(18),
        Fraction(12),
    )
    at_zero = eval_H_hasse_jet(case.cnf, case.family, Fraction(0), 2)
    assert at_zero == (Fraction(1), Fraction(6), Fraction(6))
    # Degree ≤ 2: the first Hasse coefficient is the centered difference.
    center = Fraction(1)
    step = Fraction(1)
    centered = (
        eval_H_at(case.cnf, case.family, center + step)
        - eval_H_at(case.cnf, case.family, center - step)
    ) / (2 * step)
    assert jet[1] == centered == Fraction(18)


def test_contradiction_is_a_nonzero_polynomial_with_zero_count() -> None:
    case = contradiction_case()
    assert brute_count(case.cnf) == 0
    for point in (Fraction(0), Fraction(1), Fraction(2), Fraction(3, 2)):
        assert eval_H_at(case.cnf, case.family, point) == _closed_contradiction(point)
    jet = eval_H_hasse_jet(case.cnf, case.family, Fraction(1), 2)
    assert jet == (Fraction(-2), Fraction(-3), Fraction(-1))
    assert degree_from_hasse(jet) == 2
    assert jet != (Fraction(0), Fraction(0), Fraction(0))


def test_step4_cases_match_brute_force_with_fewer_points() -> None:
    seen: set[str] = set()
    for case in step4_cases():
        assert case.name not in seen
        seen.add(case.name)
        degree_cap = degree_bound(case.cnf, R=case.family.R, k=case.k)
        full_jet = eval_H_hasse_jet(
            case.cnf, case.family, Fraction(1), degree_cap
        )
        other = eval_H_hasse_jet(case.cnf, case.family, Fraction(2), degree_cap)
        assert degree_from_hasse(full_jet) == degree_from_hasse(other)
        assert len(full_jet) == degree_cap + 1
        prefix = eval_H_hasse_jet(case.cnf, case.family, Fraction(1), 1)
        assert prefix == full_jet[:2]
        for t in range(0, degree_cap + 2):
            predicted = hasse_jet_value(full_jet, Fraction(1), Fraction(t))
            assert predicted == eval_H_at(case.cnf, case.family, Fraction(t))
        window = [
            eval_H_at(case.cnf, case.family, Fraction(t))
            for t in range(1, degree_cap + 3)
        ]
        assert degree_from_hasse(full_jet) == equispaced_degree(window)

        rows = analyze_step4_case(case)
        assert rows
        schemes = [row.scheme for row in rows]
        assert schemes[0] in {"paired", "single-point"}
        assert len(schemes) == len(set(schemes))
        for row in rows:
            problems = recovery_problems(row)
            assert problems == [], (case.name, row.scheme, problems)
            assert row.h0 == Fraction(brute_count(case.cnf))
            lagrange_h0, lagrange_cap, _fam = interpolate_H0(
                case.cnf, fam=case.family, k=case.k
            )
            assert lagrange_cap == row.degree_cap
            assert row.lagrange_h0 == lagrange_h0 == row.h0
            assert row.conditions == row.degree_cap + 1
            assert row.distinct_points < row.lagrange_points
            assert row.cube_passes == row.distinct_points
            assert row.assignments == row.distinct_points * (2**row.n)
            assert max(row.multiplicities) == row.max_derivative_order + 1
            assert row.observed_degree <= row.degree_cap
            assert row.observed_degree == degree_from_hasse(full_jet)

    assert seen == {
        "negative-unit",
        "unit-pair",
        "contradiction-n1",
        "one-clause-n3",
        "random-n4-m4-r1",
        "random-n4-m3-r2",
        "random-n6-m6-r1",
        "bounded-n6-m4-r1",
    }


def test_recorded_point_counts_and_known_counts() -> None:
    by_name = {case.name: case for case in step4_cases()}

    negative = analyze_step4_case(by_name["negative-unit"])
    assert len(negative) == 1
    assert negative[0].scheme == "single-point"
    assert negative[0].degree_cap == 1
    assert negative[0].multiplicities == (2,)
    assert negative[0].distinct_points == 1
    assert negative[0].lagrange_points == 2
    assert negative[0].brute == 1
    assert negative[0].observed_degree == 1
    assert negative[0].h0 == Fraction(1)

    pair = {row.scheme: row for row in analyze_step4_case(by_name["unit-pair"])}
    assert set(pair) == {"paired", "single-point"}
    assert pair["paired"].multiplicities == (2, 1)
    assert pair["paired"].distinct_points == 2
    assert pair["single-point"].multiplicities == (3,)
    assert pair["single-point"].distinct_points == 1
    assert pair["paired"].lagrange_points == 3
    assert pair["paired"].brute == 1
    assert pair["paired"].observed_degree == 2
    assert pair["paired"].h0 == Fraction(1)

    contradiction = {
        row.scheme: row for row in analyze_step4_case(by_name["contradiction-n1"])
    }
    assert contradiction["paired"].h0 == Fraction(0)
    assert contradiction["paired"].brute == 0
    assert contradiction["paired"].observed_degree == 2
    assert contradiction["paired"].multiplicities == (2, 1)
    assert contradiction["single-point"].distinct_points == 1

    clause = {
        row.scheme: row for row in analyze_step4_case(by_name["one-clause-n3"])
    }
    assert clause["paired"].degree_cap == 3
    assert clause["paired"].multiplicities == (2, 2)
    assert clause["paired"].distinct_points == 2
    assert clause["triple"].multiplicities == (3, 1)
    assert clause["triple"].distinct_points == 2
    assert clause["single-point"].multiplicities == (4,)
    assert clause["single-point"].distinct_points == 1
    assert clause["paired"].lagrange_points == 4
    assert clause["paired"].brute == 7
    assert clause["paired"].h0 == Fraction(7)

    bounded = {
        row.scheme: row
        for row in analyze_step4_case(by_name["bounded-n6-m4-r1"])
    }
    assert bounded["paired"].degree_cap == 12
    assert bounded["paired"].multiplicities == (2, 2, 2, 2, 2, 2, 1)
    assert bounded["paired"].distinct_points == 7
    assert bounded["triple"].multiplicities == (3, 3, 3, 3, 1)
    assert bounded["triple"].distinct_points == 5
    assert bounded["single-point"].distinct_points == 1
    assert bounded["single-point"].multiplicities == (13,)
    assert bounded["paired"].lagrange_points == 13
    assert bounded["paired"].brute == 44
    assert bounded["paired"].h0 == Fraction(44)
    assert bounded["paired"].observed_degree == 12
    assert bounded["paired"].assignments == 7 * 64
    assert bounded["single-point"].assignments == 64


def test_one_clause_matches_step0_lagrange_family() -> None:
    case = next(item for item in step4_cases() if item.name == "one-clause-n3")
    h0, degree_cap, _fam = interpolate_H0(
        case.cnf, R=1, rng=random.Random(0), k=3
    )
    # Same formula and the same replacement seed as the Step 0 single-clause test.
    assert degree_cap == 3
    assert h0 == Fraction(7)
    rows = analyze_step4_case(case)
    assert all(row.h0 == h0 for row in rows)
    rebuilt = random_replacements(case.cnf, R=1, rng=random.Random(0))
    assert rebuilt.rows == case.family.rows
