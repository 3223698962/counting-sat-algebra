"""Step 1: occurrence ≤ 3 generator/checker and a degree-cap interpolation."""

from __future__ import annotations

import random
from fractions import Fraction

import pytest

from counting_sat_algebra.bounded_occurrence import (
    BENCH_SEED,
    bench_path,
    build_bench_document,
    check_bounded_occurrence,
    cnf_from_record,
    degree_bound_3n,
    equispaced_degree,
    load_bench_document,
    random_bounded_3cnf,
    run_degree_bound_experiment,
)
from counting_sat_algebra.cnf import CNF


def test_equispaced_degree_of_a_quadratic() -> None:
    # t^2 at t = 1, 2, 3, 4.
    degree = equispaced_degree([Fraction(1), Fraction(4), Fraction(9), Fraction(16)])
    assert degree == 2


@pytest.mark.parametrize(
    ("n", "m"),
    [(3, 0), (3, 1), (3, 3), (6, 1), (6, 6), (8, 5), (10, 10)],
)
def test_generator_respects_occurrence_cap(n: int, m: int) -> None:
    cnf = random_bounded_3cnf(n, m, rng=random.Random(1000 * n + m))
    report = check_bounded_occurrence(cnf)
    assert report.ok, report.violations
    assert cnf.n == n
    assert cnf.m == m
    assert report.max_occurrence <= 3
    assert report.total_literal_occurrences == 3 * m
    assert report.total_literal_occurrences <= 3 * n


def test_generator_rejects_infeasible_density() -> None:
    with pytest.raises(ValueError):
        random_bounded_3cnf(5, 6, rng=random.Random(0))


def test_checker_rejects_four_occurrences() -> None:
    cnf = CNF(
        5,
        [
            [(0, False), (1, False), (2, False)],
            [(0, True), (1, True), (3, False)],
            [(0, False), (2, True), (3, True)],
            [(0, True), (1, False), (4, False)],
        ],
    )
    report = check_bounded_occurrence(cnf)
    assert not report.ok
    assert report.max_occurrence == 4


def test_checker_rejects_short_clause() -> None:
    cnf = CNF(3, [[(0, False), (1, True)]])
    report = check_bounded_occurrence(cnf)
    assert not report.ok
    assert report.width_max == 2


def test_degree_bound_on_one_bounded_formula() -> None:
    cnf = random_bounded_3cnf(6, 5, rng=random.Random(7))
    result = run_degree_bound_experiment(cnf, R=1, rng=random.Random(8))
    assert result["bounded_ok"]
    assert result["D_mkR"] == 5 * 3 * 1
    assert result["D_3nR"] == degree_bound_3n(6, 1)
    assert result["D_mkR"] <= result["D_3nR"]
    assert result["loose_cap_applies"]
    assert result["ok_mkR"]
    assert result["ok_3nR"]
    assert result["H0_mkR"] == Fraction(result["brute"])
    assert result["H0_3nR"] == Fraction(result["brute"])
    assert result["observed_degree"] <= result["D_mkR"]


def test_bench_file_matches_fixed_seed_and_degree_caps() -> None:
    path = bench_path()
    document = load_bench_document(path)
    assert document == build_bench_document(BENCH_SEED)
    assert document["seed"] == 42
    for record in document["instances"]:
        cnf = cnf_from_record(record)
        report = check_bounded_occurrence(cnf)
        assert report.ok, report.violations
        assert report.max_occurrence == record["max_occurrence"]
        assert report.max_occurrence <= 3
        result = run_degree_bound_experiment(
            cnf, R=int(record["R"]), rng=random.Random(42)
        )
        assert result["ok_mkR"]
        assert result["ok_3nR"]
        assert result["H0_mkR"] == result["H0_3nR"] == Fraction(result["brute"])
