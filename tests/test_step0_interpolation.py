"""Step 0: H(0) from Lagrange interpolation matches brute-force #F."""

from __future__ import annotations

import random
from fractions import Fraction

from counting_sat_algebra.cnf import CNF, brute_count, random_3cnf
from counting_sat_algebra.replacements import interpolate_H0, random_replacements


def test_single_clause_count_and_interpolation() -> None:
    # x0 ∨ x1 ∨ x2 is false only at 000.
    cnf = CNF(3, [[(0, False), (1, False), (2, False)]])
    assert brute_count(cnf) == 7
    h0, degree_cap, _fam = interpolate_H0(cnf, R=1, rng=random.Random(0), k=3)
    assert degree_cap == 3
    assert h0 == Fraction(7)


def test_replacement_at_zero_restores_the_variable() -> None:
    rng = random.Random(1)
    cnf = random_3cnf(4, 3, rng=rng)
    fam = random_replacements(cnf, R=2, rng=rng)
    assignment = (1, 0, 1, 1)
    for e, clause in enumerate(cnf.clauses):
        for j, (var, _neg) in enumerate(clause):
            value = fam.linear_form(e, j, assignment, Fraction(0))
            assert value == Fraction(assignment[var])


def test_several_random_formulas_match_brute_force() -> None:
    cases = [
        (4, 4, 1, 11),
        (5, 6, 1, 12),
        (4, 3, 2, 13),
        (6, 6, 1, 14),  # n >= 6
    ]
    for n, m, degree, seed in cases:
        cnf = random_3cnf(n, m, rng=random.Random(seed))
        brute = brute_count(cnf)
        h0, degree_cap, _fam = interpolate_H0(
            cnf, R=degree, rng=random.Random(seed + 100), k=3
        )
        assert degree_cap == m * 3 * degree
        assert h0 == Fraction(brute)
        assert h0.denominator == 1
