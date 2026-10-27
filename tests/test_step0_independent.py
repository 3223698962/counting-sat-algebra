"""Step 0: disjoint variable groups multiply the solution counts."""

from __future__ import annotations

import random

from counting_sat_algebra.cnf import CNF, brute_count, concatenate_disjoint
from counting_sat_algebra.independent import product_identity_check, run_independent_suite


def test_hand_built_product_identity() -> None:
    left = CNF(3, [[(0, False), (1, False), (2, True)]])
    right = CNF(3, [[(0, True), (1, False), (2, False)]])
    combined = concatenate_disjoint(left, right)
    assert combined.n == 6
    assert combined.m == 2
    assert brute_count(combined) == brute_count(left) * brute_count(right)


def test_random_product_identity_suite() -> None:
    rows = run_independent_suite(trials=5, seed=0)
    assert len(rows) == 5
    assert all(row["ok"] for row in rows)
    assert all(row["count"] == row["product"] for row in rows)


def test_product_identity_check_one_draw() -> None:
    n, m, count, product, ok = product_identity_check(3, 3, 4, 4, rng=random.Random(2))
    assert n == 7
    assert m == 7
    assert ok
    assert count == product
