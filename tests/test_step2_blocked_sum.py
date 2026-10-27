"""Step 2: blocked sum equals full-space H(t) on hand-built replacements.

The measurement enumeration is allowed to be exponential in the number of
blocks. These tests check agreement and the operation counts. They do not
claim a subexponential runtime, the matrix M0, or Theorem 1.1.
"""

from __future__ import annotations

import itertools
import random
from fractions import Fraction

import pytest

from counting_sat_algebra.blocked_sum import (
    blocked_sum,
    compare_blocked_to_full,
    constant_alpha,
    count_at_zero_matches_brute,
    hand_built_examples,
    hand_built_replacement,
    measurement_linear_form,
    popcount_group_instance,
    popcount_growth,
    popcount_operation_counts,
)
from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.design import make_design, popcount_design
from counting_sat_algebra.replacements import degree_bound, eval_H_at, phi_clause_numeric


def _closed_form_full_sum(cnf: CNF, structured, t: Fraction) -> Fraction:
    """Third evaluation of H(t): closed form on the full cube, no blocking."""
    if any(len(clause) == 0 for clause in cnf.clauses):
        return Fraction(0)
    total = Fraction(0)
    for bits in itertools.product((0, 1), repeat=cnf.n):
        measurement = structured.design.measure(bits)
        product = Fraction(1)
        for clause_index, clause in enumerate(cnf.clauses):
            values = [
                measurement_linear_form(
                    structured,
                    clause_index,
                    occ,
                    bits[var],
                    measurement,
                    t,
                )
                for occ, (var, _neg) in enumerate(clause)
            ]
            product *= phi_clause_numeric(values, structured.family.occ_neg[clause_index])
        total += product
    return total


def _independent_enumeration(design) -> tuple[int, int]:
    """Tuple count and local-evaluation count from the block measurement sets."""
    sizes = []
    for group, block in enumerate(design.blocks):
        seen = set()
        for bits in itertools.product((0, 1), repeat=len(block)):
            assignment = {var: bit for var, bit in zip(block, bits)}
            seen.add(design.local_measurement(group, assignment))
        sizes.append((len(block), len(seen)))
    tuples = 1
    for _size, distinct in sizes:
        tuples *= distinct
    local = 0
    for size, distinct in sizes:
        local += (2**size) * (tuples // distinct)
    return tuples, local


def _assert_measurements_add(structured) -> None:
    design = structured.design
    for bits in itertools.product((0, 1), repeat=design.n):
        expected = design.measure(bits)
        parts = []
        for group, block in enumerate(design.blocks):
            assignment = {var: bits[var] for var in block}
            parts.append(design.local_measurement(group, assignment))
        got = tuple(
            sum((part[q] for part in parts), start=Fraction(0))
            for q in range(design.tau)
        )
        assert got == expected


def _assert_rows_match_measurement_form(cnf: CNF, structured, points) -> None:
    for bits in itertools.product((0, 1), repeat=cnf.n):
        measurement = structured.design.measure(bits)
        for clause_index, clause in enumerate(cnf.clauses):
            for occ, (var, _neg) in enumerate(clause):
                for t in points:
                    point = Fraction(t)
                    from_rows = structured.family.linear_form(clause_index, occ, bits, point)
                    from_measurements = measurement_linear_form(
                        structured,
                        clause_index,
                        occ,
                        bits[var],
                        measurement,
                        point,
                    )
                    assert from_rows == from_measurements
                    if point == 0:
                        assert from_rows == Fraction(bits[var])


def test_hand_expanded_two_unit_clauses() -> None:
    # ⟨P, X⟩ = X_v + t (X_0 + X_1) on two positive unit clauses.
    # Expanding the four assignments gives H(t) = 1 + 6t + 6t^2.
    cnf = CNF(2, [[(0, False)], [(1, False)]])
    design = popcount_design(2, [(0,), (1,)], D=degree_bound(cnf, R=1, k=1), r=1)
    structured = hand_built_replacement(
        cnf,
        design,
        R=1,
        alpha=constant_alpha(cnf, R=1, tau=1, value=Fraction(1)),
    )
    assert brute_count(cnf) == 1
    expected = {
        Fraction(0): Fraction(1),
        Fraction(1): Fraction(13),
        Fraction(2): Fraction(37),
        Fraction(1, 2): Fraction(11, 2),
    }
    for t, value in expected.items():
        assert value == Fraction(1) + 6 * t + 6 * t * t
        full = eval_H_at(cnf, structured.family, t)
        blocked = blocked_sum(cnf, structured, t)
        closed = _closed_form_full_sum(cnf, structured, t)
        assert full == value
        assert blocked.value == value
        assert closed == value
    tuples, local, full_space = popcount_operation_counts(2, 1)
    blocked = blocked_sum(cnf, structured, Fraction(1))
    assert (blocked.measurement_tuples, blocked.local_evaluations, blocked.full_space) == (
        tuples,
        local,
        full_space,
    )


def test_hand_expanded_negative_unit() -> None:
    # One negative unit, ⟨P, X⟩ = X (1 + t). Then H(t) = 1 - t, and H(0) = #F = 1.
    cnf = CNF(1, [[(0, True)]])
    design = popcount_design(1, [(0,)], D=1, r=1)
    structured = hand_built_replacement(
        cnf,
        design,
        R=1,
        alpha=((((Fraction(1),),),),),
    )
    assert brute_count(cnf) == 1
    for t in (Fraction(0), Fraction(1), Fraction(2), Fraction(3, 2)):
        expected = Fraction(1) - t
        assert eval_H_at(cnf, structured.family, t) == expected
        assert blocked_sum(cnf, structured, t).value == expected


def test_examples_agree_with_full_sum_and_have_several_blocks() -> None:
    examples = hand_built_examples()
    assert len(examples) >= 4
    assert any(len(example.structured.design.blocks) > 1 for example in examples)
    expected_counts = {
        "one-block-n3": 7,
        "two-block-3cnf-n6": 42,
        "three-block-n6": 27,
        "two-block-tau2-R2-n6": 42,
    }
    expected_params = {
        "one-block-n3": (3, 1, 3, 1),
        "two-block-3cnf-n6": (3, 1, 9, 2),
        "three-block-n6": (2, 1, 6, 3),
        "two-block-tau2-R2-n6": (3, 2, 18, 2),
    }
    for example in examples:
        cnf = example.cnf
        structured = example.structured
        design = structured.design
        assert brute_count(cnf) == expected_counts[example.name]
        assert count_at_zero_matches_brute(cnf, structured)
        r, tau, degree_cap, blocks = expected_params[example.name]
        assert (design.r, design.tau, design.D, len(design.blocks)) == (
            r,
            tau,
            degree_cap,
            blocks,
        )
        _assert_measurements_add(structured)
        _assert_rows_match_measurement_form(cnf, structured, example.sample_points)
        tuples, local = _independent_enumeration(design)
        for t, full, blocked in compare_blocked_to_full(cnf, structured, example.sample_points):
            assert blocked.value == full
            assert _closed_form_full_sum(cnf, structured, t) == full
            assert blocked.measurement_tuples == tuples
            assert blocked.local_evaluations == local
            assert blocked.full_space == 2**cnf.n
            assert blocked.degree_cap == design.D
        if design.tau == 1 and all(weight == Fraction(1) for weight in design.weights[0]):
            formula = popcount_operation_counts(len(design.blocks), design.r)
            assert (tuples, local, 2**cnf.n) == formula


def test_measurement_tuples_grow_with_the_number_of_groups() -> None:
    # Equal blocks of size 2, one popcount measurement. Tuple count is 3^b.
    # b = 1 was also expanded by hand: H(0) = 3 and H(1) = 11.
    tuples_seen = []
    local_seen = []
    for groups in range(1, 5):
        cnf, structured = popcount_group_instance(groups, block_size=2, R=1)
        assert brute_count(cnf) == 3**groups
        for t, expected_count in ((Fraction(0), 3**groups), (Fraction(1), None)):
            blocked = blocked_sum(cnf, structured, t)
            full = eval_H_at(cnf, structured.family, t)
            formula = popcount_operation_counts(groups, 2)
            assert blocked.value == full
            assert (
                blocked.measurement_tuples,
                blocked.local_evaluations,
                blocked.full_space,
            ) == formula
            assert blocked.seconds >= 0
            if t == 0:
                assert blocked.value == Fraction(expected_count)
                assert blocked.value.denominator == 1
        blocked_at_one = blocked_sum(cnf, structured, Fraction(1))
        tuples_seen.append(blocked_at_one.measurement_tuples)
        local_seen.append(blocked_at_one.local_evaluations)
    assert tuples_seen == [3, 9, 27, 81]
    assert local_seen == [4, 24, 108, 432]
    assert all(tuples_seen[i] < tuples_seen[i + 1] for i in range(3))
    assert all(local_seen[i] < local_seen[i + 1] for i in range(3))

    cnf, structured = popcount_group_instance(1, block_size=2, R=1)
    assert blocked_sum(cnf, structured, Fraction(1)).value == Fraction(11)

    growth = popcount_growth((1, 2, 3, 4, 5), block_size=2, R=1, t=1)
    assert [row.groups for row in growth] == [1, 2, 3, 4, 5]
    assert all(row.ok for row in growth)
    assert [row.blocked.measurement_tuples for row in growth] == [3, 9, 27, 81, 243]
    assert [row.blocked.full_space for row in growth] == [4, 16, 64, 256, 1024]


def test_block_size_three_local_count_can_sit_below_the_cube() -> None:
    # The identity does not care which side is smaller. For three blocks of size 3
    # the popcount formula gives 384 local evaluations and a cube of size 512.
    # The tuple count is 4^3 = 64, still exponential in the number of blocks.
    groups, block_size = 3, 3
    cnf, structured = popcount_group_instance(groups, block_size=block_size, R=1)
    blocked = blocked_sum(cnf, structured, Fraction(1))
    full = eval_H_at(cnf, structured.family, Fraction(1))
    tuples, local, space = popcount_operation_counts(groups, block_size)
    assert (tuples, local, space) == (64, 384, 512)
    assert blocked.value == full
    assert (blocked.measurement_tuples, blocked.local_evaluations, blocked.full_space) == (
        tuples,
        local,
        space,
    )
    assert local < space
    assert blocked_sum(cnf, structured, Fraction(0)).value == Fraction(brute_count(cnf))


def test_block_size_one_measurement_space_is_the_full_cube() -> None:
    cnf, structured = popcount_group_instance(3, block_size=1, R=1)
    blocked = blocked_sum(cnf, structured, Fraction(1))
    assert blocked.measurement_tuples == 8
    assert blocked.full_space == 8
    assert blocked.value == eval_H_at(cnf, structured.family, Fraction(1))
    assert brute_count(cnf) == 1
    assert blocked_sum(cnf, structured, Fraction(0)).value == Fraction(1)


def test_zero_measurements_are_one_tuple() -> None:
    # τ = 0: no shared measurement. The blocked sum is one product of block sums.
    cnf = CNF(4, [[(0, False), (1, True)], [(2, False), (3, False)]])
    design = make_design(
        4,
        [(0, 1), (2, 3)],
        weights=[],
        D=degree_bound(cnf, R=1, k=2),
        r=2,
    )
    structured = hand_built_replacement(
        cnf,
        design,
        R=1,
        alpha=constant_alpha(cnf, R=1, tau=0, value=Fraction(1)),
    )
    assert brute_count(cnf) == 9
    for t in (Fraction(0), Fraction(1), Fraction(4)):
        blocked = blocked_sum(cnf, structured, t)
        assert blocked.measurement_tuples == 1
        assert blocked.local_evaluations == 8
        assert blocked.value == Fraction(9)
        assert blocked.value == eval_H_at(cnf, structured.family, t)


def test_rational_weights_agree() -> None:
    cnf = CNF(2, [[(0, False)], [(1, True)]])
    design = make_design(
        2,
        [(0,), (1,)],
        [[Fraction(1, 2), Fraction(1, 2)]],
        D=degree_bound(cnf, R=1, k=1),
        r=1,
    )
    structured = hand_built_replacement(
        cnf,
        design,
        R=1,
        alpha=constant_alpha(cnf, R=1, tau=1, value=Fraction(1)),
    )
    _assert_measurements_add(structured)
    for t, full, blocked in compare_blocked_to_full(
        cnf, structured, (0, 1, Fraction(1, 2))
    ):
        assert blocked.value == full
        assert _closed_form_full_sum(cnf, structured, t) == full
    assert blocked_sum(cnf, structured, 0).value == Fraction(brute_count(cnf))


def test_seeded_random_alpha_agrees() -> None:
    cnf = CNF(
        6,
        [
            [(0, False), (1, True), (2, False)],
            [(3, True), (4, False), (5, True)],
        ],
    )
    design = popcount_design(6, [(0, 1, 2), (3, 4, 5)], D=degree_bound(cnf, R=2, k=3), r=3)
    structured = hand_built_replacement(cnf, design, R=2, rng=random.Random(2026))
    _assert_rows_match_measurement_form(cnf, structured, (Fraction(0), Fraction(1)))
    for _t, full, blocked in compare_blocked_to_full(cnf, structured, (0, 1, 2)):
        assert blocked.value == full
    assert count_at_zero_matches_brute(cnf, structured)


def test_r0_still_enumerates_measurements() -> None:
    # Degree 0 replacements ignore μ, so H is the constant #F, but the popcount
    # tuples are still enumerated.
    cnf, structured = popcount_group_instance(2, block_size=2, R=0)
    assert structured.R == 0
    assert structured.design.D == 0
    blocked = blocked_sum(cnf, structured, Fraction(5))
    assert blocked.measurement_tuples == 9
    assert blocked.value == Fraction(brute_count(cnf))
    assert blocked.value == eval_H_at(cnf, structured.family, Fraction(5))


def test_no_clauses_sums_to_two_to_the_n() -> None:
    cnf = CNF(4, [])
    design = popcount_design(4, [(0, 1), (2, 3)], D=0, r=2)
    structured = hand_built_replacement(cnf, design, R=1, alpha=())
    blocked = blocked_sum(cnf, structured, Fraction(1))
    assert blocked.value == 16
    assert blocked.measurement_tuples == 9
    assert eval_H_at(cnf, structured.family, Fraction(1)) == 16


def test_empty_clause_is_identically_zero() -> None:
    cnf = CNF(2, [[]])
    design = popcount_design(2, [(0, 1)], D=0, r=2)
    structured = hand_built_replacement(cnf, design, R=1, alpha=((),))
    blocked = blocked_sum(cnf, structured, Fraction(1))
    assert blocked.value == 0
    assert blocked.measurement_tuples == 0
    assert eval_H_at(cnf, structured.family, Fraction(1)) == 0


def test_design_rejects_illegal_partitions_and_crossing_clauses() -> None:
    with pytest.raises(ValueError, match="repeated"):
        make_design(3, [(0, 1), (1, 2)], [[1, 1, 1]])
    with pytest.raises(ValueError, match="missing"):
        make_design(3, [(0, 1)], [[1, 1, 1]])
    with pytest.raises(ValueError, match="larger than r"):
        make_design(3, [(0, 1, 2)], [[1, 1, 1]], r=2)
    with pytest.raises(ValueError, match="D"):
        make_design(1, [(0,)], [[1]], D=-1)
    with pytest.raises(ValueError, match="Fraction"):
        make_design(2, [(0,), (1,)], [[1.5, 0.5]])

    crossing = CNF(4, [[(0, False), (2, False)]])
    design = popcount_design(4, [(0, 1), (2, 3)], D=0, r=2)
    with pytest.raises(ValueError, match="one block"):
        hand_built_replacement(
            crossing,
            design,
            R=1,
            alpha=constant_alpha(crossing, R=1, tau=1),
        )
    local = CNF(4, [[(0, False), (1, True)]])
    with pytest.raises(ValueError, match="alpha"):
        hand_built_replacement(local, design, R=1, alpha=[])


def test_enumeration_guards_refuse_large_measurement_spaces() -> None:
    cnf, structured = popcount_group_instance(2, block_size=2, R=1)
    with pytest.raises(ValueError, match="exponential"):
        blocked_sum(cnf, structured, 1, max_measurement_tuples=8)
    with pytest.raises(ValueError, match="enumeration cap"):
        blocked_sum(cnf, structured, 1, max_block_size=1)
