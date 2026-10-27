"""Step 5: counting-preserving copying down to occurrence ≤ 3.

#F' = #F on instances inside the brute-force window. The reduced formula is a
mixed-width CNF, so the exact-width-3 checker still rejects a formula that
used the cycle. These tests do not claim Theorem 1.1, a subexponential
runtime, or 2^{o(n)}.
"""

from __future__ import annotations

import random
from fractions import Fraction

import pytest

from counting_sat_algebra.bounded_occurrence import (
    check_bounded_occurrence,
    random_bounded_3cnf,
    run_degree_bound_experiment,
)
from counting_sat_algebra.cnf import CNF, brute_count, is_satisfied
from counting_sat_algebra.sparsify import (
    BRUTE_N_MAX,
    EXTENSION_N_MAX,
    MIN_OCCURRENCE_CAP,
    PIPELINE_ASSIGNMENT_BUDGET,
    RECOVERY_FORMULA,
    analyze_step5_case,
    copies_for_occurrence,
    copying_lower_bound,
    cycle_clauses,
    expand_blocks,
    extension_problems,
    heavy_3cnf,
    project_assignment,
    rewritten_source_clauses,
    run_step1_pipeline,
    sparsify,
    step5_cases,
    two_block_formula,
)


def _units(count: int, n_vars: int = 1) -> CNF:
    return CNF(n_vars, [[(0, False)] for _ in range(count)])


def test_copies_for_occurrence_matches_the_cycle_slots() -> None:
    assert MIN_OCCURRENCE_CAP == 3
    assert copies_for_occurrence(0) == 1
    assert copies_for_occurrence(3) == 1
    assert copies_for_occurrence(4) == 4
    assert copies_for_occurrence(6, cap=4) == 3
    with pytest.raises(ValueError, match="at least 3"):
        copies_for_occurrence(4, cap=2)
    with pytest.raises(ValueError):
        sparsify(_units(1), max_occurrence=True)  # type: ignore[arg-type]


def test_four_units_have_an_explicit_cycle_and_preserve_the_count() -> None:
    spec = sparsify(_units(4))
    assert spec.recovery == RECOVERY_FORMULA
    assert spec.n_after == 4
    assert spec.m_after == 8
    assert spec.cycle_clause_count == 4
    assert spec.max_occurrence_after == 3
    assert spec.literal_occurrences_after == 12
    assert rewritten_source_clauses(spec) == [[(i, False)] for i in range(4)]
    assert cycle_clauses(spec) == [
        [(0, True), (1, False)],
        [(1, True), (2, False)],
        [(2, True), (3, False)],
        [(3, True), (0, False)],
    ]
    assert is_satisfied((1, 1, 1, 1), spec.reduced)
    assert not is_satisfied((1, 1, 1, 0), spec.reduced)
    assert not is_satisfied((0, 0, 0, 0), spec.reduced)
    assert project_assignment(spec, (1, 1, 1, 0)) is None
    assert extension_problems(spec) == []
    report = check_bounded_occurrence(spec.reduced)
    assert not report.ok


def test_heavy_3cnf_copies_only_the_heavy_variable() -> None:
    source = heavy_3cnf()
    saved = [list(clause) for clause in source.clauses]
    spec = sparsify(source)
    assert source.clauses == saved
    assert spec.n_before == 4 and spec.n_after == 7
    assert spec.m_before == 4 and spec.m_after == 8
    assert spec.max_occurrence_before == 4
    assert spec.max_occurrence_after == 3
    assert [len(info.copies) for info in spec.variables] == [4, 1, 1, 1]
    assert rewritten_source_clauses(spec) == [
        [(0, False), (4, False), (5, False)],
        [(1, False), (4, False), (6, False)],
        [(2, False), (5, False), (6, False)],
        [(3, True), (4, True), (5, True)],
    ]
    assert extension_problems(spec) == []


def test_already_bounded_formula_is_unchanged() -> None:
    source = CNF(3, [[(0, False), (1, False), (2, True)]])
    spec = sparsify(source)
    assert spec.cycle_clause_count == 0
    assert spec.reduced.n == source.n
    assert spec.reduced.clauses == source.clauses
    wider = random_bounded_3cnf(6, 4, rng=random.Random(7))
    copied = sparsify(wider)
    assert copied.reduced.clauses == wider.clauses
    assert copied.max_occurrence_after <= 3
    assert check_bounded_occurrence(copied.reduced).ok
    result = run_degree_bound_experiment(copied.reduced, R=1, rng=random.Random(8))
    assert result["ok_mkR"] and result["ok_3nR"]
    assert result["H0_mkR"] == Fraction(result["brute"])


def test_occurrence_three_meets_the_lower_bound_with_equality() -> None:
    source = CNF(3, [[(0, False), (1, False), (2, False)] for _ in range(3)])
    spec = sparsify(source)
    assert spec.literal_occurrences_before == 9
    assert spec.n_after == copying_lower_bound(9, 3) == 3
    assert spec.cycle_clause_count == 0


def test_cap_four_uses_fewer_copies_and_still_preserves_the_count() -> None:
    spec = sparsify(_units(6), max_occurrence=4)
    assert [len(info.copies) for info in spec.variables] == [3]
    assert spec.n_after == 3
    assert spec.m_after == 9
    assert spec.max_occurrence_after == 4
    assert spec.cycle_clause_count == 3
    assert extension_problems(spec) == []
    pipeline = run_step1_pipeline(spec, R=1, rng=random.Random(4))
    assert pipeline.status == "agreed"
    assert pipeline.h0 == Fraction(1)
    assert pipeline.observed_degree is not None
    assert pipeline.observed_degree <= pipeline.degree_cap_mass


def test_width_four_clause_preserves_the_count() -> None:
    source = CNF(
        5,
        [
            [(0, False), (1, False), (2, False), (3, False)],
            [(0, False), (1, False), (2, True), (4, False)],
            [(0, True), (1, False), (3, False), (4, True)],
            [(0, False), (2, False), (3, True), (4, False)],
        ],
    )
    spec = sparsify(source)
    assert spec.max_occurrence_before == 4
    assert spec.max_occurrence_after <= 3
    assert spec.n_after == 8
    assert max(len(clause) for clause in spec.reduced.clauses) == 4
    # Cycle clauses have width 2, so neither exact-width checker accepts the output.
    assert not check_bounded_occurrence(spec.reduced).ok
    assert not check_bounded_occurrence(spec.reduced, width=4).ok
    assert extension_problems(spec) == []
    pipeline = run_step1_pipeline(spec, R=0, rng=random.Random(0))
    assert pipeline.status == "agreed"
    assert pipeline.degree_cap == 0
    assert pipeline.h0 == Fraction(brute_count(source))


def test_empty_clause_and_empty_formula() -> None:
    empty_clause = sparsify(CNF(2, [[]]))
    assert empty_clause.reduced.clauses == [[]]
    assert extension_problems(empty_clause) == []
    assert brute_count(empty_clause.reduced) == 0
    no_clauses = sparsify(CNF(3, []))
    assert no_clauses.reduced.n == 3
    assert no_clauses.reduced.m == 0
    assert brute_count(no_clauses.reduced) == 8
    pipeline = run_step1_pipeline(no_clauses, R=1, rng=random.Random(1))
    assert pipeline.status == "agreed"
    assert pipeline.h0 == Fraction(8)
    assert pipeline.observed_degree == 0


def test_split_contradiction_has_no_models() -> None:
    source = CNF(
        1,
        [[(0, False)], [(0, False)], [(0, False)], [(0, False)], [(0, True)]],
    )
    spec = sparsify(source)
    assert spec.n_after == 5
    assert spec.max_occurrence_after == 3
    assert extension_problems(spec) == []
    assert brute_count(spec.source) == brute_count(spec.reduced) == 0


def test_repeated_variable_and_bad_partitions_are_refused() -> None:
    with pytest.raises(ValueError, match="repeats"):
        sparsify(CNF(2, [[(0, False), (0, True), (1, False)]]))
    spec = sparsify(CNF(3, [[(0, False), (1, False), (2, False)]]))
    with pytest.raises(ValueError, match="missing"):
        expand_blocks(spec, ((0, 1),))
    with pytest.raises(ValueError, match="repeated"):
        expand_blocks(spec, ((0, 0), (1, 2)))


def test_two_block_copies_stay_inside_the_source_blocks() -> None:
    spec = sparsify(two_block_formula())
    assert expand_blocks(spec, ((0, 1), (2,))) == (
        (0, 1, 2, 3, 4),
        (5, 6, 7, 8),
    )
    assert spec.n_after == 9
    assert spec.literal_occurrences_after == 26
    for clause in cycle_clauses(spec):
        variables = {var for var, _neg in clause}
        assert variables <= {0, 1, 2, 3} or variables <= {5, 6, 7, 8}


def test_brute_wall_is_the_same_cutoff_as_brute_count() -> None:
    spec = sparsify(_units(BRUTE_N_MAX + 1))
    assert spec.n_after == 21
    assert spec.max_occurrence_after == 3
    with pytest.raises(ValueError, match="20"):
        brute_count(spec.reduced)
    with pytest.raises(ValueError, match="above"):
        extension_problems(spec)


def test_step5_cases_match_the_recorded_blowup() -> None:
    """Locks the blowup table. Expected refusals are passing outcomes."""
    expected = {
        "units-n1-m4": {
            "n": (1, 4),
            "m": (4, 8),
            "occ": (4, 3),
            "count": 1,
            "degree": 12,
            "caps": (16, 12, 12),
        },
        "free-bit-n2": {
            "n": (2, 5),
            "m": (4, 8),
            "occ": (4, 3),
            "count": 2,
            "degree": 12,
            "caps": (16, 12, 15),
        },
        "heavy-3cnf-n4": {
            "n": (4, 7),
            "m": (4, 8),
            "occ": (4, 3),
            "count": 10,
            "degree": 20,
            "caps": (24, 20, 21),
        },
        "bounded-clause-n3": {
            "n": (3, 3),
            "m": (1, 1),
            "occ": (1, 1),
            "count": 7,
            "degree": 3,
            "caps": (3, 3, 9),
        },
        "empty-clause-n2": {
            "n": (2, 2),
            "m": (1, 1),
            "occ": (0, 0),
            "count": 0,
            "degree": 0,
            "caps": (0, 0, 6),
        },
        "two-block-copies": {
            "n": (3, 9),
            "m": (8, 16),
            "occ": (4, 3),
            "count": 2,
            "tuples": 30,
            "local": 256,
            "full": 512,
        },
        "cross-block-copies": {
            "n": (3, 6),
            "m": (4, 8),
            "occ": (4, 3),
            "count": 4,
            "degree": 14,
            "caps": (24, 14, 18),
        },
        "units-n1-m16": {
            "n": (1, 16),
            "m": (16, 32),
            "occ": (16, 3),
            "count": 1,
            "assignments": 8_650_752,
        },
        "units-n1-m21": {"n": (1, 21), "m": (21, 42), "occ": (21, 3), "source": 1},
        "random-3cnf-n5-m12": {
            "n": (5, 36),
            "m": (12, 48),
            "occ": (8, 3),
            "source": 8,
            "lower": 26,
        },
    }
    cases = step5_cases()
    assert [case.name for case in cases] == list(expected)
    for case in cases:
        report = analyze_step5_case(case)
        assert report.ok, report.problems
        spec = report.sparsification
        row = expected[case.name]
        assert (spec.n_before, spec.n_after) == row["n"]
        assert (spec.m_before, spec.m_after) == row["m"]
        assert (spec.max_occurrence_before, spec.max_occurrence_after) == row["occ"]
        assert spec.max_occurrence_after <= 3
        assert spec.literal_occurrences_after <= 3 * spec.n_after
        assert spec.n_after >= spec.n_before
        if spec.cap == 3:
            assert spec.n_after >= copying_lower_bound(
                spec.literal_occurrences_before, spec.n_before
            )
        if "count" in row:
            assert report.count.status == "agreed"
            assert report.count.brute_source == report.count.brute_reduced == row["count"]
        if "source" in row:
            assert report.count.status == "refused"
            assert report.count.brute_source == row["source"]
            assert report.count.brute_reduced is None
            assert spec.n_after > BRUTE_N_MAX
            assert report.pipeline is not None
            assert report.pipeline.status == "refused"
        if "degree" in row:
            assert report.pipeline is not None and report.pipeline.status == "agreed"
            pipeline = report.pipeline
            assert pipeline.h0 == Fraction(row["count"])
            assert pipeline.observed_degree == row["degree"]
            assert (
                pipeline.degree_cap,
                pipeline.degree_cap_mass,
                pipeline.degree_cap_loose,
            ) == row["caps"]
            assert pipeline.observed_degree <= pipeline.degree_cap_mass
        if "tuples" in row:
            assert report.blocked is not None and report.blocked.status == "agreed"
            assert report.blocked.blocked_at_zero == Fraction(row["count"])
            assert report.blocked.tuples == row["tuples"]
            assert report.blocked.local_evaluations == row["local"]
            assert report.blocked.full_space == row["full"]
            assert report.blocked.blocked_at_one == Fraction(
                -47986541792015699006680248
            )
        if case.name == "cross-block-copies":
            assert report.blocked is not None and report.blocked.status == "refused"
            assert "one block" in report.blocked.reason
        if case.name == "units-n1-m16":
            assert report.pipeline is not None
            assert report.pipeline.status == "refused"
            assert report.pipeline.assignments == row["assignments"]
            assert report.pipeline.assignments > PIPELINE_ASSIGNMENT_BUDGET
        if case.name == "random-3cnf-n5-m12":
            assert spec.n_after >= row["lower"]
            assert copying_lower_bound(3 * 12, 5) == row["lower"]
            assert spec.n_after > EXTENSION_N_MAX
