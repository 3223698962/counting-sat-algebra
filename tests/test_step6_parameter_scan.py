"""Step 6: a machine-limited parameter scan against brute force.

The quick profile is what pytest runs. It locks matches on the small grid and
refusals where a construction is unfinished or over a rail. The full profile
is the bench document and ``verify --step6``. Wall times are not a complexity
result. These tests do not claim Theorem 1.1, a subexponential runtime, or
2^{o(n)}.
"""

from __future__ import annotations

import pytest

from counting_sat_algebra.blocked_sum import popcount_operation_counts
from counting_sat_algebra.parameter_scan import (
    ASSIGNMENT_BUDGET,
    DISTANCE_TO_FULL_ALGORITHM,
    TIMING_DISCLAIMER,
    degree_plan,
    hermite_paired_plan,
    lagrange_assignments,
    load_scan_document,
    run_parameter_scan,
    scan_specs,
)
from counting_sat_algebra.replacement_matrix import FULL_M0_LOG2_DIMENSION_AT_R2
from counting_sat_algebra.sparsify import BRUTE_N_MAX


def _by_name(report):
    return {cell.name: cell for cell in report.cells}


def test_quick_scan_matches_the_locked_outcomes() -> None:
    report = run_parameter_scan("quick")
    assert report.ok
    assert report.disclaimer == TIMING_DISCLAIMER
    assert report.distance == DISTANCE_TO_FULL_ALGORITHM
    cells = _by_name(report)
    assert set(cells) == {spec.name for spec in scan_specs("quick")}

    one = cells["blocked-r2-b1"]
    assert one.status == "agreed" and one.measured
    assert one.brute == 3 and one.compared == 3
    assert (one.tuples, one.local_evaluations, one.cube) == (3, 4, 4)
    assert one.seconds is not None and one.seconds_brute is not None and one.seconds_full is not None

    two = cells["blocked-r2-b2"]
    assert two.brute == 9 and two.tuples == 9 and two.local_evaluations == 24
    assert two.cube == 16 and two.local_evaluations > two.cube

    three = cells["blocked-r2-b3"]
    assert three.brute == 27 and three.tuples == 27 and three.local_evaluations == 108
    assert three.cube == 64

    wider = cells["blocked-r2-b2-R2"]
    assert wider.brute == two.brute == 9
    assert wider.tuples == two.tuples == 9
    assert wider.degree_cap == 8 and two.degree_cap == 4

    units = cells["blocked-r1-b4"]
    assert units.brute == 1 and units.tuples == units.cube == 16

    capped = cells["blocked-r1-b17-tuple-cap"]
    assert capped.status == "refused" and not capped.measured
    assert capped.n == 17 and capped.n <= BRUTE_N_MAX
    assert capped.brute == 1 and capped.seconds_brute is not None
    assert capped.seconds is None and capped.seconds_full is None
    assert capped.tuples == 2**17
    assert "above the cap" in capped.reason
    assert "not timed" in capped.reason

    wall = cells["blocked-r2-b11-brute-wall"]
    assert wall.status == "refused" and not wall.measured
    assert wall.n == 22 and wall.brute is None
    assert wall.seconds is None and wall.seconds_brute is None and wall.seconds_full is None
    assert wall.tuples == 3**11 and wall.cube == 2**22
    assert "refuses" in wall.reason and "not started" in wall.reason
    assert "No timing" in wall.reason

    for name, brute in (("lagrange-r2-b1", 3), ("lagrange-r2-b2", 9)):
        cell = cells[name]
        assert cell.status == "agreed" and cell.compared == brute == cell.brute
        assert cell.assignments == (cell.degree_cap + 2) * cell.cube

    lagrange_stop = cells["lagrange-r2-b7-budget"]
    assert lagrange_stop.status == "refused" and lagrange_stop.brute is None
    assert lagrange_stop.seconds is None and lagrange_stop.seconds_brute is None
    assert lagrange_stop.assignments > ASSIGNMENT_BUDGET
    assert "not computed" in lagrange_stop.reason

    paired = cells["hermite-paired-r2-b1"]
    assert paired.status == "agreed" and paired.brute == 3
    assert paired.cube_passes == 2 and paired.jet_width == 2
    assert paired.degree_cap == 2
    assert paired.assignments == (paired.cube_passes + 1) * paired.cube
    assert "D+1" in paired.reason
    assert "single-point" in paired.reason

    hermite_stop = cells["hermite-paired-r2-b7-budget"]
    assert hermite_stop.status == "refused" and hermite_stop.seconds is None
    assert hermite_stop.assignments > ASSIGNMENT_BUDGET
    assert "not started" in hermite_stop.reason

    jet = cells["hermite-single-point-not-run"]
    assert jet.status == "refused" and not jet.measured
    assert jet.seconds is None and jet.seconds_brute is None
    assert jet.cube_passes == 1 and jet.jet_width == 3 and jet.assignments == 4
    assert "was not started" in jet.reason
    assert "exponential cost" in jet.reason

    clause = cells["degree-clause-n3"]
    assert clause.status == "agreed" and clause.brute == 7 and clause.compared == 7
    assert clause.observed_degree is not None and clause.observed_degree <= clause.degree_cap

    bench = cells["degree-bounded-n6-m4"]
    assert bench.status == "agreed" and bench.brute == 44 and bench.compared == 44
    assert bench.observed_degree == 12 and bench.degree_cap == 12

    degree_stop = cells["degree-bounded-n12-m12-budget"]
    assert degree_stop.status == "refused" and degree_stop.brute is None
    assert degree_stop.seconds is None and degree_stop.seconds_brute is None
    _tight, _loose, _window, assignments = degree_plan(12, 12, 1)
    assert degree_stop.assignments == assignments > ASSIGNMENT_BUDGET
    assert "was not sampled" in degree_stop.reason

    copied = cells["copy-popcount-b1"]
    assert copied.status == "agreed" and copied.brute == 3
    assert copied.n_after == 2 and copied.occurrence_after == 1

    blowup = cells["copy-blowup-n5-m12"]
    assert blowup.status == "refused" and blowup.brute == 8
    assert blowup.n_after == 36 and blowup.lower_bound_n == 26
    assert blowup.occurrence_after == 3
    assert blowup.assignments == 2**36
    assert blowup.seconds is None and blowup.seconds_brute is not None
    assert "not timed" in blowup.reason

    units_wall = cells["copy-units-n21"]
    assert units_wall.status == "refused" and units_wall.brute == 1
    assert units_wall.n_after == 21 and units_wall.seconds is None

    cross = cells["copy-cross-block"]
    assert cross.status == "refused" and cross.brute == 4 and cross.n_after == 6
    assert cross.seconds is None and cross.seconds_brute is not None
    assert "block" in cross.reason and "Step 3 is not applied" in cross.reason
    assert "not timed" in cross.reason

    swap = cells["matrix-swap-4"]
    assert swap.status == "agreed" and swap.brute == 1 and swap.compared == 1
    assert swap.observed_degree == 0
    assert "side length 4" in swap.reason and "not full M0" in swap.reason

    full = cells["matrix-full-m0"]
    assert full.status == "refused" and not full.measured
    assert full.side_log2 == FULL_M0_LOG2_DIMENSION_AT_R2 == 384
    assert full.seconds is None and full.seconds_brute is None and full.seconds_full is None
    assert full.brute is None and full.compared is None
    assert "before allocating" in full.reason
    assert "No timing" in full.reason


def test_refused_pipelines_are_not_timed() -> None:
    report = run_parameter_scan("quick")
    for cell in report.cells:
        if cell.status != "agreed":
            assert cell.seconds is None
            assert cell.compared is None or cell.name == "copy-popcount-b1"
        if cell.status == "refused":
            assert cell.seconds is None
            assert "2^{o(n)}" not in cell.reason or "not" in cell.reason


def test_full_profile_extends_quick_without_running_it() -> None:
    quick = scan_specs("quick")
    full = scan_specs("full")
    quick_names = [spec.name for spec in quick]
    full_names = [spec.name for spec in full]
    assert full_names[: len(quick_names)] == quick_names
    assert len(full) > len(quick)
    by_name = {spec.name: spec for spec in full}
    assert by_name["blocked-r2-b7"].expect == "agreed"
    assert by_name["blocked-r2-b7"].n == 14
    assert by_name["lagrange-r2-b7-budget"].expect == "refused"
    assert by_name["lagrange-r2-b7-budget"].n == 14
    assert by_name["hermite-paired-r2-b7-budget"].expect == "refused"
    assert by_name["copy-random-n4-m8"].expect == "refused"
    assert by_name["blocked-r3-b3"].expect == "agreed"
    with pytest.raises(ValueError, match="profile"):
        scan_specs("everything")


def test_closed_forms_match_the_spec_predicates() -> None:
    for spec in scan_specs("full"):
        if spec.step == "blocked":
            tuples, local, cube = popcount_operation_counts(spec.groups, spec.block_size)
            assert spec.n == spec.groups * spec.block_size
            assert cube == 2**spec.n
            assert tuples == (spec.block_size + 1) ** spec.groups
            assert local == spec.groups * (2**spec.block_size) * (spec.block_size + 1) ** (spec.groups - 1)
        if spec.step == "lagrange" and spec.expect == "refused":
            _n, _cap, assignments = lagrange_assignments(spec.groups, spec.block_size, spec.R)
            assert assignments > ASSIGNMENT_BUDGET or spec.n > BRUTE_N_MAX
        if spec.step == "hermite" and spec.kind == "hermite" and spec.expect == "refused":
            _n, _cap, _passes, _width, assignments = hermite_paired_plan(
                spec.groups, spec.block_size, spec.R
            )
            assert assignments > ASSIGNMENT_BUDGET or spec.n > BRUTE_N_MAX


def test_distance_lists_unfinished_work_and_does_not_overclaim() -> None:
    text = "\n".join(DISTANCE_TO_FULL_ALGORITHM)
    for phrase in (
        "Full M0",
        "2^384",
        "scheduler",
        "n' = n + o(n)",
        "3^{n/2}",
        "Hermite",
        "quadratic",
        "cross",
        "Theorem 1.1 is not implemented",
        "does not refute ETH",
    ):
        assert phrase in text
    lowered = text.lower()
    assert "eth is false" not in lowered
    assert "refutes eth" not in lowered
    assert "runs in 2^{o(n)}" not in lowered
    assert TIMING_DISCLAIMER.startswith("Wall times are machine-local")
    assert "not a complexity result" in TIMING_DISCLAIMER


def test_bench_document_locks_the_full_scan() -> None:
    document = load_scan_document()
    specs = scan_specs("full")
    assert document["profile"] == "full"
    assert document["disclaimer"] == TIMING_DISCLAIMER
    assert document["distance"] == list(DISTANCE_TO_FULL_ALGORITHM)
    assert "not a complexity result" in document["description"]
    assert [cell["name"] for cell in document["cells"]] == [spec.name for spec in specs]
    by_spec = {spec.name: spec for spec in specs}
    for cell in document["cells"]:
        spec = by_spec[cell["name"]]
        assert cell["status"] == spec.expect == cell["expect"]
        assert cell["problems"] == []
        if not spec.times_pipeline:
            assert cell["seconds"] is None
        if not spec.times_brute:
            assert cell["seconds_brute"] is None
        if cell["status"] == "refused":
            assert cell["seconds"] is None
            assert cell["measured"] is False
        if spec.step == "blocked" and spec.block_size == 2 and cell["status"] == "agreed":
            assert cell["brute"] == 3 ** spec.groups
            assert cell["compared"] == str(cell["brute"])
            assert cell["tuples"] == 3 ** spec.groups
            assert cell["cube"] == 2 ** (2 * spec.groups)
            assert cell["measured"] is True
            assert cell["seconds"] is not None
        if spec.name == "blocked-r3-b3":
            assert cell["brute"] == 7**3 == 343
            assert cell["tuples"] == 64 and cell["local_evaluations"] == 384
            assert cell["cube"] == 512
        if spec.name == "blocked-r2-b7":
            assert cell["brute"] == 2187 and cell["n"] == 14
        if spec.name == "copy-blowup-n5-m12":
            assert cell["n_after"] == 36 and cell["brute"] == 8
            assert cell["assignments"] == 2**36 and cell["seconds"] is None
        if spec.name == "copy-random-n4-m8":
            assert cell["n_after"] == 24 and cell["seconds"] is None
        if spec.name == "matrix-full-m0":
            assert cell["side_log2"] == 384
            assert cell["seconds"] is None and cell["seconds_brute"] is None
        if spec.name == "degree-bounded-n6-m4":
            assert cell["brute"] == 44 and cell["observed_degree"] == 12
        if spec.name == "hermite-single-point-not-run":
            assert cell["jet_width"] == 3 and cell["assignments"] == 4
            assert cell["seconds"] is None
