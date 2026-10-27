"""Step 7: a higher slice cap, two new points, and a balance-check scaffold.

The cap is 256. Full M0 at R = 32, Q = 192 is still refused before allocation.
The scaffold refuses when the Lemma 4.5 inequalities fail, and it does not
choose schedules when they hold. These tests do not claim Theorem 1.1, a
refutation of ETH, or a 2^{o(n)} count.
"""

from __future__ import annotations

import ast
import inspect
from pathlib import Path

from fractions import Fraction

import pytest

from counting_sat_algebra.cnf import CNF
from counting_sat_algebra.replacement_matrix import (
    FULL_M0_LOG2_DIMENSION_AT_R2,
    FULL_M0_MIN_Q,
    FULL_M0_MIN_R,
    SLICE_MAX_DIMENSION,
    build_auxiliary_slice,
    degree_cap_report,
    groups_at,
    step3_examples,
)
from counting_sat_algebra.replacements import eval_H_at
from counting_sat_algebra.scheduler import (
    check_balance,
    placed_identity_problems,
    schedule_placement,
    scopes_from_cnf,
    step7_examples,
)


def _imported_names() -> set[str]:
    source_path = Path(__file__).resolve().parents[1] / "src" / "counting_sat_algebra" / "scheduler.py"
    tree = ast.parse(source_path.read_text())
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                names.add(node.module.split(".")[0])
            names.update(alias.name for alias in node.names)
    return names


def _units(n: int) -> CNF:
    return CNF(n, [[(i, False)] for i in range(n)])


def _example(name: str):
    found = [example for example in step7_examples(heavy=True) if example.name == name]
    assert len(found) == 1
    return found[0]


def test_cap_is_256_and_the_next_powers_stay_refused() -> None:
    assert SLICE_MAX_DIMENSION == 256
    assert 2 ** FULL_M0_LOG2_DIMENSION_AT_R2 > SLICE_MAX_DIMENSION
    assert FULL_M0_MIN_R == 32 and FULL_M0_MIN_Q == 192
    # r=2, Q=5 has side length 1024. r=3, Q=3 has side length 729.
    with pytest.raises(ValueError, match=r"1024"):
        build_auxiliary_slice(2, 5, 1)
    with pytest.raises(ValueError, match=r"729"):
        build_auxiliary_slice(3, 3, 1)
    with pytest.raises(ValueError, match=r"2\*\*384"):
        build_auxiliary_slice(2, FULL_M0_MIN_Q, FULL_M0_MIN_R)
    with pytest.raises(ValueError, match=r"2\*\*384"):
        build_auxiliary_slice(2, FULL_M0_MIN_Q, FULL_M0_MIN_R, max_dimension=10**18)


def test_scope_statistics_count_a_repeated_variable_once() -> None:
    check = check_balance(1, [(0, 0)], r=2, Q=1, R=1)
    assert (check.n, check.M, check.b, check.w, check.D) == (1, 1, 1, 1, 4)
    cnf = CNF(3, [[(0, False), (1, False), (2, True)]])
    assert scopes_from_cnf(cnf) == ((0, 1, 2),)
    wider = check_balance(cnf.n, scopes_from_cnf(cnf), r=2, Q=3, R=2)
    assert (wider.M, wider.b, wider.w, wider.D) == (3, 1, 3, 64)
    assert wider.mass_bound == 2 * (1 + 3) * 64 * 64 == 32768
    assert wider.variable_bound == 2 * 3 * 2 * 2 == 24
    assert wider.minimum_n == 32768
    assert not wider.admitted
    assert wider.degree_gap_ok is True


def test_small_formulas_fail_the_balance_inequalities() -> None:
    two = check_balance(2, scopes_from_cnf(_units(2)), r=2, Q=1, R=1)
    assert two.mass_bound == 64 and two.variable_bound == 8
    assert not two.mass_ok and not two.variable_ok and two.range_ok
    assert two.minimum_n == 64
    assert any("2(b+w)D^2" in reason for reason in two.reasons)
    assert any("2Qr^2" in reason for reason in two.reasons)
    assert any("order of D^2" in reason for reason in two.reasons)

    # n = 2 D^2, but the (b+w)/b factor still asks for n >= 4 D^2.
    edge = check_balance(32, scopes_from_cnf(_units(32)), r=2, Q=1, R=1)
    assert edge.D == 4 and edge.order_bound == 32
    assert edge.n == edge.order_bound
    assert not edge.mass_ok and edge.variable_ok and edge.range_ok
    assert edge.minimum_n == 64
    assert not any("order of D^2" in reason for reason in edge.reasons)
    assert not edge.admitted

    for example in step3_examples():
        check = check_balance(
            example.placed.cnf.n,
            scopes_from_cnf(example.placed.cnf),
            example.auxiliary.r,
            example.auxiliary.Q,
            example.auxiliary.R,
        )
        assert not check.admitted
        outcome = schedule_placement(example.placed.cnf, auxiliary=example.auxiliary)
        assert outcome.placed is None and not outcome.accepted and not outcome.built_slice
        assert "No placement was invented" in outcome.reason
        assert "balancer" in outcome.reason


def test_admitted_parameters_do_not_choose_a_schedule(monkeypatch: pytest.MonkeyPatch) -> None:
    def _boom(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("build_auxiliary_slice was called")

    monkeypatch.setattr(
        "counting_sat_algebra.scheduler.build_auxiliary_slice",
        _boom,
    )
    cnf = _units(64)
    check = check_balance(cnf.n, scopes_from_cnf(cnf), r=2, Q=1, R=1)
    assert check.admitted
    assert check.mass_bound == 64 and check.variable_bound == 8
    assert check.degree_gap_ok is False
    outcome = schedule_placement(cnf, r=2, Q=1, R=1)
    assert outcome.balance.admitted and outcome.placed is None
    assert not outcome.built_slice and not outcome.accepted
    assert "no schedules were chosen" in outcome.reason
    assert "not implemented" in outcome.reason
    assert "1 <= R < Q" in outcome.reason
    again = schedule_placement(cnf, r=2, Q=1)
    assert again.placed is None and again.balance.admitted
    assert again.balance.degree_gap_ok is None


def test_full_m0_is_refused_before_allocation(monkeypatch: pytest.MonkeyPatch) -> None:
    def _boom(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("build_auxiliary_slice was called")

    monkeypatch.setattr(
        "counting_sat_algebra.scheduler.build_auxiliary_slice",
        _boom,
    )
    cnf = CNF(1, [[(0, False)]])
    outcome = schedule_placement(cnf, r=2, Q=FULL_M0_MIN_Q, R=FULL_M0_MIN_R)
    assert outcome.balance.D == 2**FULL_M0_LOG2_DIMENSION_AT_R2
    assert outcome.allocation_refused and not outcome.built_slice
    assert outcome.placed is None and not outcome.accepted
    assert "not allocated" in outcome.reason
    assert "2**384" in outcome.reason
    assert not outcome.balance.admitted
    supplied = schedule_placement(
        cnf,
        r=2,
        Q=FULL_M0_MIN_Q,
        R=FULL_M0_MIN_R,
        variable_schedules=((0,) * FULL_M0_MIN_Q,),
        clause_schedules=((0,) * FULL_M0_MIN_Q,),
    )
    assert supplied.placed is None and supplied.allocation_refused
    assert "not allocated" in supplied.reason


def test_caller_schedules_run_design_problems_at() -> None:
    for example in step3_examples():
        outcome = schedule_placement(
            example.placed.cnf,
            auxiliary=example.auxiliary,
            variable_schedules=example.placed.variable_schedules,
            clause_schedules=example.placed.clause_schedules,
        )
        assert outcome.accepted
        assert outcome.caller_supplied and not outcome.built_slice
        assert not outcome.balance.admitted
        assert outcome.design_problems == ()
        assert outcome.restoration_problems == ()
        assert outcome.group_size_ok is True
        assert "not chosen by a balancer" in outcome.reason
        assert "balance inequalities fail" in outcome.reason
    swap = next(example for example in step3_examples() if example.name == "swap-4")
    remark = next(example for example in step3_examples() if example.name == "remark-4-4")
    swap_out = schedule_placement(
        swap.placed.cnf,
        r=2,
        Q=1,
        R=1,
        variable_schedules=((0,), (1,)),
        clause_schedules=((1,), (0,)),
    )
    assert swap_out.accepted and swap_out.built_slice
    assert swap_out.multiplicity_ok is True
    assert swap_out.multiplicity_bound == Fraction(1)
    remark_out = schedule_placement(
        remark.placed.cnf,
        auxiliary=remark.auxiliary,
        variable_schedules=remark.placed.variable_schedules,
        clause_schedules=remark.placed.clause_schedules,
    )
    assert remark_out.accepted and remark_out.multiplicity_ok is False
    assert remark_out.multiplicity_bound == Fraction(1, 4)

    crowded = build_auxiliary_slice(3, 1, 1)
    refused = schedule_placement(
        CNF(3, [[(0, False), (1, False), (2, False)]]),
        auxiliary=crowded,
        variable_schedules=((0,), (0,), (0,)),
        clause_schedules=((0,),),
    )
    assert not refused.accepted and refused.group_size_ok is False
    assert any("size" in problem for _t, problem in refused.design_problems)


def test_r3_q2_matches_the_step3_checks() -> None:
    example = next(row for row in step7_examples(heavy=False) if row.name == "r3-q2")
    sl = example.auxiliary
    assert (sl.r, sl.Q, sl.R, sl.dimension, sl.A0, sl.interval) == (3, 2, 1, 81, 0, (0,))
    assert sl.R < sl.Q and sl.used_default_interval
    assert len(sl.C_rows) == 45
    assert placed_identity_problems(example.placed) == []
    placed = example.placed
    assert groups_at(placed, 1) == ((0, 1), (2,), ())
    assert placed.clause_schedules == ((0, 0),)
    assert 2 not in groups_at(placed, 1)[0]
    report_values = {
        Fraction(0): Fraction(7),
        Fraction(1): Fraction(7),
        Fraction(2): Fraction(7),
    }
    for t, expected in report_values.items():
        assert eval_H_at(placed.cnf, placed.family, t) == expected
    outcome = schedule_placement(
        placed.cnf,
        auxiliary=sl,
        variable_schedules=placed.variable_schedules,
        clause_schedules=placed.clause_schedules,
    )
    assert outcome.accepted and not outcome.balance.admitted
    assert outcome.multiplicity == 2
    assert outcome.multiplicity_bound == Fraction(6, 81)
    assert outcome.multiplicity_ok is False
    assert outcome.group_size_ok is True


def test_r2_q4_matches_the_step3_checks() -> None:
    example = _example("r2-q4")
    sl = example.auxiliary
    assert (sl.r, sl.Q, sl.R, sl.dimension, sl.A0, sl.interval) == (
        2,
        4,
        3,
        256,
        0,
        (0, 1, 2),
    )
    assert sl.R < sl.Q and sl.used_default_interval
    assert len(sl.C_rows) == 13
    assert placed_identity_problems(example.placed) == []
    placed = example.placed
    assert groups_at(placed, 4) == ((0, 1), (2,))
    assert placed.clause_schedules[0][3] == 0
    expected = {
        Fraction(0): Fraction(7),
        Fraction(1): Fraction(29, 4),
        Fraction(2): Fraction(15, 2),
        Fraction(3): Fraction(31, 4),
        Fraction(4): Fraction(8),
    }
    for t, value in expected.items():
        assert eval_H_at(placed.cnf, placed.family, t) == value
    outcome = schedule_placement(
        placed.cnf,
        auxiliary=sl,
        variable_schedules=placed.variable_schedules,
        clause_schedules=placed.clause_schedules,
    )
    assert outcome.accepted and not outcome.balance.admitted
    assert outcome.balance.degree_gap_ok is True
    assert outcome.multiplicity == 2
    assert outcome.multiplicity_ok is False
    report = degree_cap_report(placed)
    assert report.ok
    assert report.brute == 7 and report.observed_degree == 1
    assert report.cap_occurrences == 6 and report.cap_step1 == 9


def test_scheduler_does_not_sample_time_or_claim_the_theorem() -> None:
    imported = _imported_names()
    assert "random" not in imported
    assert "time" not in imported
    assert "brute_count" not in imported
    source = inspect.getsource(schedule_placement)
    assert "import random" not in source and "Random(" not in source
    text = source.lower()
    assert "eth is false" not in text
    assert "theorem 1.1" not in text
    with pytest.raises(ValueError, match="together"):
        schedule_placement(_units(1), r=2, Q=1, R=1, variable_schedules=((0,),))
    with pytest.raises(ValueError, match="bool|int"):
        check_balance(1, [(True,)], r=2, Q=1)  # type: ignore[list-item]
