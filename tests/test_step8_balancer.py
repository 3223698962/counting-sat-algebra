"""Step 8: a machine-limited successive-choice balancer.

The search runs only when the Lemma 4.5 inequalities hold and the slice is
under the cap. It does not enumerate the schedule cube and it does not sample
a placement. The admitted toy meets the multiplicity and group bounds and is
refused before ``design_problems_at``, because that walk would be 2**32
Boolean points. These tests do not claim Theorem 1.1, a refutation of ETH, or
a 2^{o(n)} count.
"""

from __future__ import annotations

import ast
import inspect
from dataclasses import replace
from fractions import Fraction
from functools import partial

import pytest

from counting_sat_algebra.cnf import CNF
from counting_sat_algebra.replacement_matrix import (
    FULL_M0_MIN_Q,
    FULL_M0_MIN_R,
    build_auxiliary_slice,
)
from counting_sat_algebra.scheduler import (
    BALANCER_DESIGN_BOOLEAN_BUDGET,
    BALANCER_MAX_ALPHABET,
    BALANCER_MAX_TRIALS,
    _finish_chosen_placement,
    _search_block_reason,
    _search_schedules,
    check_balance,
    schedule_placement as _schedule_placement,
    scopes_from_cnf,
)


schedule_placement = partial(_schedule_placement, design_method="boolean")

def _units(n: int) -> CNF:
    return CNF(n, [[(i, False)] for i in range(n)])


def _boom(*_args: object, **_kwargs: object) -> None:
    raise AssertionError("build_auxiliary_slice was called")


def test_balance_failures_do_not_search(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("counting_sat_algebra.scheduler.build_auxiliary_slice", _boom)
    monkeypatch.setattr("counting_sat_algebra.scheduler._search_schedules", _boom)
    for cnf, kwargs in (
        (_units(2), {"r": 2, "Q": 1, "R": 1}),
        (_units(8), {"r": 2, "Q": 2, "R": 1}),
        (CNF(3, [[(0, False), (1, False), (2, True)]]), {"r": 2, "Q": 1, "R": 1}),
    ):
        outcome = schedule_placement(cnf, **kwargs)
        assert outcome.placed is None and not outcome.accepted
        assert not outcome.balancer_attempted and not outcome.balancer_built
        assert not outcome.built_slice and outcome.multiplicity is None
        assert "No placement was invented" in outcome.reason
        assert "was not run" in outcome.reason
        assert not outcome.balance.admitted


def test_admitted_toy_refuses_on_the_design_budget(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("counting_sat_algebra.scheduler.build_auxiliary_slice", _boom)
    cnf = _units(64)
    first = schedule_placement(cnf, r=2, Q=1, R=1)
    second = schedule_placement(cnf, r=2, Q=1, R=1)
    assert first.reason == second.reason
    assert first.multiplicity == 16 and first.multiplicity_bound == Fraction(32)
    assert first.trials == (64 + 64) * 2
    assert first.outside_group_size == 32
    assert "schedule cube was not enumerated" in first.reason
    search = _search_schedules(cnf, first.balance)
    assert search.variable_schedules == _search_schedules(cnf, first.balance).variable_schedules
    assert search.clause_schedules is not None and search.variable_schedules is not None
    assert len(set(search.variable_schedules)) == 2
    assert len(set(search.clause_schedules)) == 2
    assert search.chosen and search.design_boolean_points is None


def test_gap_ok_admitted_instance_also_refuses_the_design_budget(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr("counting_sat_algebra.scheduler.build_auxiliary_slice", _boom)
    cnf = _units(1024)
    outcome = schedule_placement(cnf, r=2, Q=2, R=1)
    assert outcome.balance.admitted and outcome.balance.degree_gap_ok is True
    assert outcome.balance.D == 16
    assert outcome.placed is None and not outcome.built_slice and not outcome.accepted
    assert outcome.balancer_attempted and not outcome.balancer_built
    assert outcome.multiplicity == 64 and outcome.multiplicity_bound == Fraction(128)
    assert outcome.multiplicity_ok is True and outcome.group_size_ok is True
    assert outcome.max_group_size == 512 and outcome.outside_group_size == 512
    assert outcome.trials == (1024 + 1024) * 4
    assert "2**512" in outcome.reason
    assert "1 <= R < Q" not in outcome.reason


def test_full_m0_and_the_next_power_do_not_search(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("counting_sat_algebra.scheduler._search_schedules", _boom)
    monkeypatch.setattr("counting_sat_algebra.scheduler.build_auxiliary_slice", _boom)
    full = schedule_placement(_units(1), r=2, Q=FULL_M0_MIN_Q, R=FULL_M0_MIN_R)
    assert full.allocation_refused and full.placed is None
    assert not full.balancer_attempted
    assert "not allocated" in full.reason and "2**384" in full.reason
    wide = schedule_placement(_units(1), r=2, Q=5, R=1)
    assert wide.balance.D == 1024 and wide.allocation_refused
    assert wide.placed is None and not wide.balancer_attempted
    assert "not allocated" in wide.reason
    with pytest.raises(ValueError, match=r"1024"):
        build_auxiliary_slice(2, 5, 1)


def test_direct_finish_on_an_unadmitted_formula_passes_design() -> None:
    """The checker wiring. ``schedule_placement`` does not take this branch.

    Eight units at r = 2, Q = 2 fail the balance inequalities. The same
    schedules, placed by a direct call, pass ``design_problems_at`` and the
    T = 0 restoration. That call is not a Lemma 4.5 witness.
    """
    cnf = _units(8)
    check = check_balance(cnf.n, scopes_from_cnf(cnf), r=2, Q=2, R=1)
    assert not check.admitted
    assert check.degree_gap_ok is True
    assert check.mass_bound == 2 * (1 + 1) * 16 * 16 == 1024
    search = _search_schedules(cnf, check)
    assert search.chosen and search.multiplicity == 1
    assert search.design_boolean_points == 256
    assert search.design_boolean_points <= BALANCER_DESIGN_BOOLEAN_BUDGET
    finished = _finish_chosen_placement(cnf, check, search, None)
    assert finished.accepted and finished.balancer_built and finished.design_checked
    assert finished.caller_supplied is False and finished.built_slice
    assert finished.placed is not None and finished.placed.copies == 1
    assert finished.design_problems == () and finished.restoration_problems == ()
    assert "were not the gate" in finished.reason
    assert "not supplied by the caller" in finished.reason
    refused = schedule_placement(cnf, r=2, Q=2, R=1)
    assert refused.placed is None and not refused.balancer_attempted
    assert "was not run" in refused.reason


def test_forced_admission_reaches_the_checker_inside_schedule_placement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exercise the success return. The real inequalities on this formula fail.

    The flag is forced so the gate's true branch runs. The numbers below are
    the real scope statistics, not a second balance proof.
    """
    real = check_balance

    def _force(n: int, scopes: object, r: int, Q: int, R: int | None = None):
        check = real(n, scopes, r, Q, R)  # type: ignore[arg-type]
        return replace(check, admitted=True, reasons=())

    monkeypatch.setattr("counting_sat_algebra.scheduler.check_balance", _force)
    cnf = _units(8)
    outcome = schedule_placement(cnf, r=2, Q=2, R=1)
    assert outcome.accepted and outcome.balancer_built and not outcome.caller_supplied
    assert outcome.design_checked and outcome.built_slice
    assert outcome.design_problems == () and outcome.restoration_problems == ()
    assert outcome.multiplicity == 1 and outcome.group_size_ok is True
    assert "were chosen by the successive-choice balancer" in outcome.reason
    monkeypatch.setattr(
        "counting_sat_algebra.scheduler.BALANCER_DESIGN_BOOLEAN_BUDGET",
        0,
    )
    blocked = schedule_placement(cnf, r=2, Q=2, R=1)
    assert blocked.placed is None and not blocked.balancer_built
    assert blocked.balancer_attempted and "budget" in blocked.reason
    assert "place_formula was not called" in blocked.reason


def test_multiplicity_group_alphabet_and_trial_refusals(monkeypatch: pytest.MonkeyPatch) -> None:
    triple = CNF(1, [[(0, False)], [(0, False)], [(0, False)]])
    check = check_balance(triple.n, scopes_from_cnf(triple), r=2, Q=1, R=1)
    assert not check.admitted
    search = _search_schedules(triple, check)
    assert search.attempted and not search.chosen
    assert search.multiplicity == 2 and search.multiplicity_bound == Fraction(3, 2)
    assert search.multiplicity_ok is False
    assert "exceeds 2M/D" in search.reason
    assert "No placement was invented" in search.reason
    refused = schedule_placement(triple, r=2, Q=1, R=1)
    assert refused.placed is None and not refused.balancer_attempted

    tight = check_balance(2, scopes_from_cnf(_units(2)), r=8, Q=1, R=1)
    grouped = _search_schedules(_units(2), tight)
    assert grouped.attempted and not grouped.chosen and grouped.group_size_ok is False
    assert "group" in grouped.reason and "No placement was invented" in grouped.reason

    monkeypatch.setattr("counting_sat_algebra.scheduler.BALANCER_MAX_ALPHABET", 1)
    alphabet = _search_schedules(_units(1), check_balance(1, [(0,)], r=2, Q=1))
    assert not alphabet.attempted and alphabet.variable_schedules is None
    assert "alphabet" in alphabet.reason and "schedule cube was not enumerated" in alphabet.reason
    monkeypatch.setattr(
        "counting_sat_algebra.scheduler.BALANCER_MAX_ALPHABET",
        BALANCER_MAX_ALPHABET,
    )
    monkeypatch.setattr("counting_sat_algebra.scheduler.BALANCER_MAX_TRIALS", 4)
    trials = _search_schedules(_units(3), check_balance(3, scopes_from_cnf(_units(3)), r=2, Q=1))
    assert not trials.attempted and "100000" not in trials.reason
    assert "4" in trials.reason and "schedule cube was not enumerated" in trials.reason
    assert _search_block_reason(100_000, 0, 2, 1) is not None
    assert _search_block_reason(1, 1, 2, 1) is None


def test_search_enumerates_only_the_alphabet() -> None:
    from counting_sat_algebra.conditional_balance import choose_schedules
    source = inspect.getsource(choose_schedules)
    assert source.count("product(") == 1
    assert "repeat=Q" in source
    assert "import random" not in source and "Random(" not in source
    tree = ast.parse(
        (
            __import__("pathlib").Path(__file__).resolve().parents[1]
            / "src"
            / "counting_sat_algebra"
            / "scheduler.py"
        ).read_text()
    )
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert "random" not in imported and "time" not in imported
