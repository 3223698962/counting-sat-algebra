#!/usr/bin/env python3
"""
Verification CLI: Step 0 (product identity + H(T) interpolation), a small
Step 1 bounded-occurrence degree check, Step 2 blocked-sum agreement, a
Step 3 tiny replacement-matrix slice, Step 4 Hermite recovery of H(0),
Step 5 occurrence copying with #F' = #F on instances inside the brute window,
and a Step 6 parameter scan that records matches and refusals.

Usage:
  python -m counting_sat_algebra.verify
  python -m counting_sat_algebra.verify --quick
  python -m counting_sat_algebra.verify --step2
  python -m counting_sat_algebra.verify --step3
  python -m counting_sat_algebra.verify --step4
  python -m counting_sat_algebra.verify --step5
  python -m counting_sat_algebra.verify --step6
  python -m counting_sat_algebra.verify --step7
  python -m counting_sat_algebra.verify --step8
"""

from __future__ import annotations

import argparse
import random
import sys
import time
from fractions import Fraction
from typing import List, Sequence, Tuple

from .blocked_sum import (
    blocked_sum,
    compare_blocked_to_full,
    hand_built_examples,
    popcount_group_instance,
    popcount_growth,
)
from .bounded_occurrence import random_bounded_3cnf, run_degree_bound_experiment
from .cnf import brute_count, random_3cnf
from .hermite import analyze_step4_case, recovery_problems, step4_cases
from .independent import run_independent_suite
from .replacement_matrix import (
    FULL_M0_LOG2_DIMENSION_AT_R2,
    FULL_M0_MIN_Q,
    FULL_M0_MIN_R,
    SLICE_MAX_DIMENSION,
    auxiliary_problems,
    build_auxiliary_slice,
    degree_cap_report,
    design_problems_at,
    restoration_problems,
    step3_examples,
)
from .scheduler import (
    BALANCER_DESIGN_BOOLEAN_BUDGET,
    _finish_chosen_placement,
    _search_schedules,
    check_balance,
    placed_identity_problems,
    schedule_placement,
    scopes_from_cnf,
    step7_examples,
)
from .parameter_scan import TIMING_DISCLAIMER, format_distance, format_scan_line, run_parameter_scan
from .replacements import eval_H_at, interpolate_H0
from .self_reduce import find_assignment_by_self_reduction
from .sparsify import Step5Report, analyze_step5_case, step5_cases


def run_interpolation_suite(
    trials: int = 50,
    seed: int = 42,
    R: int = 1,
    n_min: int = 4,
    n_max: int = 10,
) -> List[dict]:
    """
    Random 3-CNF: n in [n_min, n_max], m in [n, 2n].
    Compare brute #F with interpolated H(0). R=1 keeps the degree cap modest.
    """
    rng = random.Random(seed)
    rows: List[dict] = []
    for i in range(trials):
        n = rng.randint(n_min, n_max)
        m = rng.randint(n, 2 * n)
        cnf = random_3cnf(n, m, rng=rng)
        t0 = time.perf_counter()
        brute = brute_count(cnf)
        h0, D, _fam = interpolate_H0(cnf, R=R, rng=rng, k=3)
        dt = time.perf_counter() - t0
        ok = (h0 == Fraction(brute)) and h0.denominator == 1
        rows.append(
            {
                "trial": i,
                "n": n,
                "m": m,
                "brute": brute,
                "H0": h0,
                "D": D,
                "ok": ok,
                "sec": dt,
            }
        )
        status = "OK" if ok else "FAIL"
        print(
            f"  [B {i+1:02d}/{trials}] n={n:2d} m={m:2d}  "
            f"#F={brute:5d}  H(0)={h0}  deg≤{D}  {status}  ({dt:.2f}s)",
            flush=True,
        )
        if not ok:
            break
    return rows


def run_self_reduction_demo(seed: int = 7) -> None:
    rng = random.Random(seed)
    for _ in range(30):
        cnf = random_3cnf(6, 8, rng=rng)
        if brute_count(cnf) > 0:
            break
    else:
        print("  [self-reduce] skipped (no satisfiable sample)")
        return
    asg = find_assignment_by_self_reduction(cnf)
    print(f"  [self-reduce] n={cnf.n} m={cnf.m} #F={brute_count(cnf)}")
    print(f"  [self-reduce] assignment={asg}")
    if asg is None:
        print("FAIL: self-reduction returned no assignment", file=sys.stderr)
        raise SystemExit(1)
    print("  [self-reduce] OK")


def run_bounded_suite(
    specs: Sequence[Tuple[int, int, int]],
    seed: int,
) -> List[dict]:
    """
    specs: (n, m, R) bounded-occurrence 3-CNFs.
    Interpolate at D = m*k*R and at D = 3*n*R; both must recover #F.
    """
    rows: List[dict] = []
    for i, (n, m, R) in enumerate(specs):
        rng = random.Random(seed + i)
        cnf = random_bounded_3cnf(n, m, rng=rng)
        t0 = time.perf_counter()
        result = run_degree_bound_experiment(cnf, R=R, rng=random.Random(seed + 1000 + i))
        dt = time.perf_counter() - t0
        ok = bool(result["ok_mkR"] and result["ok_3nR"])
        result["ok"] = ok
        result["sec"] = dt
        rows.append(result)
        status = "OK" if ok else "FAIL"
        print(
            f"  [C {i+1:02d}/{len(specs)}] n={n} m={m} R={R}  "
            f"occ≤{result['max_occurrence']}  #F={result['brute']}  "
            f"H(0)={result['H0_mkR']}  deg_obs={result['observed_degree']}  "
            f"D={result['D_mkR']}/{result['D_3nR']}  {status}  ({dt:.2f}s)",
            flush=True,
        )
        if not ok:
            break
    return rows


def run_step2_suite(group_counts: Sequence[int]) -> tuple[list[dict], list[dict]]:
    """Hand-built blocked-sum checks, then a popcount scan over the group count.

    Returns ``(example_rows, growth_rows)``. A row with ``ok`` false stops the
    scan. Measurement-tuple counts for the popcount scan are ``3^b`` with
    block size 2; that is exponential in the number of groups.
    """
    print("=== (D) Blocked sum vs full H(t), hand-built replacements ===")
    print(
        "  Simplified (r, τ, D)-design: blocks of size ≤ r and τ shared linear measurements."
    )
    print(
        "  Measurement tuples can grow exponentially with the number of blocks. "
        "Theorem 1.1 is not implemented. This is not a 2^{o(n)} algorithm."
    )
    example_rows: list[dict] = []
    examples = hand_built_examples()
    for index, example in enumerate(examples, start=1):
        started = time.perf_counter()
        comparisons = compare_blocked_to_full(
            example.cnf,
            example.structured,
            example.sample_points,
        )
        brute = brute_count(example.cnf)
        ok = all(blocked.value == full for _t, full, blocked in comparisons)
        at_zero = [blocked.value for t, _full, blocked in comparisons if t == 0]
        ok = ok and len(at_zero) == 1 and at_zero[0] == Fraction(brute)
        elapsed = time.perf_counter() - started
        design = example.structured.design
        sample = comparisons[0][2]
        row = {
            "name": example.name,
            "n": example.cnf.n,
            "m": example.cnf.m,
            "blocks": len(design.blocks),
            "r": design.r,
            "tau": design.tau,
            "R": example.structured.R,
            "D": design.D,
            "brute": brute,
            "tuples": sample.measurement_tuples,
            "local": sample.local_evaluations,
            "full_space": sample.full_space,
            "ok": ok,
            "sec": elapsed,
            "points": len(comparisons),
        }
        example_rows.append(row)
        status = "OK" if ok else "FAIL"
        print(
            f"  [D {index:02d}/{len(examples)}] {example.name}  "
            f"n={row['n']} blocks={row['blocks']} r={row['r']} τ={row['tau']} "
            f"R={row['R']} D={row['D']}  #F={brute}  "
            f"tuples={row['tuples']} local={row['local']} full={row['full_space']}  "
            f"points={row['points']}  {status}  ({elapsed:.3f}s)",
            flush=True,
        )
        if not ok:
            for t, full, blocked in comparisons:
                if blocked.value != full or (t == 0 and blocked.value != Fraction(brute)):
                    print(
                        f"    mismatch at t={t}: blocked={blocked.value} full={full} #F={brute}",
                        flush=True,
                    )
            return example_rows, []

    print("  Popcount growth, block size 2, every α coefficient 1, compared at t = 1.")
    print("  Measurement tuples = (r+1)^b = 3^b, exponential in the number of groups.")
    growth_rows: list[dict] = []
    for growth in popcount_growth(group_counts, block_size=2, R=1, t=1):
        row = {
            "n": growth.blocked.n,
            "m": growth.groups,
            "groups": growth.groups,
            "tuples": growth.blocked.measurement_tuples,
            "local": growth.blocked.local_evaluations,
            "full_space": growth.blocked.full_space,
            "seconds": growth.blocked.seconds,
            "ok": growth.ok,
            "value": growth.blocked.value,
            "full": growth.full,
        }
        growth_rows.append(row)
        status = "OK" if growth.ok else "FAIL"
        print(
            f"    b={growth.groups} n={growth.blocked.n}  "
            f"H(1)={growth.blocked.value}  tuples={growth.blocked.measurement_tuples}  "
            f"local_evals={growth.blocked.local_evaluations}  "
            f"full_space={growth.blocked.full_space}  "
            f"enum={growth.blocked.seconds:.4f}s  {status}",
            flush=True,
        )
        if not growth.ok:
            print(f"    full H(1)={growth.full}", flush=True)
            return example_rows, growth_rows

    largest = group_counts[-1]
    cnf, structured = popcount_group_instance(largest, block_size=2, R=1)
    blocked_at_zero = blocked_sum(cnf, structured, 0)
    brute = brute_count(cnf)
    full_at_zero = eval_H_at(cnf, structured.family, Fraction(0))
    ok_zero = blocked_at_zero.value == Fraction(brute) == full_at_zero
    print(
        f"    b={largest} at t=0: blocked={blocked_at_zero.value}  #F={brute}  "
        f"{'OK' if ok_zero else 'FAIL'}",
        flush=True,
    )
    if not ok_zero:
        growth_rows[-1]["ok"] = False
    return example_rows, growth_rows


def run_step3_suite() -> list[dict]:
    """Tiny M0 slices: restoration, degree caps, and a legal design at every t.

    The Step 3 examples have side length at most 64. The slice cap is
    ``SLICE_MAX_DIMENSION``. Points above 64 are Step 7. The parameter list
    with R ≥ 32 and 6R ≤ Q is not built. Schedules are the hand-written ones
    in ``step3_examples``.
    """
    print("=== (E) Replacement matrix M0 at tiny parameters ===")
    print(
        "  M0 = V0 f_I(W0), side length D = r**(2Q). "
        "The Step 3 examples have side length at most 64. "
        f"The slice cap is {SLICE_MAX_DIMENSION}; larger checked points are Step 7."
    )
    print(
        "  The parameter list with R≥32 and 6R≤Q is not built "
        f"(at r=2 that side length is 2**{FULL_M0_LOG2_DIMENSION_AT_R2}). "
        "Schedules are written out by hand. "
        "Theorem 1.1 is not implemented. This is not a 2^{o(n)} algorithm."
    )
    rows: list[dict] = []
    examples = step3_examples()
    for index, example in enumerate(examples, start=1):
        started = time.perf_counter()
        matrix_problems = auxiliary_problems(example.auxiliary)
        restored = restoration_problems(example.placed)
        design: list[str] = []
        for t in range(1, example.auxiliary.Q + 1):
            design.extend(design_problems_at(example.placed, t))
        report = degree_cap_report(example.placed)
        elapsed = time.perf_counter() - started
        ok = not matrix_problems and not restored and not design and report.ok
        sl = example.auxiliary
        row = {
            "name": example.name,
            "n": example.placed.cnf.n,
            "m": example.placed.cnf.m,
            "r": sl.r,
            "Q": sl.Q,
            "R": sl.R,
            "D": sl.dimension,
            "measurements": len(sl.C_rows),
            "brute": report.brute,
            "degree": report.observed_degree,
            "cap_m": report.cap_occurrences,
            "cap_step1": report.cap_step1,
            "ok": ok,
            "sec": elapsed,
        }
        rows.append(row)
        status = "OK" if ok else "FAIL"
        print(
            f"  [E {index:02d}/{len(examples)}] {example.name}  "
            f"r={sl.r} Q={sl.Q} R={sl.R} D={sl.dimension}  "
            f"rows={len(sl.C_rows)}  n={row['n']} m={row['m']}  "
            f"#F={report.brute}  deg={report.observed_degree}  "
            f"caps={report.cap_occurrences}/{report.cap_step1}  "
            f"{status}  ({elapsed:.3f}s)",
            flush=True,
        )
        if not ok:
            for problem in matrix_problems + restored + design:
                print(f"    {problem}", flush=True)
            if not report.ok:
                print(
                    f"    degree report H(0)={report.direct_at_zero} "
                    f"interpolated={report.h0_occurrences} #F={report.brute} "
                    f"deg={report.observed_degree}",
                    flush=True,
                )
            return rows
    return rows


def run_step4_suite() -> list[dict]:
    """Hermite recovery of H(0) from Hasse jets, compared with Lagrange.

    The condition count stays D+1. Distinct sample points are fewer than the
    value-only path. Each point is one pass over {0,1}^n. A single node that
    carries every coefficient does quadratic arithmetic in D per assignment.
    """
    print("=== (F) Hermite / Hasse interpolation ===")
    print(
        "  Several Hasse coefficients of H at each sample t, from the formal "
        "derivative of the clause product."
    )
    print(
        "  The condition count stays D+1, packed into fewer distinct points than "
        "Lagrange on T = 1..D+1."
    )
    print(
        "  Each distinct t is one pass over 2^n assignments. Derivative orders "
        "are not separate passes. A full jet of order D does quadratic arithmetic "
        "per assignment."
    )
    print(
        "  Theorem 1.1 is not implemented. This is not a 2^{o(n)} algorithm."
    )
    rows: list[dict] = []
    cases = step4_cases()
    for index, case in enumerate(cases, start=1):
        started = time.perf_counter()
        recoveries = analyze_step4_case(case)
        elapsed = time.perf_counter() - started
        problems: list[str] = []
        for recovery in recoveries:
            problems.extend(
                f"{recovery.scheme}: {problem}"
                for problem in recovery_problems(recovery)
            )
        if not recoveries:
            rows.append(
                {
                    "name": case.name,
                    "n": case.cnf.n,
                    "m": case.cnf.m,
                    "brute": -1,
                    "degree_cap": -1,
                    "observed_degree": -1,
                    "schemes": 0,
                    "ok": False,
                    "sec": time.perf_counter() - started,
                }
            )
            print(f"  [F {index:02d}/{len(cases)}] {case.name}  no schemes  FAIL", flush=True)
            return rows
        ok = not problems
        head = recoveries[0]
        row = {
            "name": case.name,
            "n": head.n,
            "m": head.m,
            "brute": head.brute,
            "degree_cap": head.degree_cap,
            "observed_degree": head.observed_degree,
            "schemes": len(recoveries),
            "ok": ok,
            "sec": elapsed,
        }
        rows.append(row)
        status = "OK" if ok else "FAIL"
        print(
            f"  [F {index:02d}/{len(cases)}] {case.name}  "
            f"n={head.n} m={head.m} R={head.R} k={head.k} D={head.degree_cap}  "
            f"#F={head.brute}  deg={head.observed_degree}  "
            f"schemes={len(recoveries)}  {status}  ({elapsed:.3f}s)",
            flush=True,
        )
        for recovery in recoveries:
            print(
                f"      {recovery.scheme:<12} {recovery.sample_comparison()}  "
                f"passes={recovery.cube_passes}  assignments={recovery.assignments}  "
                f"H(0)={recovery.h0}",
                flush=True,
            )
        if not ok:
            for problem in problems:
                print(f"    {problem}", flush=True)
            return rows
    return rows


def run_step5_suite() -> list[dict]:
    """Copy high-occurrence variables, then check #F' and the downstream steps.

    The recovery on a counted instance is #F' = #F. A reduced n above 20, an
    interpolation plan above the assignment budget, or a Step 2 clause that
    meets two blocks is printed as REFUSED and does not fail the suite when
    that refusal is the one the case records. Step 3 is not invoked.
    """
    print("=== (G) Copy variables down to occurrence ≤ 3 ===")
    print("  Implication cycle on the copies. On a counted instance, #F' = #F.")
    print(
        "  Reduced n above 20, an interpolation budget overrun, and a clause that "
        "meets two Step 2 blocks are recorded and are not treated as counts."
    )
    print("  Theorem 1.1 is not implemented. This is not a 2^{o(n)} algorithm.")
    rows: list[dict] = []
    cases = step5_cases()
    for index, case in enumerate(cases, start=1):
        report = analyze_step5_case(case)
        row = _step5_row(report)
        rows.append(row)
        print(_format_step5_line(index, len(cases), report, row), flush=True)
        if not report.ok:
            for problem in report.problems:
                print(f"    {problem}", flush=True)
            return rows
    return rows


def run_step6_suite(profile: str) -> list[dict]:
    """Parameter scan. ``quick`` is the small grid; ``full`` adds larger n.

    A refused cell is a passing row when it is the refusal the grid records.
    Wall times are printed and are not a complexity result. Full M0 is not built.
    """
    print("=== (H) Parameter scan ===")
    print(f"  Profile {profile}.")
    print(f"  {TIMING_DISCLAIMER}")
    print(
        "  Cells past the brute-force wall, the assignment budget, or the "
        "measurement-tuple cap are refused before that enumeration."
    )
    print("  Theorem 1.1 is not implemented. This is not a 2^{o(n)} algorithm.")
    report = run_parameter_scan(profile)
    rows: list[dict] = []
    total = len(report.cells)
    for index, cell in enumerate(report.cells, start=1):
        row = _step6_row(cell)
        rows.append(row)
        print(format_scan_line(index, total, cell), flush=True)
        if not cell.ok:
            for problem in cell.problems:
                print(f"    {problem}", flush=True)
            print()
            print(format_distance(report.distance))
            return rows
    print()
    print(format_distance(report.distance))
    return rows


def _step6_row(cell) -> dict:
    if cell.status == "agreed":
        label = "AGREED"
    elif cell.status == "refused":
        label = "REFUSED"
    else:
        label = "FAIL"
    return {
        "name": cell.name,
        "n": -1 if cell.n is None else cell.n,
        "brute": -1 if cell.brute is None else cell.brute,
        "label": label,
        "ok": cell.ok,
        "status": cell.status,
    }


def _step5_row(report: Step5Report) -> dict:
    spec = report.sparsification
    brute = report.count.brute_reduced
    label = "FAIL"
    if report.ok:
        label = "RECORDED" if report.count.status == "refused" else "OK"
    return {
        "name": report.name,
        "n": spec.n_after,
        "m": spec.m_after,
        "n_before": spec.n_before,
        "m_before": spec.m_before,
        "occ_before": spec.max_occurrence_before,
        "occ_after": spec.max_occurrence_after,
        "brute": brute,
        "label": label,
        "ok": report.ok,
        "sec": report.seconds,
        "count_status": report.count.status,
        "pipeline_status": None if report.pipeline is None else report.pipeline.status,
        "blocked_status": None if report.blocked is None else report.blocked.status,
    }


def _format_step5_line(index: int, total: int, report: Step5Report, row: dict) -> str:
    spec = report.sparsification
    if report.count.status == "agreed":
        count_text = str(report.count.brute_reduced)
    else:
        count_text = "REFUSED"
    if report.pipeline is None:
        if report.blocked is not None and report.blocked.status == "agreed":
            h0_text = f"{report.blocked.blocked_at_zero} via blocked sum"
        else:
            h0_text = "—"
    elif report.pipeline.status == "agreed":
        degree = report.pipeline.observed_degree
        mass = report.pipeline.degree_cap_mass
        h0_text = f"{report.pipeline.h0} deg={degree}/{mass}"
    elif report.pipeline.status == "refused":
        h0_text = "REFUSED"
    else:
        h0_text = "FAIL"
    if report.blocked is None:
        block_text = "—"
    elif report.blocked.status == "agreed":
        block_text = f"OK tuples={report.blocked.tuples}"
    elif report.blocked.status == "refused":
        block_text = "REFUSED"
    else:
        block_text = "FAIL"
    return (
        f"  [G {index:02d}/{total}] {report.name}  "
        f"n {spec.n_before}→{spec.n_after}  m {spec.m_before}→{spec.m_after}  "
        f"occ {spec.max_occurrence_before}→{spec.max_occurrence_after}  "
        f"#F={count_text}  H(0)={h0_text}  blocks={block_text}  "
        f"{row['label']}  ({report.seconds:.3f}s)"
    )


def _units(n: int):
    from .cnf import CNF

    return CNF(n, [[(i, False)] for i in range(n)])


def run_step7_suite(*, heavy: bool) -> list[dict]:
    """Higher slice cap, explicit new points, and balance refusals.

    ``heavy=False`` skips the side-length-256 identity. Refusals are not timed.
    The clock on an identity that did run is not a complexity result.
    """
    print("=== (I) Slice cap, new points, and the balance scaffold ===")
    print(
        f"  Slice cap D ≤ {SLICE_MAX_DIMENSION}. "
        "The side-length-256 point runs only outside --quick."
        if heavy
        else f"  Slice cap D ≤ {SLICE_MAX_DIMENSION}. This quick run skips side length 256."
    )
    print(
        "  Balance checks are 2(b+w)D^2 <= M and 2Qr^2 <= n. "
        "When they fail, no schedule is chosen."
    )
    print(
        "  When they hold, Step 8's successive-choice search runs inside a "
        "machine budget. This suite records the admitted point and does not "
        "time the refusal."
    )
    print("  Theorem 1.1 is not implemented. This is not a 2^{o(n)} algorithm.")
    rows: list[dict] = []

    def add(name: str, n: int, m: int, brute: int, ok: bool, label: str) -> None:
        rows.append(
            {
                "name": name,
                "n": n,
                "m": m,
                "brute": brute,
                "ok": ok,
                "label": label,
            }
        )

    examples = step7_examples(heavy=heavy)
    for index, example in enumerate(examples, start=1):
        started = time.perf_counter()
        problems = placed_identity_problems(example.placed)
        outcome = schedule_placement(
            example.placed.cnf,
            auxiliary=example.auxiliary,
            variable_schedules=example.placed.variable_schedules,
            clause_schedules=example.placed.clause_schedules,
        )
        elapsed = time.perf_counter() - started
        ok = not problems and outcome.accepted and not outcome.balance.admitted
        sl = example.auxiliary
        add(example.name, example.placed.cnf.n, example.placed.cnf.m, brute_count(example.placed.cnf), ok, "OK" if ok else "FAIL")
        status = "OK" if ok else "FAIL"
        print(
            f"  [I {index:02d}] {example.name}  "
            f"r={sl.r} Q={sl.Q} R={sl.R} D={sl.dimension}  "
            f"rows={len(sl.C_rows)}  "
            f"{status}  ({elapsed:.3f}s)",
            flush=True,
        )
        if not ok:
            for problem in problems:
                print(f"    {problem}", flush=True)
            if not outcome.accepted:
                print(f"    schedule: {outcome.reason}", flush=True)
            return rows

    for example in step3_examples():
        outcome = schedule_placement(
            example.placed.cnf,
            auxiliary=example.auxiliary,
            variable_schedules=example.placed.variable_schedules,
            clause_schedules=example.placed.clause_schedules,
        )
        ok = outcome.accepted and not outcome.balance.admitted and not outcome.built_slice
        add(f"supplied-{example.name}", example.placed.cnf.n, example.placed.cnf.m, brute_count(example.placed.cnf), ok, "OK" if ok else "FAIL")
        print(
            f"  [I] supplied {example.name}  "
            f"{'OK' if ok else 'FAIL'}  balance admitted={outcome.balance.admitted}",
            flush=True,
        )
        if not ok:
            print(f"    {outcome.reason}", flush=True)
            return rows

    small = schedule_placement(_units(2), r=2, Q=1, R=1)
    small_ok = (
        small.placed is None
        and not small.accepted
        and not small.balance.admitted
        and not small.built_slice
        and "No placement was invented" in small.reason
    )
    add("refuse-n2", 2, 2, -1, small_ok, "REFUSED" if small_ok else "FAIL")
    print(f"  [I] refuse n=2  {'REFUSED' if small_ok else 'FAIL'}", flush=True)
    if not small_ok:
        print(f"    {small.reason}", flush=True)
        return rows

    admitted = schedule_placement(_units(64), r=2, Q=1, R=1)
    admitted_ok = (
        admitted.balance.admitted
        and admitted.placed is None
        and not admitted.built_slice
        and admitted.balancer_attempted
        and not admitted.balancer_built
        and admitted.multiplicity == 16
        and admitted.multiplicity_ok is True
        and admitted.outside_group_size == 32
        and "2**32" in admitted.reason
        and "place_formula was not called" in admitted.reason
    )
    add("admit-n64-budget", 64, 64, -1, admitted_ok, "REFUSED" if admitted_ok else "FAIL")
    print(
        f"  [I] n=64 meets the inequalities; design budget refuses the placement  "
        f"{'REFUSED' if admitted_ok else 'FAIL'}",
        flush=True,
    )
    if not admitted_ok:
        print(f"    {admitted.reason}", flush=True)
        return rows

    try:
        build_auxiliary_slice(2, FULL_M0_MIN_Q, FULL_M0_MIN_R)
    except ValueError as exc:
        full_message = str(exc)
    else:
        full_message = ""
    full = schedule_placement(_units(1), r=2, Q=FULL_M0_MIN_Q, R=FULL_M0_MIN_R)
    full_ok = (
        full.placed is None
        and full.allocation_refused
        and not full.built_slice
        and "not allocated" in full.reason
        and "2**384" in full_message
        and "not built" in full_message
    )
    add("refuse-full-m0", 1, 1, -1, full_ok, "REFUSED" if full_ok else "FAIL")
    print(
        f"  [I] full M0 refused before allocation  {'REFUSED' if full_ok else 'FAIL'}",
        flush=True,
    )
    if not full_ok:
        print(f"    build: {full_message}", flush=True)
        print(f"    schedule: {full.reason}", flush=True)
        return rows
    return rows


def run_step8_suite() -> list[dict]:
    """Successive-choice balancer. Refusals are not timed.

    Every row here either refuses before a matrix is allocated, or is the
    direct checker call on a formula the inequalities do not admit. That
    direct call is timed because the design check did run. The clock is one
    machine's measurement, not a complexity result.
    """
    print("=== (J) Successive-choice balancer ===")
    print(
        "  The search runs only when 2(b+w)D^2 <= M and 2Qr^2 <= n hold and "
        f"the slice is under the cap {SLICE_MAX_DIMENSION}."
    )
    print(
        "  It tries symbols in {0, ..., r-1}^Q. It does not enumerate the "
        "schedule cube and it does not sample a placement."
    )
    print(
        f"  design_problems_at is called only when its Boolean walk is at most "
        f"{BALANCER_DESIGN_BOOLEAN_BUDGET} points. Above that, the call refuses."
    )
    print("  Theorem 1.1 is not implemented. This is not a 2^{o(n)} algorithm.")
    rows: list[dict] = []

    def add(name: str, n: int, m: int, brute: int, ok: bool, label: str) -> None:
        rows.append(
            {
                "name": name,
                "n": n,
                "m": m,
                "brute": brute,
                "ok": ok,
                "label": label,
            }
        )

    small = schedule_placement(_units(2), r=2, Q=1, R=1)
    small_ok = (
        small.placed is None
        and not small.accepted
        and not small.balance.admitted
        and not small.balancer_attempted
        and "No placement was invented" in small.reason
        and "was not run" in small.reason
    )
    add("refuse-n2", 2, 2, -1, small_ok, "REFUSED" if small_ok else "FAIL")
    print(f"  [J] refuse n=2  {'REFUSED' if small_ok else 'FAIL'}", flush=True)
    if not small_ok:
        print(f"    {small.reason}", flush=True)
        return rows

    toy = schedule_placement(_units(64), r=2, Q=1, R=1)
    toy_ok = (
        toy.balance.admitted
        and toy.placed is None
        and not toy.built_slice
        and not toy.accepted
        and toy.balancer_attempted
        and not toy.balancer_built
        and toy.multiplicity == 16
        and toy.multiplicity_ok is True
        and toy.group_size_ok is True
        and toy.outside_group_size == 32
        and toy.trials == 256
        and "2**32" in toy.reason
        and str(BALANCER_DESIGN_BOOLEAN_BUDGET) in toy.reason
        and "place_formula was not called" in toy.reason
    )
    add("admit-n64-budget", 64, 64, -1, toy_ok, "REFUSED" if toy_ok else "FAIL")
    print(
        f"  [J] n=64 admitted; multiplicity 16; design budget refuses  "
        f"{'REFUSED' if toy_ok else 'FAIL'}",
        flush=True,
    )
    if not toy_ok:
        print(f"    {toy.reason}", flush=True)
        return rows

    wide = schedule_placement(_units(1024), r=2, Q=2, R=1)
    wide_ok = (
        wide.balance.admitted
        and wide.balance.degree_gap_ok is True
        and wide.placed is None
        and not wide.built_slice
        and wide.multiplicity == 64
        and wide.multiplicity_ok is True
        and wide.outside_group_size == 512
        and "2**512" in wide.reason
        and "place_formula was not called" in wide.reason
    )
    add("admit-n1024-budget", 1024, 1024, -1, wide_ok, "REFUSED" if wide_ok else "FAIL")
    print(
        f"  [J] n=1024, R<Q, admitted; design budget refuses  "
        f"{'REFUSED' if wide_ok else 'FAIL'}",
        flush=True,
    )
    if not wide_ok:
        print(f"    {wide.reason}", flush=True)
        return rows

    over = schedule_placement(_units(1), r=2, Q=5, R=1)
    over_ok = (
        over.balance.D == 1024
        and over.allocation_refused
        and over.placed is None
        and not over.balancer_attempted
        and "not allocated" in over.reason
    )
    add("refuse-d1024", 1, 1, -1, over_ok, "REFUSED" if over_ok else "FAIL")
    print(f"  [J] side length 1024 refused before allocation  {'REFUSED' if over_ok else 'FAIL'}", flush=True)
    if not over_ok:
        print(f"    {over.reason}", flush=True)
        return rows

    try:
        build_auxiliary_slice(2, FULL_M0_MIN_Q, FULL_M0_MIN_R)
    except ValueError as exc:
        full_message = str(exc)
    else:
        full_message = ""
    full = schedule_placement(_units(1), r=2, Q=FULL_M0_MIN_Q, R=FULL_M0_MIN_R)
    full_ok = (
        full.placed is None
        and full.allocation_refused
        and not full.balancer_attempted
        and "not allocated" in full.reason
        and "2**384" in full_message
    )
    add("refuse-full-m0", 1, 1, -1, full_ok, "REFUSED" if full_ok else "FAIL")
    print(
        f"  [J] full M0 refused before allocation  {'REFUSED' if full_ok else 'FAIL'}",
        flush=True,
    )
    if not full_ok:
        print(f"    build: {full_message}", flush=True)
        print(f"    schedule: {full.reason}", flush=True)
        return rows

    wiring_cnf = _units(8)
    wiring_check = check_balance(wiring_cnf.n, scopes_from_cnf(wiring_cnf), 2, 2, 1)
    wiring_search = _search_schedules(wiring_cnf, wiring_check)
    started = time.perf_counter()
    wiring = _finish_chosen_placement(wiring_cnf, wiring_check, wiring_search, None)
    elapsed = time.perf_counter() - started
    gate = schedule_placement(wiring_cnf, r=2, Q=2, R=1)
    wiring_ok = (
        not wiring_check.admitted
        and wiring.accepted
        and wiring.balancer_built
        and not wiring.caller_supplied
        and wiring.design_problems == ()
        and wiring.restoration_problems == ()
        and wiring.multiplicity == 1
        and gate.placed is None
        and not gate.balancer_attempted
        and "were not the gate" in wiring.reason
    )
    add("wiring-n8-not-admitted", 8, 8, -1, wiring_ok, "OK" if wiring_ok else "FAIL")
    print(
        f"  [J] n=8 is not admitted; a direct check passed design and T=0  "
        f"{'OK' if wiring_ok else 'FAIL'}  ({elapsed:.3f}s)",
        flush=True,
    )
    print(
        "  [J] schedule_placement on that formula refuses before the search",
        flush=True,
    )
    if not wiring_ok:
        print(f"    direct: {wiring.reason}", flush=True)
        print(f"    gate: {gate.reason}", flush=True)
        return rows
    return rows


def print_summary_table(
    interp_rows: List[dict],
    indep_rows: List[dict],
    bounded_rows: List[dict],
    blocked_rows: List[dict] | None = None,
    matrix_rows: List[dict] | None = None,
    hermite_rows: List[dict] | None = None,
    sparse_rows: List[dict] | None = None,
    scan_rows: List[dict] | None = None,
    step7_rows: List[dict] | None = None,
    step8_rows: List[dict] | None = None,
) -> None:
    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    numeric = bool(
        interp_rows
        or indep_rows
        or bounded_rows
        or blocked_rows
        or matrix_rows
        or hermite_rows
        or sparse_rows
        or step7_rows
        or step8_rows
    )
    if numeric:
        print(f"{'kind':<8} {'n':>4} {'m':>4} {'#F':>8} {'match':>6}")
        print("-" * 72)
    for r in interp_rows[:15]:
        print(
            f"{'interp':<8} {r['n']:4d} {r['m']:4d} {r['brute']:8d} "
            f"{'OK' if r['ok'] else 'FAIL':>6}"
        )
    if len(interp_rows) > 15:
        print(f"  ... ({len(interp_rows) - 15} more interp rows)")
    for r in indep_rows[:10]:
        print(
            f"{'indep':<8} {r['n']:4d} {r['m']:4d} {r['count']:8d} "
            f"{'OK' if r['ok'] else 'FAIL':>6}"
        )
    if len(indep_rows) > 10:
        print(f"  ... ({len(indep_rows) - 10} more indep rows)")
    for r in bounded_rows:
        print(
            f"{'bounded':<8} {r['n']:4d} {r['m']:4d} {r['brute']:8d} "
            f"{'OK' if r['ok'] else 'FAIL':>6}"
        )
    for r in blocked_rows or []:
        print(
            f"{'blocked':<8} {r['n']:4d} {r['m']:4d} {r['brute']:8d} "
            f"{'OK' if r['ok'] else 'FAIL':>6}"
        )
    for r in matrix_rows or []:
        print(
            f"{'matrix':<8} {r['n']:4d} {r['m']:4d} {r['brute']:8d} "
            f"{'OK' if r['ok'] else 'FAIL':>6}"
        )
    for r in hermite_rows or []:
        print(
            f"{'hermite':<8} {r['n']:4d} {r['m']:4d} {r['brute']:8d} "
            f"{'OK' if r['ok'] else 'FAIL':>6}"
        )
    for r in sparse_rows or []:
        brute = r["brute"]
        brute_cell = f"{brute:8d}" if isinstance(brute, int) and brute >= 0 else f"{'—':>8}"
        print(
            f"{'sparse':<8} {r['n']:4d} {r['m']:4d} {brute_cell} "
            f"{r['label']:>6}"
        )
    for r in step7_rows or []:
        brute = r["brute"]
        brute_cell = f"{brute:8d}" if isinstance(brute, int) and brute >= 0 else f"{'—':>8}"
        print(
            f"{'step7':<8} {r['n']:4d} {r['m']:4d} {brute_cell} "
            f"{r['label']:>6}"
        )
    for r in step8_rows or []:
        brute = r["brute"]
        brute_cell = f"{brute:8d}" if isinstance(brute, int) and brute >= 0 else f"{'—':>8}"
        print(
            f"{'step8':<8} {r['n']:4d} {r['m']:4d} {brute_cell} "
            f"{r['label']:>6}"
        )
    if numeric:
        print("-" * 72)
    if interp_rows:
        print(
            f"Interpolation (B): {sum(1 for r in interp_rows if r['ok'])}/{len(interp_rows)} passed"
        )
    if indep_rows:
        print(
            f"Independent   (A): {sum(1 for r in indep_rows if r['ok'])}/{len(indep_rows)} passed"
        )
    if bounded_rows:
        print(
            f"Bounded occ.  (C): {sum(1 for r in bounded_rows if r['ok'])}/{len(bounded_rows)} passed"
        )
    if blocked_rows:
        print(
            f"Blocked sum   (D): {sum(1 for r in blocked_rows if r['ok'])}/{len(blocked_rows)} passed"
        )
    if matrix_rows:
        print(
            f"Matrix slice  (E): {sum(1 for r in matrix_rows if r['ok'])}/{len(matrix_rows)} passed"
        )
    if hermite_rows:
        print(
            f"Hermite       (F): {sum(1 for r in hermite_rows if r['ok'])}/{len(hermite_rows)} passed"
        )
    if sparse_rows:
        recorded = sum(1 for r in sparse_rows if r["label"] == "RECORDED")
        print(
            f"Sparsify      (G): {sum(1 for r in sparse_rows if r['ok'])}/{len(sparse_rows)} passed"
            + (f" ({recorded} recorded past the brute wall)" if recorded else "")
        )
    if scan_rows:
        agreed = sum(1 for row in scan_rows if row["status"] == "agreed")
        refused = sum(1 for row in scan_rows if row["status"] == "refused")
        print(
            f"Scan          (H): {sum(1 for row in scan_rows if row['ok'])}/{len(scan_rows)} passed"
            f" ({agreed} agreed, {refused} refused as recorded)"
        )
        for row in scan_rows:
            print(f"  {row['label']:<8} {row['name']}")
    if step7_rows:
        refused = sum(1 for row in step7_rows if row["label"] == "REFUSED")
        print(
            f"Balance/slice (I): {sum(1 for row in step7_rows if row['ok'])}/{len(step7_rows)} passed"
            + (f" ({refused} refused as recorded)" if refused else "")
        )
    if step8_rows:
        refused = sum(1 for row in step8_rows if row["label"] == "REFUSED")
        print(
            f"Balancer      (J): {sum(1 for row in step8_rows if row['ok'])}/{len(step8_rows)} passed"
            + (f" ({refused} refused as recorded)" if refused else "")
        )
    print("=" * 72)


def _selected_scope(steps: Sequence[int]) -> str:
    parts = ["Scope: the requested steps only."]
    if 2 in steps:
        parts.append(
            "Blocked sums match full-space H(t) on the hand-built families. "
            "Popcount measurement tuples grow as (r+1)^b."
        )
    if 3 in steps:
        parts.append(
            f"Side length is at most {SLICE_MAX_DIMENSION}. "
            "The parameter list with R≥32 and 6R≤Q is not built."
        )
    if 4 in steps:
        parts.append(
            "Hermite recovery matches brute force at H(0). The condition count "
            "stays D+1. Each Hermite node still enumerates 2^n assignments."
        )
    if 5 in steps:
        parts.append(
            "Where Step 5 ran, #F' = #F inside the brute-force window, and rows "
            "past that window are recorded rather than counted. Step 3 is not applied."
        )
    if 6 in steps:
        parts.append(
            "The parameter scan records agreement with brute force where a "
            "construction ran, and a refusal where the cube, the tuple cap, or an "
            "unfinished construction stopped the cell. Refused pipelines were not "
            "timed. The scan does not close the distance to a full counting algorithm."
        )
    if 7 in steps:
        parts.append(
            f"Step 7 checks slices up to side length {SLICE_MAX_DIMENSION}. "
            "Balance inequalities that fail are refused without choosing a schedule. "
            "Caller-supplied schedules, including the hand-written Step 3 placements, "
            "are checked with design_problems_at. Full M0 was refused before allocation. "
            "The admitted point is handed to the Step 8 search."
        )
    if 8 in steps:
        parts.append(
            "Step 8 runs successive choice only when the balance inequalities "
            "hold and the slice is under the cap. On the admitted rows the "
            "multiplicity and group bounds hold, and design_problems_at is not "
            "called because its Boolean walk is above the budget. Formulas that "
            "fail the inequalities still refuse with no placement. Full M0 and "
            "side length 1024 were refused before allocation. A direct check on "
            "eight units, which are not admitted, passed design and T=0; "
            "schedule_placement does not make that call. The schedule cube was "
            "not enumerated."
        )
    parts.append("No Theorem 1.1, no 2^{o(n)} claim.")
    return " ".join(parts)


def _run_selected(
    steps: Sequence[int],
    growth_groups: Sequence[int],
    scan_profile: str,
) -> int:
    """Run a subset of Steps 2–8. The full path stays in ``main``."""
    label = ", ".join(str(number) for number in steps)
    noun = "Step" if len(steps) == 1 else "Steps"
    print(f"counting_sat_algebra.verify — {noun} {label} only")
    print("Theorem 1.1 is not implemented. This program does not refute ETH.")
    print()

    example_rows: list[dict] = []
    matrix_rows: list[dict] = []
    hermite_rows: list[dict] = []
    sparse_rows: list[dict] = []
    scan_rows: list[dict] = []
    step7_rows: list[dict] = []
    step8_rows: list[dict] = []

    if 2 in steps:
        example_rows, growth_rows = run_step2_suite(growth_groups)
        if (
            not example_rows
            or not all(row["ok"] for row in example_rows)
            or len(growth_rows) != len(growth_groups)
            or not all(row["ok"] for row in growth_rows)
        ):
            print("FAIL: blocked sum disagreed with the full sum", file=sys.stderr)
            return 1
        print()
    if 3 in steps:
        matrix_rows = run_step3_suite()
        print()
        if not matrix_rows or not all(row["ok"] for row in matrix_rows):
            print("FAIL: replacement-matrix slice disagreed with its checks", file=sys.stderr)
            return 1
        print()
    if 4 in steps:
        hermite_rows = run_step4_suite()
        print()
        if not hermite_rows or not all(row["ok"] for row in hermite_rows):
            print("FAIL: Hermite recovery disagreed with brute force", file=sys.stderr)
            return 1
        print()
    if 5 in steps:
        sparse_rows = run_step5_suite()
        print()
        if not sparse_rows or not all(row["ok"] for row in sparse_rows):
            print("FAIL: occurrence copying disagreed with its checks", file=sys.stderr)
            return 1
        print()
    if 6 in steps:
        scan_rows = run_step6_suite(scan_profile)
        print()
        if not scan_rows or not all(row["ok"] for row in scan_rows):
            print("FAIL: parameter scan disagreed with its checks", file=sys.stderr)
            return 1
        print()
    if 7 in steps:
        step7_rows = run_step7_suite(heavy=scan_profile == "full")
        print()
        if not step7_rows or not all(row["ok"] for row in step7_rows):
            print("FAIL: Step 7 slice or balance check disagreed", file=sys.stderr)
            return 1
        print()
    if 8 in steps:
        step8_rows = run_step8_suite()
        print()
        if not step8_rows or not all(row["ok"] for row in step8_rows):
            print("FAIL: Step 8 balancer disagreed with its checks", file=sys.stderr)
            return 1
        print()

    print_summary_table(
        [],
        [],
        [],
        example_rows or None,
        matrix_rows or None,
        hermite_rows or None,
        sparse_rows or None,
        scan_rows or None,
        step7_rows or None,
        step8_rows or None,
    )
    print()
    print("ALL TESTS PASSED")
    print(_selected_scope(steps))
    return 0


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Counting-polynomial skeleton verifier (Steps 0 through 8)"
    )
    parser.add_argument("--trials-b", type=int, default=50, help="interpolation trials")
    parser.add_argument("--trials-a", type=int, default=20, help="independent-component trials")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--R", type=int, default=1, help="replacement degree for suite B")
    parser.add_argument("--quick", action="store_true", help="fewer and smaller tests")
    parser.add_argument(
        "--step2",
        action="store_true",
        help="run only the Step 2 blocked-sum checks",
    )
    parser.add_argument(
        "--step3",
        action="store_true",
        help="run only the Step 3 tiny replacement-matrix checks",
    )
    parser.add_argument(
        "--step4",
        action="store_true",
        help="run the Step 4 Hermite / Hasse checks",
    )
    parser.add_argument(
        "--step5",
        action="store_true",
        help="run the Step 5 occurrence-copying checks",
    )
    parser.add_argument(
        "--step6",
        action="store_true",
        help="run the Step 6 parameter scan (full grid; --quick keeps the small grid)",
    )
    parser.add_argument(
        "--step7",
        action="store_true",
        help=(
            "run the Step 7 slice and balance checks "
            "(--quick skips the side-length-256 point)"
        ),
    )
    parser.add_argument(
        "--step8",
        action="store_true",
        help=(
            "run the Step 8 successive-choice balancer "
            "(refusals are not timed; the design walk above the budget is not started)"
        ),
    )
    args = parser.parse_args(argv)

    n_min, n_max = 4, 10
    bounded_specs: List[Tuple[int, int, int]] = [(6, 4, 1), (6, 6, 1), (7, 5, 1)]
    # Block size 2 popcount scan. Tuple count is 3^b, so b=5 is still a 243-tuple sum.
    growth_groups: tuple[int, ...] = (1, 2, 3, 4, 5)
    if args.quick:
        args.trials_b = min(args.trials_b, 6)
        args.trials_a = min(args.trials_a, 6)
        n_max = 7
        bounded_specs = [(6, 4, 1), (6, 6, 1)]
        growth_groups = (1, 2, 3, 4)

    selected = [
        number
        for number, enabled in (
            (2, args.step2),
            (3, args.step3),
            (4, args.step4),
            (5, args.step5),
            (6, args.step6),
            (7, args.step7),
            (8, args.step8),
        )
        if enabled
    ]
    if selected:
        return _run_selected(
            selected,
            growth_groups,
            "quick" if args.quick else "full",
        )

    print("counting_sat_algebra.verify — Steps 0, 1, 2, 3, 4, 5, 6, 7, and 8")
    print("Counting polynomial, univariate interpolation, bounded-occurrence checks,")
    print("a blocked sum on hand-built replacements, a tiny replacement matrix,")
    print("Hermite recovery of H(0) from Hasse coefficients,")
    print("occurrence copying with #F' = #F inside the brute-force window,")
    print("a parameter scan of those pieces against brute force,")
    print("a balance-check scaffold for placements into a larger slice,")
    print("and a machine-limited successive-choice balancer.")
    print("Theorem 1.1 is not implemented. This program does not refute ETH.")
    print()

    print("=== (A) Independent components ===")
    indep_rows = run_independent_suite(trials=args.trials_a, seed=args.seed + 1)
    for r in indep_rows:
        status = "OK" if r["ok"] else "FAIL"
        print(
            f"  [A {r['trial']+1:02d}/{args.trials_a}] "
            f"n={r['n']} (={r['n1']}+{r['n2']}) m={r['m']}  "
            f"#F={r['count']}  product={r['product']}  {status}"
        )
        if not r["ok"]:
            print("FAIL: independent-component identity broken", file=sys.stderr)
            return 1
    print(f"  (A) all {len(indep_rows)} trials passed")
    print()

    print("=== (B) Univariate H(T) interpolation ===")
    print(f"  R={args.R}, trials={args.trials_b}, seed={args.seed}, n={n_min}..{n_max}")
    interp_rows = run_interpolation_suite(
        trials=args.trials_b,
        seed=args.seed,
        R=args.R,
        n_min=n_min,
        n_max=n_max,
    )
    if not all(r["ok"] for r in interp_rows) or len(interp_rows) < args.trials_b:
        print("FAIL: interpolation mismatch", file=sys.stderr)
        return 1
    print(f"  (B) all {len(interp_rows)} trials passed")
    print()

    print("=== (C) Bounded occurrence, degree caps m*k*R and 3*n*R ===")
    bounded_rows = run_bounded_suite(bounded_specs, seed=args.seed + 7)
    if not bounded_rows or not all(r["ok"] for r in bounded_rows):
        print("FAIL: bounded-occurrence degree check", file=sys.stderr)
        return 1
    print(f"  (C) all {len(bounded_rows)} instances passed")
    print()

    example_rows, growth_rows = run_step2_suite(growth_groups)
    if (
        not example_rows
        or not all(row["ok"] for row in example_rows)
        or len(growth_rows) != len(growth_groups)
        or not all(row["ok"] for row in growth_rows)
    ):
        print("FAIL: blocked sum disagreed with the full sum", file=sys.stderr)
        return 1
    print(f"  (D) all {len(example_rows)} hand-built instances passed")
    print()

    matrix_rows = run_step3_suite()
    if not matrix_rows or not all(row["ok"] for row in matrix_rows):
        print("FAIL: replacement-matrix slice disagreed with its checks", file=sys.stderr)
        return 1
    print(f"  (E) all {len(matrix_rows)} replacement-matrix slices passed")
    print()

    hermite_rows = run_step4_suite()
    if not hermite_rows or not all(row["ok"] for row in hermite_rows):
        print("FAIL: Hermite recovery disagreed with brute force", file=sys.stderr)
        return 1
    print(f"  (F) all {len(hermite_rows)} Hermite instances passed")
    print()

    sparse_rows = run_step5_suite()
    if not sparse_rows or not all(row["ok"] for row in sparse_rows):
        print("FAIL: occurrence copying disagreed with its checks", file=sys.stderr)
        return 1
    print(f"  (G) all {len(sparse_rows)} occurrence-copying instances passed")
    print()

    scan_rows = run_step6_suite("quick" if args.quick else "full")
    if not scan_rows or not all(row["ok"] for row in scan_rows):
        print("FAIL: parameter scan disagreed with its checks", file=sys.stderr)
        return 1
    print(f"  (H) all {len(scan_rows)} parameter-scan cells passed")
    print()

    step7_rows = run_step7_suite(heavy=not args.quick)
    if not step7_rows or not all(row["ok"] for row in step7_rows):
        print("FAIL: Step 7 slice or balance check disagreed", file=sys.stderr)
        return 1
    print(f"  (I) all {len(step7_rows)} Step 7 rows passed")
    print()

    step8_rows = run_step8_suite()
    if not step8_rows or not all(row["ok"] for row in step8_rows):
        print("FAIL: Step 8 balancer disagreed with its checks", file=sys.stderr)
        return 1
    print(f"  (J) all {len(step8_rows)} Step 8 rows passed")
    print()

    print("=== Optional: self-reduction demo ===")
    run_self_reduction_demo(seed=args.seed + 99)
    print()

    print_summary_table(
        interp_rows,
        indep_rows,
        bounded_rows,
        example_rows,
        matrix_rows,
        hermite_rows,
        sparse_rows,
        scan_rows,
        step7_rows,
        step8_rows,
    )
    print()
    print("ALL TESTS PASSED")
    print(
        "Scope: H(0)=#F by brute-force sampling plus Lagrange, the disjoint-product "
        "identity, degree-cap checks on occurrence-≤3 3-CNFs, blocked-sum agreement "
        "with full-space H(t) on hand-built replacements, a replacement "
        f"matrix of side length at most {SLICE_MAX_DIMENSION}, Hermite recovery "
        "of H(0) from Hasse coefficients at fewer distinct sample points than D+1, "
        "occurrence copying with #F' = #F inside the brute-force window, a "
        "parameter scan that records matches and refusals on a machine-limited grid, "
        "and a balance scaffold that refuses to choose schedules when either of "
        "2(b+w)D^2 <= M and 2Qr^2 <= n fails. "
        "The Hermite condition count stays D+1, and each of those points still "
        "enumerates 2^n assignments. Copying can push the reduced n past 20; those "
        "rows are recorded and not counted. Popcount measurement tuples grow as (r+1)^b. "
        "Refused scan cells were not timed. Caller-supplied schedules are checked. "
        "When the balance inequalities hold, successive choice aims at multiplicity "
        "at most 2M/D and group size at most 2n/r. On the admitted rows that walk "
        "is above the design budget, so no placement is accepted. The schedule "
        "cube was not enumerated. "
        f"The parameter list with R≥32 and 6R≤Q is not built (side length 2**{FULL_M0_LOG2_DIMENSION_AT_R2} "
        "at r=2). No Theorem 1.1, no 2^{o(n)} claim."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
