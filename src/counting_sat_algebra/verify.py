#!/usr/bin/env python3
"""
Verification CLI: Step 0 (product identity + H(T) interpolation), a small
Step 1 bounded-occurrence degree check, and Step 2 blocked-sum agreement.

Usage:
  python -m counting_sat_algebra.verify
  python -m counting_sat_algebra.verify --quick
  python -m counting_sat_algebra.verify --step2
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
from .independent import run_independent_suite
from .replacements import eval_H_at, interpolate_H0
from .self_reduce import find_assignment_by_self_reduction


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


def print_summary_table(
    interp_rows: List[dict],
    indep_rows: List[dict],
    bounded_rows: List[dict],
    blocked_rows: List[dict] | None = None,
) -> None:
    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
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
    print("=" * 72)


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Counting-polynomial skeleton verifier (Steps 0, 1, and 2)"
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

    if args.step2:
        print("counting_sat_algebra.verify — Step 2 only")
    else:
        print("counting_sat_algebra.verify — Steps 0, 1, and 2")
        print("Counting polynomial, univariate interpolation, bounded-occurrence checks,")
        print("and a blocked sum on hand-built replacements.")
    print("Theorem 1.1 is not implemented. This program does not refute ETH.")
    print()

    if args.step2:
        example_rows, growth_rows = run_step2_suite(growth_groups)
        print()
        if (
            not example_rows
            or not all(row["ok"] for row in example_rows)
            or len(growth_rows) != len(growth_groups)
            or not all(row["ok"] for row in growth_rows)
        ):
            print("FAIL: blocked sum disagreed with the full sum", file=sys.stderr)
            return 1
        print_summary_table([], [], [], example_rows)
        print()
        print("ALL TESTS PASSED")
        print(
            "Scope: Step 2 only. On these hand-built families the blocked sum equals "
            "full-space H(t), including H(0) = #F. For popcount blocks of size 2 the "
            "measurement tuple count is 3^b, exponential in the number of groups. "
            "No matrix M0, no Theorem 1.1, no 2^{o(n)} claim."
        )
        return 0

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

    print("=== Optional: self-reduction demo ===")
    run_self_reduction_demo(seed=args.seed + 99)
    print()

    print_summary_table(interp_rows, indep_rows, bounded_rows, example_rows)
    print()
    print("ALL TESTS PASSED")
    print(
        "Scope: H(0)=#F by brute-force sampling plus Lagrange, the disjoint-product "
        "identity, degree-cap checks on occurrence-≤3 3-CNFs, and blocked-sum "
        "agreement with full-space H(t) on hand-built replacements. "
        "Popcount measurement tuples grow as (r+1)^b. "
        "No matrix M0, no Theorem 1.1, no 2^{o(n)} claim."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
