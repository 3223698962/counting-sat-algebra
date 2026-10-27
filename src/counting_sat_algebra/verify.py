#!/usr/bin/env python3
"""
Verification CLI for the Step 0 skeleton.

  python -m counting_sat_algebra.verify
  python -m counting_sat_algebra.verify --quick
"""

from __future__ import annotations

import argparse
import random
import sys
import time
from fractions import Fraction
from typing import List

from .cnf import brute_count, random_3cnf
from .independent import run_independent_suite
from .replacements import interpolate_H0
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
    Compare brute-force #F with interpolated H(0).
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


def print_summary_table(interp_rows: List[dict], indep_rows: List[dict]) -> None:
    print()
    print("=" * 64)
    print("SUMMARY")
    print("=" * 64)
    print(f"{'kind':<8} {'n':>4} {'m':>4} {'#F':>8} {'match':>6}")
    print("-" * 64)
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
    print("-" * 64)
    n_ok_b = sum(1 for r in interp_rows if r["ok"])
    n_ok_a = sum(1 for r in indep_rows if r["ok"])
    print(f"Interpolation (B): {n_ok_b}/{len(interp_rows)} passed")
    print(f"Independent   (A): {n_ok_a}/{len(indep_rows)} passed")
    print("=" * 64)


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Counting-polynomial skeleton verifier (Step 0)"
    )
    parser.add_argument("--trials-b", type=int, default=50, help="interpolation trials")
    parser.add_argument("--trials-a", type=int, default=20, help="independent-component trials")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--R", type=int, default=1, help="replacement degree")
    parser.add_argument("--quick", action="store_true", help="fewer and smaller tests")
    args = parser.parse_args(argv)

    n_min, n_max = 4, 10
    if args.quick:
        args.trials_b = min(args.trials_b, 6)
        args.trials_a = min(args.trials_a, 6)
        n_max = 7

    print("counting_sat_algebra.verify — counting polynomial + interpolation skeleton")
    print("Step 0 only. Does not implement a subexponential algorithm and does not refute ETH.")
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

    print("=== Optional: self-reduction demo ===")
    run_self_reduction_demo(seed=args.seed + 99)
    print()

    print_summary_table(interp_rows, indep_rows)
    print()
    print("ALL TESTS PASSED")
    print(
        "Scope: H(0)=#F recovered by Lagrange from brute-force samples H(1),…,H(D+1), "
        "and the disjoint-variable product identity. No design, no blocked-sum, no M0, "
        "no 2^{o(n)} claim."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
