"""
Bounded-occurrence 3-CNF testbed (paper §2.1, experimental only).

A feasible exact 3-CNF has 3m literal occurrences and at most 3 per variable,
so 3m ≤ 3n, i.e. m ≤ n. The generator builds m consecutive width-3 windows
on a random relabeling of the variables. Each variable then sits in at most
three windows.

Degree caps checked here, for replacement degree R and clause width k = 3:

  D_tight = m * k * R = M * R,   M = sum_e |I_e|
  D_loose = 3 * n * R

On this family M = 3m ≤ 3n, so D_tight ≤ D_loose. Both are upper bounds on
deg H when each replacement has degree ≤ R. The experiment still evaluates
H by enumerating {0,1}^n. It checks the bound; it does not speed up the sum.
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import List, Sequence

from .cnf import CNF, brute_count
from .replacements import (
    degree_bound,
    eval_H_at,
    lagrange_eval,
    lagrange_eval_at_zero,
    random_replacements,
)

# Fixed-seed bench shipped under benches/. Order is part of the file format.
BENCH_SEED = 42
BENCH_SPECS: tuple[tuple[str, int, int, int], ...] = (
    ("n6-m4-r1", 6, 4, 1),
    ("n6-m6-r1", 6, 6, 1),
    ("n7-m5-r1", 7, 5, 1),
)


@dataclass(frozen=True)
class OccurrenceCheck:
    """Result of the bounded-occurrence checker."""

    ok: bool
    n: int
    m: int
    width_max: int
    max_occurrence: int
    total_literal_occurrences: int
    violations: tuple[str, ...]


def occurrence_counts(cnf: CNF) -> List[int]:
    """How many clauses contain each variable. One count per clause, not per literal."""
    counts = [0] * cnf.n
    for clause in cnf.clauses:
        for var in {v for v, _neg in clause}:
            counts[var] += 1
    return counts


def check_bounded_occurrence(
    cnf: CNF,
    max_occ: int = 3,
    width: int = 3,
) -> OccurrenceCheck:
    """
    Accept an exact width-`width` CNF whose variables each occur in at most
    `max_occ` clauses. Clauses must use `width` distinct variables.
    """
    violations: List[str] = []
    width_max = 0
    for index, clause in enumerate(cnf.clauses):
        width_max = max(width_max, len(clause))
        variables = [v for v, _neg in clause]
        if len(clause) != width:
            violations.append(f"clause {index} has width {len(clause)}, expected {width}")
        if len(set(variables)) != len(variables):
            violations.append(f"clause {index} repeats a variable")
    counts = occurrence_counts(cnf)
    max_seen = max(counts) if counts else 0
    for var, count in enumerate(counts):
        if count > max_occ:
            violations.append(f"variable {var} occurs in {count} clauses (cap {max_occ})")
    total = sum(len(clause) for clause in cnf.clauses)
    return OccurrenceCheck(
        ok=not violations,
        n=cnf.n,
        m=cnf.m,
        width_max=width_max,
        max_occurrence=max_seen,
        total_literal_occurrences=total,
        violations=tuple(violations),
    )


def degree_bound_3n(n: int, R: int) -> int:
    """Loose cap deg H ≤ 3 n R from M ≤ 3n and deg(replacement) ≤ R."""
    if n < 0 or R < 0:
        raise ValueError("n and R must be nonnegative")
    return 3 * n * R


def random_bounded_3cnf(
    n: int,
    m: int,
    rng: random.Random | None = None,
    max_occ: int = 3,
) -> CNF:
    """
    Exact 3-CNF on n variables with m clauses and at most 3 occurrences each.

    Feasible only when m ≤ n (because 3m ≤ 3n) and n ≥ 3. `max_occ` must be
    at least 3; the construction guarantees occurrence ≤ 3, which also meets
    any larger cap.
    """
    if n < 3:
        raise ValueError("need n >= 3 for a 3-CNF")
    if m < 0:
        raise ValueError("m must be nonnegative")
    if max_occ < 3:
        raise ValueError("this generator guarantees occurrence <= 3; max_occ must be >= 3")
    if m > n:
        raise ValueError(
            f"infeasible bounded 3-CNF: m={m} > n={n} "
            f"(3m literal slots cannot fit in {max_occ} occurrences each)"
        )
    rng = rng or random.Random()
    labels = list(range(n))
    rng.shuffle(labels)
    start = rng.randrange(n) if n else 0
    clauses: List[list[tuple[int, bool]]] = []
    for step in range(m):
        base = (start + step) % n
        variables = [labels[(base + offset) % n] for offset in range(3)]
        clause = [(var, bool(rng.getrandbits(1))) for var in variables]
        clauses.append(clause)
    cnf = CNF(n, clauses)
    report = check_bounded_occurrence(cnf, max_occ=3, width=3)
    if not report.ok:
        raise RuntimeError(f"generator invariant failed: {report.violations}")
    return cnf


def equispaced_degree(samples: Sequence[Fraction]) -> int:
    """
    Degree of the polynomial through samples at t = 1, 2, ..., len(samples).

    Uses exact forward differences. A constant polynomial, including the zero
    polynomial, is reported as degree 0 when there is at least one sample.
    """
    if not samples:
        raise ValueError("need at least one sample")
    layer = [Fraction(value) for value in samples]
    degree = 0
    nonzero = any(value != 0 for value in layer)
    if not nonzero:
        return 0
    order = 0
    while len(layer) > 1:
        layer = [layer[i + 1] - layer[i] for i in range(len(layer) - 1)]
        order += 1
        if any(value != 0 for value in layer):
            degree = order
            nonzero = True
    return degree


def _sample_H(cnf: CNF, family, t_max: int) -> List[Fraction]:
    """Return [H(1), H(2), ..., H(t_max)]."""
    if t_max < 1:
        raise ValueError("t_max must be positive")
    return [eval_H_at(cnf, family, Fraction(t)) for t in range(1, t_max + 1)]


def _recover_at_cap(
    samples: Sequence[Fraction],
    degree_cap: int,
) -> tuple[Fraction, bool, int]:
    """
    Interpolate deg ≤ degree_cap from H(1)..H(D+1).

    Returns (H(0), extra point H(D+2) matches the interpolant, observed degree
    on H(1)..H(D+2)).
    """
    if degree_cap < 0:
        raise ValueError("degree cap must be nonnegative")
    needed = degree_cap + 2
    if len(samples) < needed:
        raise ValueError(f"need {needed} samples, have {len(samples)}")
    xs = [Fraction(t) for t in range(1, degree_cap + 2)]
    ys = list(samples[: degree_cap + 1])
    h0 = lagrange_eval_at_zero(xs, ys)
    predicted = lagrange_eval(xs, ys, Fraction(degree_cap + 2))
    agrees = predicted == samples[degree_cap + 1]
    observed = equispaced_degree(samples[:needed])
    return h0, agrees, observed


def run_degree_bound_experiment(
    cnf: CNF,
    R: int = 1,
    rng: random.Random | None = None,
) -> dict:
    """
    Build degree-≤ R replacements and recover #F at both degree caps.

    `ok_mkR` / `ok_3nR` require H(0) = #F, agreement at one extra node, and
    observed degree within the cap. The observed degree is measured on the
    tighter cap's sample window.
    """
    if R < 0:
        raise ValueError("R must be nonnegative")
    rng = rng or random.Random()
    report = check_bounded_occurrence(cnf)
    family = random_replacements(cnf, R=R, rng=rng)
    if family.R != R:
        raise RuntimeError("replacement family degree drifted")
    brute = brute_count(cnf)
    tight = degree_bound(cnf, R=R, k=3)
    loose = degree_bound_3n(cnf.n, R)
    if report.ok and report.total_literal_occurrences != 3 * cnf.m:
        raise RuntimeError("exact 3-CNF must have M = 3m")
    samples = _sample_H(cnf, family, max(tight, loose) + 2)
    h0_tight, extra_tight, observed = _recover_at_cap(samples, tight)
    h0_loose, extra_loose, observed_loose = _recover_at_cap(samples, loose)
    # The tighter window already determines the degree when it is ≤ tight.
    # observed_loose is computed on a longer window and must agree.
    if observed != observed_loose:
        # Longer window can only confirm the same degree when deg ≤ tight.
        # If the loose window sees a higher degree, keep the larger figure.
        observed_report = max(observed, observed_loose)
    else:
        observed_report = observed

    def _matches(h0: Fraction, extra: bool, cap: int) -> bool:
        return (
            h0 == Fraction(brute)
            and h0.denominator == 1
            and extra
            and observed_report <= cap
        )

    mass = report.total_literal_occurrences
    loose_applies = mass <= 3 * cnf.n
    return {
        "n": cnf.n,
        "m": cnf.m,
        "R": R,
        "bounded_ok": report.ok,
        "max_occurrence": report.max_occurrence,
        "M": mass,
        "loose_cap_applies": loose_applies,
        "D_mkR": tight,
        "D_3nR": loose,
        "brute": brute,
        "H0_mkR": h0_tight,
        "H0_3nR": h0_loose,
        "extra_ok_mkR": extra_tight,
        "extra_ok_3nR": extra_loose,
        "observed_degree": observed_report,
        "ok_mkR": _matches(h0_tight, extra_tight, tight),
        "ok_3nR": loose_applies and _matches(h0_loose, extra_loose, loose),
    }


def cnf_to_record(name: str, cnf: CNF, R: int) -> dict:
    report = check_bounded_occurrence(cnf)
    return {
        "name": name,
        "n": cnf.n,
        "m": cnf.m,
        "R": R,
        "max_occurrence": report.max_occurrence,
        "M": report.total_literal_occurrences,
        "clauses": [[[var, neg] for var, neg in clause] for clause in cnf.clauses],
    }


def cnf_from_record(record: dict) -> CNF:
    clauses = [
        [(int(var), bool(neg)) for var, neg in clause] for clause in record["clauses"]
    ]
    return CNF(int(record["n"]), clauses)


def build_bench_document(seed: int = BENCH_SEED) -> dict:
    """Deterministic bench document. Same seed and spec list ⇒ same JSON object."""
    instances = []
    for index, (name, n, m, R) in enumerate(BENCH_SPECS):
        cnf = random_bounded_3cnf(n, m, rng=random.Random(seed + index))
        record = cnf_to_record(name, cnf, R)
        instances.append(record)
    return {
        "seed": seed,
        "description": (
            "Fixed-seed exact 3-CNFs with at most 3 occurrences of each variable. "
            "Regenerate with counting_sat_algebra.bounded_occurrence.build_bench_document."
        ),
        "instances": instances,
    }


def bench_path(root: Path | None = None) -> Path:
    base = root or Path(__file__).resolve().parents[2]
    return base / "benches" / "bounded_occurrence_seed42.json"


def write_bench_document(path: Path, seed: int = BENCH_SEED) -> dict:
    document = build_bench_document(seed)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return document


def load_bench_document(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))
