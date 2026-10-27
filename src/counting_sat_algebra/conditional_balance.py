"""Lemma 4.5: exact conditional expectation of Phi/(M(b+w)) + Xi/(Qn).

Collision probabilities factor into variable-label and factor-label equality.
The unfixed labels are uniform and independent. Aggregating those pair
probabilities avoids enumerating their completions and the schedule cube.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence


@dataclass(frozen=True)
class ConditionalSchedules:
    variables: tuple[tuple[int, ...], ...]
    factors: tuple[tuple[int, ...], ...]
    initial_potential: Fraction
    final_potential: Fraction
    expectation_path: tuple[Fraction, ...]
    trials: int


def choose_schedules(
    n: int, scopes: Sequence[Sequence[int]], r: int, Q: int
) -> ConditionalSchedules:
    """Choose lexicographic minimizers of the reference's exact expectation.

    This finite construction also works outside the balance hypotheses; its
    load bounds are guaranteed only when the caller checks those hypotheses.
    Empty scopes get the all-zero label (the reference uses one-based labels).
    """
    if n < 1 or r < 2 or Q < 1:
        raise ValueError("need n >= 1, r >= 2, Q >= 1")
    if any(len(set(scope)) != len(scope) or any(v < 0 or v >= n for v in scope) for scope in scopes):
        raise ValueError("scopes must be subsets of the declared variables")
    counts = [0] * n
    neighbors: list[dict[int, int]] = [{} for _ in range(n)]
    for scope in scopes:
        for v in scope:
            counts[v] += 1
            for u in scope:
                if u != v:
                    neighbors[v][u] = neighbors[v].get(u, 0) + 1
    M = sum(counts)
    b = max(counts)
    w = max((len(scope) for scope in scopes), default=0)
    if not M:
        raise ValueError("the normalized potential requires positive occurrence mass")
    alphabet = tuple(itertools.product(range(r), repeat=Q))
    K = len(alphabet)
    norm = M * (b + w)
    same_v = sum(c * c for c in counts) - M
    same_f = sum(len(scope) ** 2 for scope in scopes) - M
    separate = M * M - M - same_v - same_f
    phi = Fraction(M) + Fraction(same_v + same_f, K) + Fraction(separate - M * M, K * K)
    xi = Fraction(Q * n * (r - 1), r)
    initial = phi / norm + xi / (Q * n)
    path = [initial]

    variables: list[tuple[int, ...]] = []
    marginal = [[0] * r for _ in range(Q)]
    fixed_mass: dict[tuple[int, ...], int] = {}
    total_mass = 0
    for v in range(n):
        shared_mass: dict[tuple[int, ...], int] = {}
        for u, amount in neighbors[v].items():
            if u < v:
                sigma = variables[u]
                shared_mass[sigma] = shared_mass.get(sigma, 0) + amount
        total_shared = sum(shared_mass.values())

        def delta(sigma: tuple[int, ...]) -> Fraction:
            collision = (
                counts[v] * (Fraction(fixed_mass.get(sigma, 0)) - Fraction(total_mass, K))
                + (K - 1) * (Fraction(shared_mass.get(sigma, 0)) - Fraction(total_shared, K))
            ) * Fraction(2, K)
            groups = 2 * (Fraction(sum(marginal[t][sigma[t]] for t in range(Q))) - Fraction(Q * v, r))
            return collision / norm + groups / (Q * n)

        sigma = min(alphabet, key=delta)
        path.append(path[-1] + delta(sigma))
        variables.append(sigma)
        total_mass += counts[v]
        fixed_mass[sigma] = fixed_mass.get(sigma, 0) + counts[v]
        for t, symbol in enumerate(sigma):
            marginal[t][symbol] += 1

    labels: dict[tuple[tuple[int, ...], tuple[int, ...]], int] = {}
    prior_mass: dict[tuple[int, ...], int] = {}
    factors: list[tuple[int, ...]] = []
    for scope in scopes:
        masses: dict[tuple[int, ...], int] = {}
        for v in scope:
            sigma = variables[v]
            masses[sigma] = masses.get(sigma, 0) + 1
        baseline = Fraction(sum(amount * prior_mass.get(sigma, 0) for sigma, amount in masses.items()), K)

        def delta(tau: tuple[int, ...]) -> Fraction:
            matched = sum(amount * labels.get((sigma, tau), 0) for sigma, amount in masses.items())
            return 2 * (matched - baseline) / norm

        tau = min(alphabet, key=delta) if scope else alphabet[0]
        path.append(path[-1] + delta(tau))
        factors.append(tau)
        for sigma, amount in masses.items():
            labels[sigma, tau] = labels.get((sigma, tau), 0) + amount
            prior_mass[sigma] = prior_mass.get(sigma, 0) + amount

    # Independent direct evaluation of the final Phi and Xi.
    phi_final = Fraction(sum(c * c for c in labels.values())) - Fraction(M * M, K * K)
    xi_final = Fraction(sum(size * size for row in marginal for size in row)) - Fraction(Q * n * n, r)
    final = phi_final / norm + xi_final / (Q * n)
    if path[-1] != final or any(right > left for left, right in zip(path, path[1:])):
        raise ArithmeticError("conditional expectation failed its exact potential check")
    return ConditionalSchedules(
        tuple(variables), tuple(factors), initial, final, tuple(path), (n + len(scopes)) * K
    )
