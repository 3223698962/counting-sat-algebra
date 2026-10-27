# Progress

Steps 0 and 1 are implemented and checked against brute-force counting on small instances. This repository does not refute ETH. It does not implement an (r, τ, D)-design, a blocked sum, the matrix M0, sparsification, or Theorem 1.1.

## Step 0

Shipped in `cnf.py`, `replacements.py`, `independent.py`, `self_reduce.py`, and `verify.py`.

- Clause polynomial Φ_e and brute-force #F.
- Disjoint-variable product: #F = #F1 · #F2 when no clause crosses the two groups.
- Random replacements of degree ≤ R with value X_v at T = 0. Lagrange interpolation on T = 1, …, D+1 with D = m·k·R recovers H(0) = #F.
- Pytest covers a fixed 3-clause on n = 3, several random formulas including one with n = 6, an R = 2 formula on n = 4, and the product identity.
- Each sample H(t) still sums 2^n assignments.

## Step 1 — bounded occurrence

Shipped in `bounded_occurrence.py`.

The generator builds m consecutive width-3 windows on a random relabeling. That is feasible only for m ≤ n, and every variable occurs in at most three clauses. The checker rejects a larger width, a repeated variable inside a clause, or a fourth occurrence.

Degree-cap experiment, replacement degree R = 1, fixed seeds. Both caps are interpolated in exact rational arithmetic. "deg" is the degree read from forward differences of H(1), H(2), …. "ok" means H(0) equals the brute-force count at both caps, the next node matches the interpolant, and the measured degree is at most the cap.

Bench file `benches/bounded_occurrence_seed42.json` (formula seed 42, 43, 44; replacement seed 42):

| instance | n | m | occ | M | D = m·3·R | D = 3·n·R | #F | H(0) | measured deg |
|----------|---|---|-----|---|-----------|-----------|----|------|--------------|
| n6-m4-r1 | 6 | 4 | 3 | 12 | 12 | 18 | 44 | 44 | 12 |
| n6-m6-r1 | 6 | 6 | 3 | 18 | 18 | 18 | 35 | 35 | 18 |
| n7-m5-r1 | 7 | 5 | 3 | 15 | 15 | 21 | 58 | 58 | 15 |

An additional pytest draw (`random_bounded_3cnf(6, 5, seed 7)`, replacements seed 8) gave #F = 33, measured degree 15, caps 15 and 18, and the same agreement.

On these four draws the measured degree equalled m·3·R, so that cap was tight for the sampled replacement family, and 3·n·R was strictly looser except on the instance with m = n. Both caps recovered #F. This matches the §2.1–§2.3 bookkeeping M = Σ |I_e| ≤ 3n and deg H ≤ M R on these instances. It is one family per formula, not a proof for every replacement, and it is not a faster way to evaluate the sum.

## Next: Step 2

(r, τ, D)-design, then a blocked sum compared with the full sum H(t) on a hand-built replacement family. The check is agreement with the full sum on small n. A subexponential runtime is not the goal of that step.

## Known gaps

- No design, no blocked sum, no M0, no Hermite interpolation, no sparsification.
- Interpolation cost is still Θ(2^n · D).
- The degree table is the bench set plus one pytest draw, all at R = 1.
- Self-reduction deletes fixed variables from clauses and then calls brute force on all 2^n bits, so the integer it sees is not the number of extensions. The CLI only needs a yes/no answer and checks the witness it returns.
