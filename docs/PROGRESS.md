# Progress

Steps 0, 1, and 2 are implemented and checked against brute-force counting, or against the full-space sum H(t), on small instances. This repository does not refute ETH. It does not implement the matrix M0, Hermite interpolation, sparsification, or Theorem 1.1. The Step 2 measurement loop is exponential in the number of groups.

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

## Step 2 — simplified (r, τ, D)-design and blocked sum

Shipped in `design.py` and `blocked_sum.py`.

What was implemented is a simplification, not a paper-style design and not a slice of M0. An (r, τ, D)-design here is:

- a partition of the n variables into blocks of size at most r,
- τ explicit rational linear measurements, each a weight vector in Q^n, shared by adding the blocks' contributions,
- a nonnegative integer D stored as a degree cap and not read by the sum.

Every nonempty clause must lie inside one block. The hand-built replacement is

⟨P_{e,v}(T), X⟩ = X_v + Σ_{j=1}^{R} T^j ⟨α_{e,v,j}, μ(X)⟩,

so the value at T = 0 is still X_v. With that support, fixing the local measurements (λ_1, …, λ_b) and setting μ = λ_1 + … + λ_b factors H(t) into a product of sums inside the blocks. Pytest checks that product against `eval_H_at` (the full cube) and, on the same points, against the measurement formula summed over the cube without blocking.

### Agreement

All values below are exact rationals. Blocked sum and full-space H(t) matched at every listed point. H(0) also matched the brute-force model count. For t ≠ 0, H(t) is a rational and need not be a nonnegative integer; the three-block row is negative at t = 1 and at t = 2.

| instance | n | blocks | r | τ | R | D | #F | tuples | local evals | cube 2^n | points |
|----------|---|--------|---|---|---|---|----|--------|-------------|----------|--------|
| one-block-n3 | 3 | 1 | 3 | 1 | 1 | 3 | 7 | 4 | 8 | 8 | 0, 1, 2 |
| two-block-3cnf-n6 | 6 | 2 | 3 | 1 | 1 | 9 | 42 | 16 | 64 | 64 | 0, 1, 2, 3/2 |
| three-block-n6 | 6 | 3 | 2 | 1 | 1 | 6 | 27 | 27 | 108 | 64 | 0, 1, 2 |
| two-block-tau2-R2-n6 | 6 | 2 | 3 | 2 | 2 | 18 | 42 | 36 | 96 | 64 | 0, 1, 2 |

Sample values (blocked = full):

| instance | H(0) | H(1) | H(2) | other |
|----------|------|------|------|-------|
| one-block-n3 | 7 | 43/2 | 162 | |
| two-block-3cnf-n6 | 42 | 41002342 | 17038929522 | H(3/2) = 5506822263/4 |
| three-block-n6 | 27 | -17379/4 | -854028 | |
| two-block-tau2-R2-n6 | 42 | 46112915/16 | 198388412940 | |

`one-block-n3` is one 3-clause with occurrence coefficients 1, 1/2, and -1 on the popcount measurement. `two-block-3cnf-n6` is three 3-clauses, two blocks of size 3, every α coefficient 1. `three-block-n6` is three width-2 clauses, α alternating 1 and -1/2. `two-block-tau2-R2-n6` uses the same 3-CNF as the two-block row, with τ = 2 (popcount, plus a flag on the first variable of each block) and degree R = 2, coefficients written out per occurrence.

Two expansions were also checked by hand, not only by comparing code paths:

- two positive unit clauses on two blocks, ⟨P, X⟩ = X_v + T(X_0 + X_1), gives H(T) = 1 + 6T + 6T^2 (checked at 0, 1/2, 1, 2);
- one negative unit clause, ⟨P, X⟩ = X(1 + T), gives H(T) = 1 - T.

Further pytest rows, same agreement: τ = 0 (one measurement tuple, H constant and equal to #F); a weight of 1/2; seeded random α (seed 2026, R = 2, two blocks); R = 0, which still enumerates every popcount tuple; an empty clause (both sums are 0); a formula with no clauses (both sums are 2^n). Illegal partitions, a clause that meets two blocks, a bad α table, and measurement spaces above the safety cap are rejected.

### How the cost grows with the number of groups

Popcount measurement (τ = 1, every weight 1), equal blocks of size r, b groups. This implementation walks every realizable local measurement, including when a coefficient is 0:

- measurement tuples = (r + 1)^b
- local clause-product evaluations = b · 2^r · (r + 1)^{b-1}
- full cube = 2^{r b}

For r = 2 the tuple count is 3^b and the local-evaluation count is b · 4 · 3^{b-1}. Pytest locks b = 1..5. One run of `python -m counting_sat_algebra.verify --step2` on the agent machine compared blocked H(1) with the full sum and recorded the enumeration time inside `blocked_sum` (one process, not a benchmark suite):

| b | n | H(0) = #F | H(1) | tuples | local evals | cube | enum seconds |
|---|---|-----------|------|--------|-------------|------|--------------|
| 1 | 2 | 3 | 11 | 3 | 4 | 4 | 0.0001 |
| 2 | 4 | 9 | 1065 | 9 | 24 | 16 | 0.0003 |
| 3 | 6 | 27 | 330507 | 27 | 108 | 64 | 0.0013 |
| 4 | 8 | 81 | 205056049 | 81 | 432 | 256 | 0.0051 |
| 5 | 10 | 243 | 211822967883 | 243 | 1620 | 1024 | 0.0195 |

H(1) for b = 1 was also expanded by hand (the four assignments sum to 11). From b = 1 to b = 5 the local-evaluation count grows by a factor of 405 (4 to 1620). On that same run the measured enumeration time rose with b, from about 0.1 ms to about 20 ms. These times are too small to be a complexity result. The reproducible cost is the operation count, and it is exponential in b.

On this r = 2 family the local-evaluation count is already larger than the cube for every b ≥ 2 (24 > 16, 108 > 64, and so on). The three-block example in the agreement table shows the same pattern: 108 local evaluations versus a cube of 64. A different parameter point goes the other way and was checked in pytest: three blocks of size 3 give 64 measurement tuples, 384 local evaluations, and a cube of 512, and the sums still match. Either way the tuple count (r + 1)^b is exponential in the number of groups. For block size 1 it equals 2^n. Nothing in this step runs in 2^{o(n)}.

On the width-2 popcount family, #F = 3^b happens to equal the tuple count 3^b. That is one forbidden local assignment per block of size 2, set next to popcounts {0, 1, 2}. It is not the blocked sum counting models. At t = 1 the two columns separate (11 versus 3, 1065 versus 9, and so on).

### Obstacle

A design that could feed a subexponential counting algorithm would also have to keep the measurement space small at large n and hand a low-degree matrix to Step 3. This step does neither. D is bookkeeping. Clauses that cross blocks are out of scope for the factoring identity used here. Theorem 1.1 is not implemented.

## Next: Step 3

One runnable slice of the low-degree replacement / matrix construction. If the full matrix M0 is too heavy, record the obstacle and the simplification. That step is not in this tree.

## Known gaps

- No M0, no Hermite interpolation, no sparsification, no Theorem 1.1.
- Step 2 matches the full sum on the instances above. Measurement enumeration for the popcount design is (r + 1)^b, exponential in the number of groups, and on the recorded r = 2 scan it does more local work than summing the cube.
- Interpolation cost on the Step 0/1 path is still Θ(2^n · D).
- The degree table is the bench set plus one pytest draw, all at R = 1.
- Self-reduction deletes fixed variables from clauses and then calls brute force on all 2^n bits, so the integer it sees is not the number of extensions. The CLI only needs a yes/no answer and checks the witness it returns.
