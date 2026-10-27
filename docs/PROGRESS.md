# Progress

Steps 0, 1, 2, a small-parameter Step 3, an exact Step 4 Hermite check, a Step 5 occurrence reduction, and a Step 6 parameter scan are implemented and compared with brute-force counting, or with the full-space sum H(t), on small instances. This repository does not refute ETH. It does not implement Theorem 1.1. Step 5 copies variables so that each occurs in at most three clauses and checks #F' = #F where both cubes fit in the brute-force window; the reduced length n' can be larger than the source length, and past n' = 20 the equality is recorded as unchecked. The matrix M0 is built only at side length r**(2Q) ≤ 64. The Step 2 measurement loop is exponential in the number of groups. Step 4 uses fewer distinct sample points than D+1 and still enumerates the cube at each point. Step 6 runs those pieces on one fixed grid and records a match or a refusal. A refusal is not timed. The scan does not close the distance listed at the end of this file.

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

## Step 3 — small-parameter replacement matrix

Shipped in `replacement_matrix.py`. This is the auxiliary matrix

M0 = V0 f_I(W0)

and the rows built from it, at side lengths 4, 9, 16, and 64. It is not the matrix at the parameter list R ≥ 32, 6R ≤ Q. Schedules are written out by hand. Hermite interpolation is not part of this matrix. Theorem 1.1 is not implemented.

Coordinates are pairs of length-Q schedules with entries in {0, …, r−1}. There are D = r**(2Q) of them. V0 exchanges the two schedules, W0 = (1/2) Σ_s (I − Π_s), and f_I is the polynomial of degree at most R−1 with f_I(j) = (−1)^j on R consecutive integers. The row is

⟨P_u(T), ζ⟩ = ((I − M0) ζ)_u + Σ_y M0[u, y] φ_{u,y}(T) ζ_y,

φ_{u,y}(T) = Π_{s ∈ E(u,y)} (1 − T/s).

At T = 0 every φ equals 1, so the row restores ζ_u. A nonzero M0 entry has at most R−1 mismatches, so the degree is at most R−1. Measurement rows are the product basis vectors whose minus-count lies outside the interval. Their entries are in {−1, 0, 1}. Pulling a row back along an explicit placement gives a replacement of the original variables. Unoccupied coordinates contribute 0.

The point r = 2, Q = 2, R = 2 has R = Q. That sits outside the hypothesis R < Q used for a runtime bound. The identities at that point are exact and were checked anyway.

### What was checked

On every slice below, the auxiliary checks passed: f_I(j) = (−1)^j on the interval, T = 0 restores each coordinate, degree ≤ R−1, (R−1)! M0 is an integer matrix, (R−1)! Q! clears the row denominators, a nonzero M0 entry has at most R−1 mismatches, M0 fixes each basis vector inside the interval, and each row of I−M0 lies in the span of the measurement rows. The measurement-row count equals the binomial count of basis vectors outside the interval.

On each placed formula, the coefficient row of occurrence (e, v) is e_v at T = 0, H(0) equals the brute-force count both by direct summation and by Lagrange at both caps, and the measured degree is at most M·(R−1) and at most m·k·R. At every t ∈ {1, …, Q} the schedules partition the variables, each group has size s with s·r ≤ 2n, each pulled-back measurement row is integral with absolute sum at most D, and every coefficient outside the clause's assigned group lies in the measurement span (over Q, and on Boolean points).

For r = 2 the size bound is s ≤ n, so every partition meets it. For `cross-r3`, r = 3 and n = 3, the bound is s ≤ 2, and each group has size 1.

| instance | r | Q | R | D | C rows | n | m | #F | deg H | M(R−1) | m·k·R | t checked |
|----------|---|---|---|---|--------|---|---|----|-------|--------|-------|-----------|
| swap-4 | 2 | 1 | 1 | 4 | 1 | 2 | 2 | 1 | 0 | 0 | 2 | 1 |
| remark-4-4 | 2 | 2 | 2 | 16 | 1 | 2 | 2 | 1 | 2 | 2 | 4 | 1, 2 |
| cross-r3 | 3 | 1 | 1 | 9 | 3 | 3 | 1 | 7 | 0 | 0 | 3 | 1 |
| r2-q3 | 2 | 3 | 2 | 64 | 10 | 3 | 1 | 7 | 2 | 3 | 6 | 1, 2, 3 |

`remark-4-4` meets the occurrence cap: the measured degree is 2, and M·(R−1) = 2. `r2-q3` has measured degree 2 under a cap of 3. The two degree-0 rows meet a cap of 0, so H is constant.

Sample values, all equal to the full-cube sum:

| instance | H(0) | H(1) | H(2) | other |
|----------|------|------|------|-------|
| swap-4 | 1 | 1 | 1 | |
| remark-4-4 | 1 | -3/4 | -2 | H(1/2) = 1/16 |
| cross-r3 | 7 | 7 | 7 | |
| r2-q3 | 7 | 8 | 29/3 | |

`swap-4` and `cross-r3` and `r2-q3` all have H(0) equal to a nonnegative integer model count. `remark-4-4` is negative at t = 1 and at t = 2. That is a rational value of H(t) for t ≠ 0, not a count.

### Hand checks

`swap-4`. Order the four coordinates (a, d) as (0,0), (0,1), (1,0), (1,1). f_I is the constant 1 and M0 is the exchange

```
1 0 0 0
0 0 1 0
0 1 0 0
0 0 0 1
```

The single measurement row is (0, 1, −1, 0). Each auxiliary row is the corresponding basis vector, so the degree is 0. The two positive unit clauses put variable 0 in group 0 and variable 1 in group 1, and assign each clause to the other group. The pulled-back measurement is X_0 − X_1. H(T) = 1, which is the model count of two positive units.

`remark-4-4`. f_I(z) = 1 − 2z on the interval {0, 1}. The M0 row of coordinate (00, 11) is

+ζ_(10,01) + ζ_(01,10) − ζ_(11,00),

and the replacement row is

ζ_(00,11) − (T/2) ζ_(10,01) − T ζ_(01,10).

The measurement row is +1, −1, −1, +1 on those four coordinates and 0 elsewhere. The placement puts X_0 at (00, 11) and X_1 at (10, 01), with the other two coordinates empty, so the two positive unit clauses become

X_0 − (T/2) X_1    and    X_1 − (T/2) X_0.

Summing the four assignments gives H(T) = 1 − 2T + T²/4. The same polynomial was summed by `eval_H_at`. At t = 2 both variables sit in group 0 while both clauses are assigned to group 1, which is empty, and the forms reduce to the measurement X_0 − X_1. Step 2's block design would refuse a clause that meets two blocks. This placement is legal for the kernel condition.

`cross-r3`. One positive 3-clause on three variables, each in its own group, the clause assigned to group 0. M0 is the exchange of (a, d) with (d, a) on all 9 coordinates. Variables 1 and 2 lie outside the assigned group. The degree is 0, so each replacement is still X_v. The three pulled-back measurement rows are −X_1, −X_2, and a zero row, so the two variables outside the assigned group are carried by integer measurements. H(T) = 7. The same clause is rejected by `validate_design` on a three-block popcount design, which is the Step 2 restriction this slice does not keep.

`r2-q3` uses the same formulas at r = 2, Q = 3, R = 2, where R < Q and the side length is 64. The schedules are ((0,1,0), (1,0,0), (0,0,1)) for the variables and (1,0,1) for the single clause x_0 ∨ x_1 ∨ ¬x_2. It was checked by the machine against brute force. It was not expanded by hand.

### Obstacle

Full M0 is still absent. The parameter list that accompanies that matrix asks for R ≥ 32 and 6R ≤ Q before a tail bound is invoked. The smallest such integers are R = 32 and Q = 192. At r = 2 the side length is 2^384. `build_auxiliary_slice` refuses every side length above 64 and says so in the error. A scheduler that would also keep groups of size at most 2n/r and label multiplicity at most 2M/D needs n on the order of D². That inequality is already false for every formula this repository compares with brute force, so the scheduler is not run. The four placements above are the simplification: the same matrix formulas, hand-written schedules, and direct checks. Nothing in this step runs in 2^{o(n)}.

## Step 4 — Hermite / multi-coefficient interpolation

Shipped in `hermite.py`. H is treated as a univariate polynomial of degree at most D = m·k·R, the same cap `degree_bound` and `interpolate_H0` use. The Lagrange path reads the value H(t) at t = 1, …, D+1. This step reads Hasse coefficients

H^{[k]}(t) = H^{(k)}(t) / k!

at fewer distinct t. The Newton form with repeated nodes takes those coefficients as its divided differences. Taylor coefficients are k! times the Hasse coefficients. The number of conditions stays D+1. A degree-D polynomial is not recovered from fewer conditions than that.

The jet is exact. For one assignment the linear forms are polynomials in T. Each clause factor Φ_e = 1 − ∏(1 − L_v) is evaluated in Q[u]/(u^{order+1}), u = T − t, and the jets are summed over the cube. One pass returns every requested coefficient at that t. Derivative orders are not separate passes. Arithmetic per assignment grows quadratically with the jet width, because each product is a truncated convolution.

Three patterns are run whenever they produce different multiplicity tuples:

- `paired`: at most two coefficients per point (the value and the first Hasse derivative). Point count ceil((D+1)/2). The last point may be value-only.
- `triple`: at most three coefficients per point (derivative orders 0, 1, and 2, and a shorter tail if the cap does not divide evenly). Point count ceil((D+1)/3).
- `single-point`: all D+1 coefficients at t = 1. One point. When D+1 ≤ 2 this is the same pattern as `paired`, and only one row is recorded.

Sample points are t = 1, 2, …, so 0 is not a node. H(0) is the interpolant. On width-1 hand formulas k = 1, so the cap is m·1·R. On every 3-CNF k = 3.

### What was checked

Every row below matched three ways: Hermite H(0), the direct sum `eval_H_at` at 0, and `interpolate_H0` on the same family and the same k. That common value equalled the brute-force model count, and it was an integer. The jet of order D at t = 1 reproduced `eval_H_at` at every integer from 0 through D+1. The degree read from that jet agreed with the degree read from forward differences of H(1), …, H(D+2), and the same degree was read from the jet at t = 2. The leading Hasse coefficient does not depend on the base point.

Three polynomials were expanded by hand, then compared with both `eval_H_at` and the jet:

- one negative unit, ⟨P, X⟩ = X(1+T), H(T) = 1 − T. At t = 1 the jet is (0, −1). #F = 1.
- two positive units, ⟨P, X⟩ = X_v + T(X_0 + X_1), H(T) = 1 + 6T + 6T^2. At t = 1 the jet is (13, 18, 6). The ordinary derivatives are (13, 18, 12). #F = 1. The first Hasse coefficient also equals the centered difference (H(2) − H(0))/2 = 18, which is exact because the degree is 2.
- (x) ∧ (¬x) with ⟨P, X⟩ = X(1+T) on both occurrences, H(T) = −T − T^2. At t = 1 the jet is (−2, −3, −1). #F = 0, and the jet is not zero, so the recovered count is not a constant-zero polynomial.

On every recorded instance the measured degree equalled the cap. The cap was tight for these families. No condition was dropped for a slack bound. That is one family per formula, not a proof for every replacement.

| instance | n | m | R | k | D | deg | #F | Lagrange points |
|----------|---|---|---|---|---|-----|----|-----------------|
| negative-unit | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 |
| unit-pair | 2 | 2 | 1 | 1 | 2 | 2 | 1 | 3 |
| contradiction-n1 | 1 | 2 | 1 | 1 | 2 | 2 | 0 | 3 |
| one-clause-n3 | 3 | 1 | 1 | 3 | 3 | 3 | 7 | 4 |
| random-n4-m4-r1 | 4 | 4 | 1 | 3 | 12 | 12 | 9 | 13 |
| random-n4-m3-r2 | 4 | 3 | 2 | 3 | 18 | 18 | 10 | 19 |
| random-n6-m6-r1 | 6 | 6 | 1 | 3 | 18 | 18 | 28 | 19 |
| bounded-n6-m4-r1 | 6 | 4 | 1 | 3 | 12 | 12 | 44 | 13 |

Multiplicities are the number of Hasse coefficients at t = 1, then t = 2, and so on. Their sum is the condition count, and it equals the Lagrange point count D+1. A blank triple cell means that pattern repeated single-point.

| instance | paired multiplicities | triple multiplicities | single-point |
|----------|-----------------------|-------------------------|--------------|
| negative-unit | | | (2,) — 1 point |
| unit-pair | (2, 1) — 2 points | | (3,) — 1 point |
| contradiction-n1 | (2, 1) — 2 points | | (3,) — 1 point |
| one-clause-n3 | (2, 2) — 2 points | (3, 1) — 2 points | (4,) — 1 point |
| random-n4-m4-r1 | (2, 2, 2, 2, 2, 2, 1) — 7 points | (3, 3, 3, 3, 1) — 5 points | (13,) — 1 point |
| random-n4-m3-r2 | (2, 2, 2, 2, 2, 2, 2, 2, 2, 1) — 10 points | (3, 3, 3, 3, 3, 3, 1) — 7 points | (19,) — 1 point |
| random-n6-m6-r1 | (2, 2, 2, 2, 2, 2, 2, 2, 2, 1) — 10 points | (3, 3, 3, 3, 3, 3, 1) — 7 points | (19,) — 1 point |
| bounded-n6-m4-r1 | (2, 2, 2, 2, 2, 2, 1) — 7 points | (3, 3, 3, 3, 1) — 5 points | (13,) — 1 point |

`one-clause-n3` is the positive 3-clause on n = 3 with replacement seed 0, the same family as the Step 0 single-clause test. The two `random-n4` rows and `random-n6-m6-r1` use the Step 0 seeds (formula seed 11, 13, 14 and replacement seed 111, 113, 114). `bounded-n6-m4-r1` is the Step 1 bench formula, replacement seed 42. Its model count 44 and degree 12 are the same figures as the Step 1 table. For `negative-unit`, D+1 = 2, so the paired pattern and the single-point pattern are the same one node.

Cube passes equal the number of distinct points. Assignments enumerated for that scheme are passes times 2^n. The Lagrange path enumerates (D+1)·2^n assignments and reads one value from each pass.

| instance | 2^n | Lagrange assignments | paired | triple | single-point | single-point jet width |
|----------|-----|----------------------|--------|----------|--------------|------------------------|
| negative-unit | 2 | 4 | | | 2 | 2 |
| unit-pair | 4 | 12 | 8 | | 4 | 3 |
| contradiction-n1 | 2 | 6 | 4 | | 2 | 3 |
| one-clause-n3 | 8 | 32 | 16 | 16 | 8 | 4 |
| random-n4-m4-r1 | 16 | 208 | 112 | 80 | 16 | 13 |
| random-n4-m3-r2 | 16 | 304 | 160 | 112 | 16 | 19 |
| random-n6-m6-r1 | 64 | 1216 | 640 | 448 | 64 | 19 |
| bounded-n6-m4-r1 | 64 | 832 | 448 | 320 | 64 | 13 |

`paired` roughly halves the number of cube passes and keeps the jet width at 2, so the extra arithmetic per assignment is a small constant. `single-point` uses one pass and a jet of width D+1. That pass does Θ(D²) rational operations per assignment. It is not a factor-of-D saving over the Lagrange path, and on a large cap it can be more arithmetic than D+1 scalar evaluations. The reproducible cost is the pair (passes, jet width). One run of `python -m counting_sat_algebra.verify --step4` marked every row OK. The slowest case on that run was `random-n6-m6-r1` at about half a second, and that figure includes the Lagrange comparison. It is not a complexity result.

### Obstacle

Producing the jet still scans {0,1}^n. Packing derivatives into a node reduces the number of distinct sample points and, for a short jet, the number of cube passes. It does not evaluate H below exponential cost in n, and it does not reduce the number of conditions below D+1. On these instances the degree sat on the cap, so there was no unused room in the bound either. Copying variables is Step 5, not part of this Hermite check. Theorem 1.1 is not implemented. Nothing in this step runs in 2^{o(n)}.

## Step 5 — occurrence copying

Shipped in `sparsify.py`. A variable that already occurs in at most B clauses stays one variable. A variable that occurs d > B times is replaced by c = ceil(d / (B − 2)) copies. Occurrences are dealt out in clause order, at most B − 2 per copy, and the copies are tied by the implication cycle

(¬y_1 ∨ y_2) ∧ (¬y_2 ∨ y_3) ∧ … ∧ (¬y_c ∨ y_1).

Each copy sits in two cycle clauses and at most B − 2 original clauses, so its occurrence is at most B. The smallest cap the gadget can meet is B = 3, which is the Step 1 cap and the default used below. A cap of 2 is refused: the cycle would leave no slot for an original clause. One pytest row uses B = 4 on six positive units (three copies, occurrence 4, #F' = #F = 1) to check that the same cycle with two original slots per copy still preserves the count. The CLI rows all use B = 3.

The cycle holds exactly when the copies are equal. Extending a source assignment by writing its bit on every copy is a bijection between source assignments and reduced assignments that satisfy the cycles, and it preserves each original clause. Therefore #F' = #F. There is no extra factor of 2. A copy is not a free variable. Pytest checks that bijection by enumerating both cubes whenever both lengths are at most 12, and checks the two cardinalities whenever both lengths are at most 20. `brute_count` refuses n > 20 before it loops. The same wall is used here.

For B = 3 the reduced length satisfies n' ≥ n and n' ≥ M − 2n, where M is the number of literal occurrences in the source. An exact 3-CNF therefore has n' ≥ 3m − 2n. The output is a mixed-width CNF: source clauses keep their width, and every cycle clause has width 2. `check_bounded_occurrence` still demands exact width 3, so it rejects a formula that used the cycle. The occurrence cap is the same number 3. Clause lists get longer by the number of cycle clauses. This is occurrence reduction, not a deletion of clauses, and it is one formula rather than a disjunction of many sparse formulas.

### Blowup

"Checked" means both model counts were computed and were equal. "Source only" means the source cube was counted and the reduced cube was not. M' is the number of literal occurrences after copying. Every checked row has M' ≤ 3n'.

| instance | n | n' | m | m' | occ | occ' | cycles | M' | #F | what was counted |
|----------|---|----|----|----|-----|------|--------|----|----|------------------|
| units-n1-m4 | 1 | 4 | 4 | 8 | 4 | 3 | 4 | 12 | 1 | both cubes |
| free-bit-n2 | 2 | 5 | 4 | 8 | 4 | 3 | 4 | 12 | 2 | both cubes |
| heavy-3cnf-n4 | 4 | 7 | 4 | 8 | 4 | 3 | 4 | 20 | 10 | both cubes |
| bounded-clause-n3 | 3 | 3 | 1 | 1 | 1 | 1 | 0 | 3 | 7 | both cubes; clauses unchanged |
| empty-clause-n2 | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | both cubes |
| two-block-copies | 3 | 9 | 8 | 16 | 4 | 3 | 8 | 26 | 2 | both cubes, then the blocked sum |
| cross-block-copies | 3 | 6 | 4 | 8 | 4 | 3 | 4 | 14 | 4 | both cubes; Step 2 refused |
| units-n1-m16 | 1 | 16 | 16 | 32 | 16 | 3 | 16 | 48 | 1 | both cubes; interpolation refused |
| units-n1-m21 | 1 | 21 | 21 | 42 | 21 | 3 | 21 | 63 | source 1 | source only |
| random-3cnf-n5-m12 | 5 | 36 | 12 | 48 | 8 | 3 | 36 | 108 | source 8 | source only |

`bounded-clause-n3` is one positive 3-clause. Occurrence is already 1, so the copying map is the identity and the cycle is empty. A separate pytest draw, `random_bounded_3cnf(6, 4, seed 7)`, is likewise unchanged, and `run_degree_bound_experiment` on that image agrees with brute force at both Step 1 caps. `free-bit-n2` is the four-unit formula plus a variable that occurs in no clause. Both sides have two models. The free bit is not removed and it is not duplicated.

`units-n1-m21` has reduced n = 21. `brute_count` raises before enumerating, and the equality is not checked. The source count is 1 because the single variable is forced true by the twenty-one positive units. The reduced side, twenty-one copies, was not enumerated. `random-3cnf-n5-m12` is `random_3cnf(5, 12, seed 0)`. Its source count is 8. The cap-3 bound gives n' ≥ 3·12 − 2·5 = 26, already past the wall, and the construction produced n' = 36. Every source variable occurred more than three times (maximum 8), so each occurrence became its own variable and n' = M = 36.

### Degree pipeline

Rows whose sample window fits the budget call `random_replacements` at R = 1 and `interpolate_H0` with k equal to the maximum clause width. The same samples are interpolated at the loose cap 3n'R, which is `degree_bound_3n` at B = 3. H(0) from both caps, the direct sum at 0, and the brute-force count of the reduced formula agreed, and one extra node agreed with the interpolant. "deg" is the degree read from forward differences of the longer window.

| instance | k | D = m'·k·R | M'·R | 3n'·R | deg | #F |
|----------|---|------------|------|-------|-----|----|
| units-n1-m4 | 2 | 16 | 12 | 12 | 12 | 1 |
| free-bit-n2 | 2 | 16 | 12 | 15 | 12 | 2 |
| heavy-3cnf-n4 | 3 | 24 | 20 | 21 | 20 | 10 |
| bounded-clause-n3 | 3 | 3 | 3 | 9 | 3 | 7 |
| empty-clause-n2 | 0 | 0 | 0 | 6 | 0 | 0 |
| cross-block-copies | 3 | 24 | 14 | 18 | 14 | 4 |

On these six draws the measured degree equalled M'·R. That is the bookkeeping bound deg H ≤ (sum of clause widths)·R. The width cap m'·k·R is at least that large, and it was strictly larger except on `bounded-clause-n3`, where the only clause has width 3, and on `empty-clause-n2`, where both the width cap and M'·R are 0 because the only clause is empty and H is the zero polynomial. The loose cap 3n'·R equalled M'·R on the four-unit row, where every reduced variable occurs three times, and was strictly larger on the other five. One family per formula. This is not a proof for every replacement.

`units-n1-m16` has n' = 16, which brute force still counts (#F' = 1). The interpolation plan is (66 + 64 + 2)·2^16 = 8,650,752 cube assignments: D = m'·k·R = 64 and the loose cap is 48, so the window has length 66. The budget is 100,000. The pipeline records the refusal and does not evaluate H. `two-block-copies` is counted and passed to Step 2. Its degree window would be D = 32 on n' = 9, which is inside the budget, and it is not part of this table. The blocked sum already checks H(0) on that formula.

A width-4 formula is in pytest and not in the CLI table. Four clauses on five variables, one variable in every clause and the others in three clauses each, become n' = 8 and m' = 8 with occurrence at most 3. The extension bijection holds. `interpolate_H0` was run at R = 0, where the cap is 0 and H is the constant #F'. The R = 1 window for that formula was not run.

One run of `python -m counting_sat_algebra.verify --step5` marked every row OK or RECORDED. The slowest interpolated row on that run was `heavy-3cnf-n4`, about 1.2 seconds, for 6,656 cube assignments. These times are one machine's measurement, not a complexity result. The reproducible cost is the pair (n', D) and the assignment count.

### Step 2, when the clauses already sit in blocks

`two-block-copies` partitions the source variables into {0, 1} and {2}. Clauses use only one of those sets. After copying, the blocks are (0, 1, 2, 3, 4) and (5, 6, 7, 8): four copies of variable 0, the unsplit variable 1, and four copies of variable 2. Cycle clauses mention copies of a single source variable, so they stay inside that block. The design is the Step 2 popcount design, τ = 1, R = 1, every higher coefficient 1. Blocked H(t) matched the full-cube sum at t = 0 and at t = 1, and H(0) = #F' = 2.

|  | value |
|--|-------|
| r | 5 |
| blocks | 2 |
| measurement tuples | 30 |
| local evaluations | 256 |
| full cube | 512 |
| H(0) | 2 |
| H(1) | −47986541792015699006680248 |

H(1) is a value of the counting polynomial at t ≠ 0. It is not a model count. The tuple count is (5 + 1)·(4 + 1) = 30, the same exponential-in-the-number-of-blocks enumeration Step 2 already records. On this row the local-evaluation count, 256, is below the cube, 512. That comparison is one parameter point. It is not a subexponential algorithm.

`cross-block-copies` uses blocks {0, 1} and {2} and a single 3-clause on variables 0, 1, and 2, plus three positive units on variable 0. Copying produces n' = 6 and #F' = #F = 4, and `interpolate_H0` recovers 4. The 3-clause still meets both blocks. `validate_design` refuses it with the Step 2 message that a clause must sit inside one block. The cycle does not repair a clause that already crossed. Step 2 is not run on a general copied formula.

### Obstacle

Step 3 is not applied. Its schedules are written out by hand for four placements, and this reduction does not produce a schedule. The cycle clauses have width 2. The side-length cap remains 64, and full M0 at R ≥ 32 and 6R ≤ Q is still not built.

The reduction does not keep n' close to n. At B = 3 a variable in d clauses becomes d variables whenever d > 3, and a 3-CNF with 3m − 2n > 20 leaves the brute-force window no matter how the clauses are arranged. The n = 5, m = 12 row is that situation at a tiny source. Nothing in this step runs in 2^{o(n)}. Theorem 1.1 is not implemented.

## Step 6 — parameter scan

Shipped in `parameter_scan.py`. The grid is fixed. Pytest and `verify --quick` run the `quick` profile (22 cells). `verify --step6` and `benches/parameter_scan.json` run the `full` profile (33 cells), which includes every quick cell and then the larger n. Nothing on the grid allocates full M0. Nothing on the grid enumerates a cube with n > 20. A cell whose pipeline was not started has no pipeline timing. The clocks below are one run of `python -m counting_sat_algebra.verify --step6` on the agent machine. The bench file is another run of the same grid. The two clocks are not required to agree. The operation counts are.

The popcount family is the one from Step 2: b blocks of size r, one clause per block, every higher coefficient 1. For r = 2 the clause forbids one of the four local assignments, the blocks are disjoint, and #F = 3^b. For r = 3, #F = 7^b. For r = 1 the clauses are units of alternating sign and #F = 1. H(0) does not depend on R, because the replacement restores X_v at T = 0. D = b · min(r, 3) · R is stored on the design and is not read by the blocked sum.

### Blocked sum against brute force

"Agreed" means blocked H(0), the full-cube sum `eval_H_at` at 0, and `brute_count` were equal. "Refused" means that comparison was not finished. The plan columns are the closed forms (r+1)^b, b · 2^r · (r+1)^{b-1}, and 2^n. On a refused row those counts were not performed, except the brute-force column of the tuple-cap row, which did walk the cube.

| cell | n | r | b | R | #F | tuples | local | cube | what ran |
|------|---|---|---|---|----|--------|-------|------|----------|
| blocked-r2-b1 | 2 | 2 | 1 | 1 | 3 | 3 | 4 | 4 | agreed |
| blocked-r2-b2 | 4 | 2 | 2 | 1 | 9 | 9 | 24 | 16 | agreed |
| blocked-r2-b3 | 6 | 2 | 3 | 1 | 27 | 27 | 108 | 64 | agreed |
| blocked-r2-b4 | 8 | 2 | 4 | 1 | 81 | 81 | 432 | 256 | agreed |
| blocked-r2-b5 | 10 | 2 | 5 | 1 | 243 | 243 | 1620 | 1024 | agreed |
| blocked-r2-b6 | 12 | 2 | 6 | 1 | 729 | 729 | 5832 | 4096 | agreed |
| blocked-r2-b7 | 14 | 2 | 7 | 1 | 2187 | 2187 | 20412 | 16384 | agreed |
| blocked-r2-b2-R2 | 4 | 2 | 2 | 2 | 9 | 9 | 24 | 16 | agreed; D = 8, tuples unchanged |
| blocked-r1-b4 | 4 | 1 | 4 | 1 | 1 | 16 | 64 | 16 | agreed; tuples = 2^n |
| blocked-r3-b3 | 9 | 3 | 3 | 1 | 343 | 64 | 384 | 512 | agreed; local evaluations below the cube |
| blocked-r1-b17-tuple-cap | 17 | 1 | 17 | 1 | 1 | 131072 | 2228224 | 131072 | brute force only; blocked sum refused |
| blocked-r2-b11-brute-wall | 22 | 2 | 11 | 1 | — | 177147 | 2598156 | 4194304 | neither enumeration started |

From b = 2 upward on the r = 2 rows, the local-evaluation count is already larger than the cube (24 > 16, 108 > 64, and so on, through 20412 > 16384). The r = 3 row goes the other way: 384 local evaluations and a cube of 512. Either way the tuple count is exponential in the number of groups. For r = 2 it is 3^{n/2}.

Wall times from that one `--step6` run, in seconds. "Blocked" is the `blocked_sum` call. "Full" is `eval_H_at` at t = 0. "Brute" is `brute_count`. A dash means that call was not made.

| cell | brute | blocked | full sum |
|------|-------|---------|----------|
| blocked-r2-b1 | <0.0001 | 0.0001 | <0.0001 |
| blocked-r2-b2 | <0.0001 | 0.0004 | 0.0004 |
| blocked-r2-b3 | <0.0001 | 0.0013 | 0.0025 |
| blocked-r2-b4 | 0.0001 | 0.0052 | 0.0146 |
| blocked-r2-b5 | 0.0006 | 0.0190 | 0.0789 |
| blocked-r2-b6 | 0.0028 | 0.0698 | 0.3892 |
| blocked-r2-b7 | 0.0109 | 0.2448 | 1.9374 |
| blocked-r3-b3 | 0.0004 | 0.0063 | 0.0489 |
| blocked-r1-b17-tuple-cap | 0.0557 | — | — |
| blocked-r2-b11-brute-wall | — | — | — |

On this run the blocked loop's clock sat below the rational cube sum, and brute force sat below both. That is not an ordering of the operation counts. Boolean model counting, the blocked rational loop, and `eval_H_at` do different work per assignment. The local-evaluation count is still larger than 2^n on the r = 2 rows with b ≥ 2. These times are not a complexity result.

The tuple-cap row is inside the brute-force wall and outside the measurement cap of 100000. `blocked_sum` raises before the measurement loop. The closed-form local-evaluation count was not performed. The brute-wall row is past both rails: `brute_count` refuses n = 22 before it loops, and the blocked sum is not called, because a parameter point with n > 20 and a small tuple count would otherwise start the measurement loop with nothing to compare it to.

### Where the other steps match, and where they stop

Same popcount family, R = 1, r = 2, unless the cell says otherwise. The assignment budget is 100000 cube assignments, the rail from Step 5. Lagrange's plan is (D+2)·2^n: the D+1 nodes of `interpolate_H0`, plus one direct evaluation at 0. Paired Hermite's plan is (passes + 1)·2^n, with passes = ceil((D+1)/2) and jet width 2, plus the direct evaluation. The single-point jet of order D is a separate cell and is not run.

| cell | n | D | plan | #F | what ran |
|------|---|---|------|----|----------|
| lagrange-r2-b1 | 2 | 2 | 16 | 3 | agreed |
| lagrange-r2-b2 | 4 | 4 | 96 | 9 | agreed |
| lagrange-r2-b3 | 6 | 6 | 512 | 27 | agreed, 0.0276s |
| lagrange-r2-b4 | 8 | 8 | 2560 | 81 | agreed, 0.2212s |
| lagrange-r2-b5 | 10 | 10 | 12288 | 243 | agreed, 1.5041s |
| lagrange-r2-b7-budget | 14 | 14 | 262144 | — | refused, not timed |
| hermite-paired-r2-b1 | 2 | 2 | 12 | 3 | agreed; 2 points, jet width 2, 3 conditions |
| hermite-paired-r2-b2 | 4 | 4 | 64 | 9 | agreed; 3 points, 0.0025s |
| hermite-paired-r2-b3 | 6 | 6 | 320 | 27 | agreed; 4 points, 0.0220s |
| hermite-paired-r2-b7-budget | 14 | 14 | 147456 | — | refused, not timed |
| hermite-single-point-not-run | 2 | 2 | 4 | — | not started; jet width would be 3 |

At n = 14 the blocked sum still runs (2187 tuples, under the measurement cap, n ≤ 20) and both interpolation plans are over the budget. That split is the scan: the steps do not share a cost. Paired Hermite on the rows that ran kept D+1 conditions and fewer distinct points than the value-only path. Each of those points is still one pass over the cube. The single-point cell records the plan for the b = 1 formula, one pass of 4 assignments with jet width D+1 = 3, and does not run it. Packing every coefficient into t = 1 does not evaluate H below exponential cost in n.

The degree cells use the Step 1 sampler, window max(m·3·R, 3·n·R)+2.

| cell | n | m | caps | plan | #F | degree | what ran |
|------|---|---|------|------|----|--------|----------|
| degree-clause-n3 | 3 | 1 | 3 and 9 | 88 | 7 | 3 | agreed, 0.0020s |
| degree-bounded-n6-m4 | 6 | 4 | 12 and 18 | 1280 | 44 | 12 | agreed, 0.1329s |
| degree-bounded-n12-m12-budget | 12 | 12 | 36 and 36 | 155648 | — | — | formula built, H not sampled, not timed |

`degree-clause-n3` is one positive 3-clause. The measured degree equalled the tight cap. `degree-bounded-n6-m4` is `random_bounded_3cnf(6, 4, seed 42)` with replacement seed 42, the same counts as the Step 1 bench row n6-m4-r1: #F = 44 and degree 12. The n = 12 row is a legal occurrence-≤3 3-CNF. The window is past the budget, so brute force was not run either.

Copying, cap B = 3. "Source #F" is the cube that was counted. The reduced cube is counted only when n' ≤ 20.

| cell | n | n' | source #F | what ran |
|------|---|----|-----------|----------|
| copy-popcount-b1 | 2 | 2 | 3 | agreed; occurrence 1, the map is the identity, #F' = #F |
| copy-blowup-n5-m12 | 5 | 36 | 8 | reduced cube not enumerated; lower bound n' ≥ 26; occurrence after copying is 3 |
| copy-units-n21 | 1 | 21 | 1 | reduced cube not enumerated; `brute_count` refuses n = 21 |
| copy-random-n4-m8 | 4 | 24 | 4 | reduced cube not enumerated; lower bound n' ≥ 16 |
| copy-cross-block | 3 | 6 | 4 | #F' = #F = 4 was counted; the blocked sum was refused |

`copy-blowup-n5-m12` is `random_3cnf(5, 12, seed 0)`, the same formula as Step 5. The plan for the reduced cube is 2^36 assignments. That loop was not started. `copy-cross-block` is the Step 5 crossing clause. After copying, n' = 6 and the 3-clause still meets both blocks. `run_step2_pipeline` returns the Step 2 refusal. The scan does not call `place_formula` and does not build a schedule. The blocked sum was not timed. Step 3's kernel condition, on its four hand-written placements, is a different construction. This reduction does not produce one of those placements.

The matrix cells:

| cell | what ran |
|------|----------|
| matrix-swap-4 | agreed. Side length 4, the hand-written exchange, n = 2, #F = 1, observed degree 0. The degree report walked 20 cube assignments. This is not full M0 |
| matrix-full-m0 | refused before allocation. `build_auxiliary_slice(2, 192, 32)` raises. The side length is 2^384, above the slice cap 64. No scheduler is invoked. No timing was recorded |

### Distance to a full counting algorithm

The scan is a record of the pieces that exist. It is not a claim that the missing pieces are close, and it is not a 2^{o(n)} algorithm. The unfinished list:

- Full M0 at R ≥ 32 and 6R ≤ Q is not built. The smallest such integers are R = 32 and Q = 192. At r = 2 the side length is 2^384. `build_auxiliary_slice` refuses that call before it allocates the matrix. The scan records the refusal and does not time a matrix.
- A scheduler that would place variables and clauses into those coordinates, with group size at most 2n/r and label multiplicity at most 2M/D, is not implemented. The four Step 3 placements are hand-written. That size bound needs n on the order of D^2, which is already false for every formula this repository compares with brute force, so the scheduler is not run.
- Sparsification that keeps n' = n + o(n) is not implemented. Implication-cycle copying at B = 3 replaces a variable of occurrence d > 3 by d variables. The scan records `random_3cnf(5, 12, seed 0)` at n' = 36, past the brute-force wall n ≤ 20, and does not enumerate that reduced cube.
- A subexponential measurement space is not implemented. The popcount design enumerates (r+1)^b tuples. For blocks of size 2 that is 3^{n/2}. On the rows this scan runs, the local-evaluation count is often larger than the cube. The scan does not claim a 2^{o(n)} measurement loop.
- Hermite interpolation does not remove the 2^n pass. The paired scheme packs D+1 conditions into fewer distinct t, and each t still sums the cube. A full jet of order D does quadratic arithmetic per assignment. The scan runs the paired scheme only inside the assignment budget. It does not run the single-point jet as a substitute for the cube, and it does not treat either scheme as a subexponential evaluation.
- Clauses that cross a Step 2 block partition are not factored by the blocked sum. Copying does not pull a crossing clause into one block, and it does not produce a Step 3 placement. The scan records that refusal and does not time a blocked sum on it.
- Theorem 1.1 is not implemented. This repository does not refute ETH. No cell of the scan is a 2^{o(n)} counting algorithm.

There is no further step in `ROADMAP.md`. The list above is the remaining work.

## Known gaps

Step 6 records these on its grid. Recording them does not remove them.

- Full M0, at the parameter list R ≥ 32 and 6R ≤ Q, is not built. The matrix in this tree has side length at most 64. Step 5 does not hand the copied formula to that matrix. The scan refuses the call and does not time it.
- The occurrence reduction checks #F' = #F inside the brute-force window and records the rows that leave it. It is not a reduction with n' = n + o(n), and it does not remove the 2^{n'} pass inside `interpolate_H0` or the 2^n pass inside each Hermite sample.
- The copied formula is mixed width. `check_bounded_occurrence` accepts it only when the reduction added no cycle, which is the case already inside the Step 1 testbed. A general copied instance is not an exact 3-CNF with m ≤ n.
- Step 2 matches the full sum on the instances in the Step 2 section and on `two-block-copies`. It refuses a copied clause that still crosses blocks. Measurement enumeration for the popcount design is (r + 1)^b, exponential in the number of groups, and on the recorded r = 2 scan it does more local work than summing the cube.
- The Step 0/1 Lagrange path is still Θ(2^n · D) scalar evaluations. The Step 4 paired path is about half as many cube passes, with a jet of width 2. The single-point path is one pass with Θ(D²) arithmetic per assignment. The Step 5 degree rows are the same shape of cost on n' and on m'·k·R. None of these is subexponential.
- The Step 1 degree table is the bench set plus one pytest draw, all at R = 1. Step 3 has its own four-row degree table. The Step 4 degree column is the eight rows in that section; on each of them the measured degree equalled the cap. The Step 5 degree column is the six rows above; on each of them the measured degree equalled M'·R.
- Self-reduction deletes fixed variables from clauses and then calls brute force on all 2^n bits, so the integer it sees is not the number of extensions. The CLI only needs a yes/no answer and checks the witness it returns.
