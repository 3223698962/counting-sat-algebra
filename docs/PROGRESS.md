# Progress

Steps 0, 1, 2, a small-parameter Step 3, and an exact Step 4 Hermite check are implemented and compared with brute-force counting, or with the full-space sum H(t), on small instances. This repository does not refute ETH. It does not implement sparsification or Theorem 1.1. The matrix M0 is built only at side length r**(2Q) ≤ 64. The Step 2 measurement loop is exponential in the number of groups. Step 4 uses fewer distinct sample points than D+1 and still enumerates the cube at each point.

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

Producing the jet still scans {0,1}^n. Packing derivatives into a node reduces the number of distinct sample points and, for a short jet, the number of cube passes. It does not evaluate H below exponential cost in n, and it does not reduce the number of conditions below D+1. On these instances the degree sat on the cap, so there was no unused room in the bound either. Sparsification is not implemented. Theorem 1.1 is not implemented. Nothing in this step runs in 2^{o(n)}.

## Next: Step 5

Counting-preserving sparsification, or variable copying, at small n, and a path from a general small k-CNF into a bounded-occurrence instance. That step is not in this tree.

## Known gaps

- Full M0, at the parameter list R ≥ 32 and 6R ≤ Q, is not built. The matrix in this tree has side length at most 64.
- No sparsification and no Theorem 1.1. Step 4 does not remove the 2^n pass inside each sample point.
- Step 2 matches the full sum on the instances above. Measurement enumeration for the popcount design is (r + 1)^b, exponential in the number of groups, and on the recorded r = 2 scan it does more local work than summing the cube.
- The Step 0/1 Lagrange path is still Θ(2^n · D) scalar evaluations. The Step 4 paired path is about half as many cube passes, with a jet of width 2. The single-point path is one pass with Θ(D²) arithmetic per assignment. Neither is subexponential.
- The Step 1 degree table is the bench set plus one pytest draw, all at R = 1. Step 3 has its own four-row degree table. The Step 4 degree column is the eight rows above; on each of them the measured degree equalled the cap.
- Self-reduction deletes fixed variables from clauses and then calls brute force on all 2^n bits, so the integer it sees is not the number of extensions. The CLI only needs a yes/no answer and checks the witness it returns.
