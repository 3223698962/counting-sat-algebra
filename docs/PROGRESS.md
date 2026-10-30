# Progress

Version **0.28.0** follows the newly supplied Lean SqrtLogSpace exact-count program. Step 23 implements compatible binary syntax, exact dyadic parameter selection and a streamed disjoint counting frontier on the original declared n. Step 24 constructs ternary indicator digits and exact accepted/marked queries, with a native algebraic leaf connected to that frontier. Step 25 builds the scheduled sparse ternary code with regenerated rows and selects history rows by violation mass scored with real marked queries. Step 26 streams the high-incidence rank loop and lets the hash leaf route every query through it. The finite Sections 3–6 and Appendix D milestones through Step 22 remain available below. Step 27 constructs the counting extension ring, the fixed nodes and the Q*J > degree Hermite placement. Step 28 streams the correction pairs and re-reads the per-node Taylor coefficient blocks by Newton Hermite. Step 29 streams the pointwise correction-pair leaf without any squarefree support table. Step 30 assembles and places one leaf at the recorded ports and threads it through one run.

The supplied Lean library contains final theorem proof source and passes the static local-import audit. Local kernel compilation and Python refinement are pending. Its time exponent retains the square of log log n; Appendix D's stronger target is a later route. The complete Python asymptotic pipeline is not implemented. `theorem_1_1_complete`, `lean_proof_complete` and `recursive_contraction` stay false. Current obligations are in [REFERENCE_ALIGNMENT.md](REFERENCE_ALIGNMENT.md) and [ROADMAP.md](../ROADMAP.md).

The Step 0–8 sections below preserve historical measurements. Step 7/8 Boolean-budget refusals are reproduced through the explicit legacy mode. They no longer describe the default scheduler. The older Step 6 distance list is a historical snapshot, updated below by the new reference alignment.

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

Full M0 is still absent. The parameter list that accompanies that matrix asks for R ≥ 32 and 6R ≤ Q before a tail bound is invoked. The smallest such integers are R = 32 and Q = 192. At r = 2 the side length is 2^384. `build_auxiliary_slice` refuses that call before it allocates the matrix. The four placements above are the Step 3 simplification: the same matrix formulas, hand-written schedules, and direct checks. Step 7 raises the cap from 64 to 256 and adds the balance-check scaffold. It still does not choose those schedules. Nothing in this step runs in 2^{o(n)}.

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

Step 3 is not applied. Its schedules are written out by hand for four placements, and this reduction does not produce a schedule. The cycle clauses have width 2. Step 7 later raises the side-length cap from 64 to 256; this reduction still does not place a formula on that matrix. Full M0 at R ≥ 32 and 6R ≤ Q is still not built.

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
| matrix-full-m0 | refused before allocation. `build_auxiliary_slice(2, 192, 32)` raises. The side length is 2^384, above the slice cap. Step 7 later set that cap to 256; the refusal is unchanged. No timing was recorded |

### Historical distance recorded by the Step 6 scan

The scan is a record of the pieces that exist. It is not a claim that the missing pieces are close, and it is not a 2^{o(n)} algorithm. The unfinished list:

- Full M0 at R ≥ 32 and 6R ≤ Q is not built. The smallest such integers are R = 32 and Q = 192. At r = 2 the side length is 2^384. `build_auxiliary_slice` refuses that call before it allocates the matrix. The scan records the refusal and does not time a matrix.
- A Step 7 scheduler scaffold checks the balance inequalities 2(b+w)D² ≤ M and 2Qr² ≤ n. They need n on the order of D². They fail for every formula this repository compares with brute force, and the scaffold then refuses: it does not choose schedules and it does not invent a random placement. Step 8 runs a successive-choice search when the inequalities hold and the slice is under the cap. On 64 disjoint units at r = 2, Q = 1, D = 4 the search meets label multiplicity 16 ≤ 2M/D = 32 and group size 32, then refuses: `design_problems_at` would enumerate 2^32 Boolean points, above the budget 4096, so `place_formula` is not called. The schedule cube is not enumerated. That search is not a proof of Lemma 4.5. Caller-supplied schedules can still be checked, including the hand-written Step 3 placements, against group size at most 2n/r and, as a recorded bound, label multiplicity at most 2M/D.
- Sparsification that keeps n' = n + o(n) is not implemented. Implication-cycle copying at B = 3 replaces a variable of occurrence d > 3 by d variables. The scan records `random_3cnf(5, 12, seed 0)` at n' = 36, past the brute-force wall n ≤ 20, and does not enumerate that reduced cube.
- A subexponential measurement space is not implemented. The popcount design enumerates (r+1)^b tuples. For blocks of size 2 that is 3^{n/2}. On the rows this scan runs, the local-evaluation count is often larger than the cube. The scan does not claim a 2^{o(n)} measurement loop.
- Hermite interpolation does not remove the 2^n pass. The paired scheme packs D+1 conditions into fewer distinct t, and each t still sums the cube. A full jet of order D does quadratic arithmetic per assignment. The scan runs the paired scheme only inside the assignment budget. It does not run the single-point jet as a substitute for the cube, and it does not treat either scheme as a subexponential evaluation.
- Clauses that cross a Step 2 block partition are not factored by the blocked sum. Copying does not pull a crossing clause into one block, and it does not produce a Step 3 placement. The scan records that refusal and does not time a blocked sum on it.
- Theorem 1.1 is not implemented. This repository does not refute ETH. No cell of the scan is a 2^{o(n)} counting algorithm.

The original Step 7/8 paths preserve these stops in legacy Boolean mode. Steps 9–15 update the default design, common evaluator, auxiliary representation, supplied-CNF entry point, replacement-support contraction and finite connected Taylor elimination; their remaining obligations are listed below.

## Step 7 — larger slice and a balance-check scaffold

Shipped in `scheduler.py`, with the cap constant in `replacement_matrix.py`. The slice cap is 256. Two explicit placements that the cap of 64 refused were checked with the same identities as Step 3: `auxiliary_problems`, T = 0 restoration of X_v, degree at most R−1, the support bound on nonzero M0 entries, the measurement rows, the kernel condition, H(0) against brute force, and both degree caps. Schedules are frozen data. They are not the output of a balancer. Theorem 1.1 is not implemented.

The cap is 256 because that is r**(2Q) for r = 2, Q = 4, and the eigenbasis walk in `auxiliary_problems` finishes on one machine. The next pure powers are 729 (r = 3, Q = 3) and 1024 (r = 2, Q = 5). Both are refused before a matrix is stored. Full M0 stays at 2^384 and is refused before allocation, including when the caller passes a larger `max_dimension`. `verify --quick` runs the side-length-81 point and the refusals. `verify --step7` also runs side length 256. The clock on that check is not a complexity result.

### What was checked

Both points use one width-3 clause on n = 3. `r3-q2` is the positive clause x_0 ∨ x_1 ∨ x_2. `r2-q4` is x_0 ∨ x_1 ∨ ¬x_2, the same clause as `r2-q3`. R < Q holds at both points. The interval is the default one.

| instance | r | Q | R | D | C rows | n | m | #F | deg H | M(R−1) | m·k·R | t checked |
|----------|---|---|---|---|--------|---|---|----|-------|--------|-------|-----------|
| r3-q2 | 3 | 2 | 1 | 81 | 45 | 3 | 1 | 7 | 0 | 0 | 3 | 1, 2 |
| r2-q4 | 2 | 4 | 3 | 256 | 13 | 3 | 1 | 7 | 1 | 6 | 9 | 1, 2, 3, 4 |

`r3-q2` has degree 0, so H is the constant 7. `r2-q4` has measured degree 1 under a cap of 6. At t = 1 the `r3-q2` groups are (0, 1), (2,), and an empty third group, and the clause is assigned to group 0, so variable 2 sits outside that group. At t = 4 the `r2-q4` groups are (0, 1) and (2,), and the clause is assigned to group 0. The kernel condition was checked at every t, not only at those two.

Sample values, all equal to the full-cube sum:

| instance | H(0) | H(1) | H(2) | other |
|----------|------|------|------|-------|
| r3-q2 | 7 | 7 | 7 | |
| r2-q4 | 7 | 29/4 | 15/2 | H(3) = 31/4, H(4) = 8 |

`r3-q2` schedules are ((0, 0), (0, 0), (1, 1)) for the variables and (0, 0) for the clause. `r2-q4` schedules are ((0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 1)) for the variables and (0, 0, 0, 0) for the clause. Label multiplicity is 2 on both. The bound 2M/D is 6/81 and 6/256, so both placements sit above it. The scaffold records that and still accepts the caller-supplied schedules, because `design_problems_at` and the T = 0 restoration passed. It does not call them a balancer output.

### Balance checks

The inequalities are the Lemma 4.5 package used in the obstacle above: n ≤ M ≤ bn, 2(b+w)D² ≤ M, and 2Qr² ≤ n. Here b is the largest number of clause-scopes containing one variable, w is the largest scope, and M is the sum of the scope sizes. A repeated variable inside one clause counts once. With b ≥ 1 the two displayed inequalities force n ≥ 2(b+w)D²/b ≥ 2D², which is the “n on the order of D²” bound. Lemma 4.7 also asks for 1 ≤ R < Q. That gap is recorded. The binomial tail in that lemma is not checked.

On the two unit clauses of `swap-4`, r = 2, Q = 1, D = 4, one has b = w = 1, M = 2, n = 2. Then 2(b+w)D² = 64 > 2 and 2Qr² = 8 > 2, and the weak order bound is n ≥ 32. The call with no schedules returns no placement. The same refusal is what the four Step 3 formulas produce. Passing the hand-written schedules anyway places them and `design_problems_at` passes. `swap-4` meets multiplicity 2M/D = 1. `remark-4-4` has 2M/D = 1/4 and multiplicity 1, so that bound fails and the design check still passes.

Thirty-two disjoint unit clauses give n = M = 32, b = w = 1, and the same r, Q, D. Then n = 2D², so the weak order bound holds, and 2(b+w)D² = 64 > 32, so the mass inequality still fails. The factor (b+w)/b is 2, and the tight necessary n is 64. Sixty-four disjoint unit clauses meet every Lemma 4.5 inequality at these parameters: 2(b+w)D² = 64 ≤ M and 2Qr² = 8 ≤ n. R = 1 is not strictly below Q = 1, so the Lemma 4.7 degree gap fails. Step 8 runs the successive-choice search on that call. The search meets the multiplicity and group bounds and does not call `build_auxiliary_slice`: the design walk would be 2^32 Boolean points, above the budget. It does not call `brute_count`. n = 64 is past the brute-force wall of 20. That cube was not enumerated. A passing numerical check is not a placement.

`schedule_placement` on r = 2, Q = 192, R = 32 refuses before allocation. The side length it reports is 2^384. `build_auxiliary_slice` raises on the same parameters, including when `max_dimension` is set far above 256, because Q > 30 is refused before the matrix is stored. No timing is attached to either refusal.

### Obstacle

Step 7 records the balance refusal and the two larger slices. Choosing schedules when the inequalities hold is Step 8. The binomial tail bound that Lemma 4.7 uses to control the number of measurement rows is not checked here. Full M0 is not built. Sparsification with n' = n + o(n) is not implemented. Nothing in this step runs in 2^{o(n)}.

## Step 8 — conditional-expectation balancer scaffold

Shipped in `scheduler.py`. When `check_balance` admits a formula and the slice is under the cap, `schedule_placement` chooses variable schedules and clause schedules by successive choice over the alphabet {0, …, r−1}^Q. At each variable, and then at each clause, every alphabet symbol is scored and the minimum is kept. The variable score spreads symbols and refuses a symbol that would put a group above size floor(2n/r), read as `size * r <= 2n`. The clause score is the increase in the sum of squares of the label loads. That increase is the one-step conditional expectation of the quadratic potential. Ties break toward the lexicographically first symbol. Nothing is sampled. The schedule cube, of size (r^Q)^(n+m), is not enumerated.

The search then compares the schedules with label multiplicity at most 2M/D and with the group bound. If either fails, the call stops and does not call `place_formula`. If both hold and the Boolean walk inside `design_problems_at` is at most 4096 points, the schedules are placed and `design_problems_at` and the T = 0 restoration are run. A supplied R that fails 1 ≤ R < Q blocks acceptance, and so does Q = 1, where no positive R is strictly below Q. If the Boolean walk is larger than the budget, `design_problems_at` is not called and `place_formula` is not called. The alphabet is not enumerated when r^Q is above 32, and the successive loop is not started when the candidate count (n+m)·r^Q is above 100000. Those two rails are not hit by an admitted slice under the cap of 256: there r^Q is at most 16. The large-n rail is what refuses an admitted formula whose n makes the candidate count too big, before the loop.

This is a machine-limited scaffold. It is not a proof of Lemma 4.5. It is not Theorem 1.1. It does not refute ETH. It is not a 2^{o(n)} count. Refusals are not timed.

### What was checked

`verify --step8` and the same rows inside `verify --quick` recorded the following. The five refusals were not timed. The direct check was.

Sixty-four disjoint positive units, r = 2, Q = 1, D = 4. This is the admitted point from Step 7. b = w = 1, M = n = 64, so 2(b+w)D² = 64 ≤ M and 2Qr² = 8 ≤ n. The alphabet has 2 symbols. The search considered 256 candidate symbols. Label multiplicity is 16, and 2M/D = 32. The largest group has size 32, and 32·r = 64 ≤ 2n = 128. Both clause symbols are used, so an outside group also has size 32. `design_problems_at` would enumerate 2^32 Boolean points for that group on each occurrence that reads it. That is above the budget 4096. `place_formula` was not called. No placement was accepted. The call with R = 1 also records that Lemma 4.7 asks for 1 ≤ R < Q. The call with R omitted records that no integer R satisfies that gap at Q = 1. Either sentence is extra. The stop is the design budget. n = 64 is past the brute-force wall of 20. That cube was not enumerated.

One thousand twenty-four disjoint positive units, r = 2, Q = 2, R = 1, D = 16. The degree gap holds. b = w = 1, M = n = 1024, so 2(b+w)D² = 1024 ≤ M and 2Qr² = 16 ≤ n. The alphabet has 4 symbols. The search considered 8192 candidate symbols. Label multiplicity is 64, and 2M/D = 128. The largest group has size 512, and 512·r = 1024 ≤ 2n = 2048. An outside group has size 512. The design walk would be 2^512 Boolean points. The budget refuses before `place_formula`. The side length 16 is under the cap, and the matrix was still not allocated, because the design check was not going to run.

Two disjoint units at the same r, Q, and R = 1 are not admitted: 2(b+w)D² = 64 > M = 2. The balancer was not run. No placement was invented. The same refusal is what a single positive 3-clause produces. Side length 1024, from r = 2, Q = 5, and full M0 at R = 32, Q = 192, side length 2^384, both refuse before allocation. The search was not started. No timing was attached.

Eight disjoint units at r = 2, Q = 2, R = 1 are not admitted either. D = 16 and 2(b+w)D² = 1024 > M = 8, and 2Qr² = 16 > n = 8. `schedule_placement` does not run the search. A direct call to the same placement check, which is not a gate for the inequalities, chose schedules with multiplicity 1 ≤ 2M/D = 1 and an outside group of size 4. The Boolean-point total is 256, inside the budget. `design_problems_at` and the T = 0 restoration passed, and the call built a slice of side length 16. One run of `verify --step8` timed that direct check at about 0.01 seconds. The clock is one machine's measurement. It is not a complexity result, and the row is not an admitted instance. The reason on that call says the balance inequalities were not the gate.

### Obstacle

On the admitted rows this repository runs, the outside group is large enough that `design_problems_at` would walk more than 4096 Boolean points. The documented toy is the smallest disjoint-unit point at r = 2, Q = 1, and its outside group has size 32. On that point, meeting multiplicity at most 2M/D forces a split: putting every variable in one group would put more than 2M/D occurrences on one label, and any clause symbol then leaves an outside group of size at least 32. The search meets the Lemma 4.5 bounds it was aimed at, and then refuses the design check instead of starting the walk. A smaller unadmitted formula can pass the design check. That does not make the inequalities hold, and `schedule_placement` does not place it.

The binomial tail bound is not checked. Full M0 is not built. Sparsification with n' = n + o(n) is not implemented. The measurement enumeration of Step 2 is still exponential in the number of groups. Hermite interpolation still walks the cube. Nothing in this step runs in 2^{o(n)}.

## Historical gaps through Step 8

Step 6 records these on its grid. Step 7 raises the slice cap and adds the balance refusal. Step 8 runs the successive-choice search on the admitted rows and refuses them before `design_problems_at`. Recording them does not remove them.

- Full M0, at the parameter list R ≥ 32 and 6R ≤ Q, is not built. The matrix in this tree has side length at most 256. Step 5 does not hand the copied formula to that matrix. The scan, the Step 7 scaffold, and the Step 8 balancer all refuse the call and do not time it. The balancer meets the multiplicity and group bounds on the admitted rows it runs and does not place them: the design walk is above its budget. It is not a proof of Lemma 4.5, and the binomial tail is not checked.
- The occurrence reduction checks #F' = #F inside the brute-force window and records the rows that leave it. It is not a reduction with n' = n + o(n), and it does not remove the 2^{n'} pass inside `interpolate_H0` or the 2^n pass inside each Hermite sample.
- The copied formula is mixed width. `check_bounded_occurrence` accepts it only when the reduction added no cycle, which is the case already inside the Step 1 testbed. A general copied instance is not an exact 3-CNF with m ≤ n.
- Step 2 matches the full sum on the instances in the Step 2 section and on `two-block-copies`. It refuses a copied clause that still crosses blocks. Measurement enumeration for the popcount design is (r + 1)^b, exponential in the number of groups, and on the recorded r = 2 scan it does more local work than summing the cube.
- The Step 0/1 Lagrange path is still Θ(2^n · D) scalar evaluations. The Step 4 paired path is about half as many cube passes, with a jet of width 2. The single-point path is one pass with Θ(D²) arithmetic per assignment. The Step 5 degree rows are the same shape of cost on n' and on m'·k·R. None of these is subexponential.
- The Step 1 degree table is the bench set plus one pytest draw, all at R = 1. Step 3 has its own four-row degree table. The Step 4 degree column is the eight rows in that section; on each of them the measured degree equalled the cap. The Step 5 degree column is the six rows above; on each of them the measured degree equalled M'·R.
- Self-reduction deletes fixed variables from clauses and then calls brute force on all 2^n bits, so the integer it sees is not the number of extensions. The CLI only needs a yes/no answer and checks the witness it returns.


## Step 9 — exact reference schedules and rational design

`conditional_balance.py` minimizes the exact normalized conditional expectation from Lemma 4.5. The pair probabilities include incidences sharing a variable or a factor. Tiny exhaustive-completion tests check the expectation at every prefix and the lexicographic minimum chosen at every step. The final potential is also computed directly from label and group loads.

`design_problems_at` defaults to a cached exact rational row span. Membership implies invariance on the complete rational kernel, hence on Boolean inputs. Optional `boolean_check=True` retains the small-instance comparison. `schedule_placement(..., design_method="boolean")` retains the Step 8 budget gate.

The admitted n=1024 unit-clause instance at r=2, Q=2, R=1 completes placement: multiplicity 64 <= 128, largest group 512, exact design and T=0 restoration pass. Its count is not evaluated. The tail probability is 7/16 > 1/4, so this parameter point fails Lemma 4.7. At n=16384, r=2, Q=3, R=2, tau=1/2, every numerical hypothesis passes and the tail is 5/32 <= 1/4. Step 12 below completes sparse placement at that point and records the global counting budget refusal.

## Step 10 — streamed shared-measurement Taylor recovery

`common_count.py` implements Lemmas 3.2–3.3 and equation (3.13), using the safe degree bound dH=M(R-1), q=floor(dH/Q), J=q+1. Independent integer rows are selected from each restricted measurement matrix. Measurement ranges use the sums of negative and positive row coefficients. Representatives are chosen before imposing E's fixed-bit restriction, as required by Lemma 3.2.

Each requested jet matches the independent full-cube Hasse evaluator; H(0) matches brute force. The local scans are deliberately counted, not presented as a speedup:

| example | count | substitution pairs | local vectors actually examined | planned upper bound |
| --- | --- | --- | --- | --- |
| swap-4 | 1 | 1 | 28 | 32 |
| remark-4-4 | 1 | 5 | 260 | 310 |
| cross-r3 | 7 | 1 | 40 | 48 |
| r2-q3 | 7 | 7 | 371 | 504 |

The evaluator streams tuples and representatives rather than storing an exponential table. It does not call brute force or full-cube Hasse evaluation as its backend. Those functions are independent test oracles. Second-order substitution recovery, negative literals, empty inputs, free variables, malformed designs and refusals before enumeration are tested.

## Step 11 — structured auxiliary entries

`structured_auxiliary.py` evaluates a coordinate of V0 fI(W0) by the Fourier eigenvalues on its commuting-swap orbit. Krawtchouk sums require polynomially many operations in Q and R per queried entry. The replacement coefficient is then equation (4.10).

Every matrix and replacement entry agrees with dense side lengths 4, 9, 16, 64 and 81. Integer measurement descriptors agree with the original dense eigenbasis. At r=2, Q=192, R=32, D=2^384, coordinate restoration, support, sample-point cancellation and the interpolation values on I are checked without enumerating coordinates or measurements.

## Step 12 — sparse placement and bounded rational kernels

`sparse_placement.py` uses only occupied auxiliary coordinates and coalesces contributions that belong to the same original variable. It stores sparse polynomial and measurement rows, with predictive limits on placement entries, replacement support, measurement-mode queries, pullback candidates and elimination entries. `SparseRationalRowSpan` checks rational kernel membership without allocating dense group matrices. The checks also cover T=0 restoration, common denominator, absolute coefficient norms and all sample-point designs. `schedule_sparse_placement` requires every numerical hypothesis of Lemma 4.7.

The n=16384 alternating-polarity unit-clause formula at r=2, Q=3, R=2, tau=1/2 passes those checks. Its concrete storage and query record is:

| quantity | value |
| --- | --- |
| label copies | 256 |
| largest group | 8192 |
| group ranks at each of three nodes | (2048, 2048) |
| sparse polynomial entries | 59392 |
| corresponding dense polynomial slots | 536870912 |
| sparse measurement entries | 14336 |
| mode-coordinate queries | 163840 |
| pullback candidates | 34816 |
| elimination entries | 10240 |

`common_count.py` accepts this sparse placement. Its global plan has q=5461 and J=5462, with a Boolean-vector upper bound of 36697 bits. The shared-measurement evaluator refuses before enumeration, returns no count, and examines zero Boolean vectors. Step 14 below completes the same placement's count through disjoint replacement supports. This record establishes finite storage usage, not an asymptotic space or time theorem.

## Step 13 — supplied CNFs and independent component counts

`dimacs.py` preserves the declared variable count and parses comments, multiline clauses and empty clauses with header validation. `component_count.py` removes tautologies and duplicate literals/clauses, finds independent components by union-find, renumbers each component and places it afresh. Every nontrivial local count comes from sparse shared-measurement recovery. Identical renumbered components optionally reuse a completed count in a capped cache. Declared unused variables contribute 2^free. A refusal returns no partial count.

```bash
python -m counting_sat_algebra count examples/count_demo.cnf --output benches/count_demo.json
counting-sat-algebra count INPUT.cnf --json
```

The example declares eight variables and has two disjoint three-literal clauses. The result is 196 = 7·7·2², with two components, two free variables and 148 local Boolean vectors examined. The saved JSON includes the exact count and resource counters. Huge results use exact hex encoding to avoid Python's decimal-digit limit.

On 10923 disjoint three-literal clauses and 32769 active variables, the result is exactly 7^10923. One local evaluation and 10922 memo hits produce this count; adding two declared free variables gives 4·7^10923. Tests disable the full-cube oracles while executing this path. The large input passes the Q=3, R=2 numerical hypotheses, but its global placement was not executed. Its count follows from independent components and does not establish a general connected-instance bound.

The Step 13 milestone passed 163 regression tests and 14 reference checks. Tests include disabled memoization, a local R=2 Taylor backend, normalization, random small formulas, declared free bits, exact DIMACS/JSON output and component/global budget refusals. Step 14 extends the saved reference report with the checks below.

## Step 14 — local Taylor contraction of the same global placement

`support_count.py` forms a dependency graph from the actual replacement-polynomial coefficients in each clause. This includes couplings that were absent from the original CNF. For disjoint support sets V_i, H(T)=2^free·product_i H_i(T). Each local cube is summed at T=1 with a Taylor jet covering the component's full degree bound, then Hermite evaluation recovers H_i(0). No original-CNF counter or global-cube oracle is used as the backend.

The same admitted n=16384, r=2, Q=3, R=2 placement from Step 12 has 10240 support components, no free variables, largest component size four, largest component degree four and total degree bound 14336. Both modes restore exactly 1, independently checked against the unique assignment of alternating-polarity unit clauses:

| quantity | memoization enabled | memoization disabled |
| --- | --- | --- |
| exact count | 1 | 1 |
| completed support components | 10240 | 10240 |
| local Taylor evaluations | 13 | 10240 |
| memo hits | 10227 | 0 |
| actual local Boolean vectors | 96 | 49152 |
| planned vector upper bound | 49152 | 49152 |
| factor evaluations | 288 | 128000 |

The vector plan is sum_i 2^|V_i| and assumes no reuse. A budget of 49152 admits both modes; 49151 refuses before local enumeration even with caching enabled. A connected thirteen-variable replacement support exceeds the default component limit twelve and also refuses with no count and zero scanned vectors.

Independent checks compare the original four small placements with shared-measurement recovery and brute force. A nonconstant coupled replacement is checked coefficient by coefficient against the full-cube Hasse oracle. Tests also cover free variables, empty clauses, invalid designs, disabled memoization and budget checks before contraction.

The Step 14 milestone passed 181 pytest tests and 17 reference checks. `--step14` selects its three rows. Reports retain the exact counts, operation statistics, executed steps and failure status without wall times. Empty or failing suites cannot be marked as successful because earlier suites passed. The incomplete theorem and Lean proof flags remain false. Step 15 below records the next historical checkpoint.

## Step 15 — finite Taylor elimination on connected supports

`tensor_count.py` eliminates Boolean variables of the actual replaced factors over \(\mathbb Q[z]\). Each live table is a Taylor jet at T=1 through that component's full degree bound. Eliminating one Boolean variable preserves the sum-product polynomial. Every component is planned before any numerical table is built. A refused call returns no count. This backend does not memoize. Slot accounting excludes input polynomial storage, planner metadata, arithmetic temporaries and integer bit lengths.

The same 13-variable connected path that Step 14 refuses, because its support exceeds the 12-variable component cap, is counted by Step 15:

| quantity | value | label |
| --- | --- | --- |
| exact count | 610 | `COUNTED` |
| Step 14 status on the same placement | refused | `REFUSED` |
| components | 1 | |
| largest component / max bucket width | 13 / 2 | |
| local Boolean vectors | 98 | |
| factor evaluations | 48 | |
| jet products | 190 | |
| peak table entries | 50 | |

`examples/connected_path.cnf` through `count --backend tensor` returns the same 610. The closed form is the number of length-13 bit strings with no adjacent zeros, Fibonacci(15).

A 1024-variable path of the same shape, plus two declared free variables, counts as 4·Fibonacci(1026). The plan has one component of 1024 variables, bucket width 2, 8186 local Boolean vectors and a peak of 4094 table entries. `count_cnf_components(..., backend='tensor', max_component_variables=1024)` returns the same count and the same vector total. The long-path test disables the global-cube, shared-measurement, component-count and Step 14 jet oracles while this contraction runs.

Two coupled replacement factors on two variables recover count 1. The jet at T=1 is 7/4, 1/2, −1/4, matching the independent full-cube Hasse jet coefficient by coefficient. That coefficient comparison is `CHECKED`; the count is `COUNTED`. A single width-13 factor is `REFUSED` before any numerical table: the factor support bound 13 exceeds the budget 12, and the Boolean-vector counter stays at zero. An induced elimination bucket wider than 12, and a later oversized component, also refuse before numerical work.

DIMACS accepts `--backend common|tensor`. `common` remains the shared-measurement default. `tensor` places each independent CNF component afresh and then calls this contraction, with `--max-tensor-width` defaulting to 12 and `--max-tensor-degree` defaulting to 64.

All 205 pytest tests and all 21 reference checks pass at version 0.16.0. `--reference` writes Steps 9–15 to `benches/reference_alignment.json`. `--step15` selects the four rows also saved in `artifacts/step15-verification.json`. Both reports keep `theorem_1_1_complete=false` and `lean_proof_complete=false`.

## Step 16 — Lemma 4.8's schedule and one finite admitted member

Shipped in `reference_parameters.py`. `lemma_4_8_proposal` follows the proof's order. λ = ⌈log₂ n⌉, ν = ⌈log₂ λ⌉, r is the greatest integer with r²ν⁴ ≤ λ, Q = ⌊λ/(128ν)⌋, and R = (β+3)⌈√λ⌉, where β ≥ 1 is the least integer with b, w ≤ λ^β. τ = 1/(rλ), ε = 1/(br²λ), D = r^{2Q} and the centered interval are formed only after r ≥ 2 and Q > R. A failure of that guard does not form D and does not allocate a matrix.

At 16384 alternating unit clauses, λ=14, ν=4, β=1, so the closed form is r=0, Q=0, R=16. D is not formed. The same guard fails at λ=2¹⁹, where the integers are ν=19, r=2, Q=215, R=2900 and Q ≤ R. Neither check is an executed counter. Label `CHECKED`.

`choose_asymptotic_parameters` then searches r, Q and R inside explicit budgets and calls `check_reference_parameters` before accepting a tuple. The interval is the centered one from the proof of Lemma 4.8. τ is the least value in (0, 1/2] that meets the exact binomial tail; equality is allowed. Among accepted placeable tuples, the order is smaller τ, then smaller R/Q, then smaller 1/r. That comparison does not prove the O() bound on Lemma 4.8's exponent.

| instance | what the chooser reports | label |
| --- | --- | --- |
| 8192 alternating units | no tuple. Tail 7/16 at r=2, Q=2, R=1 exceeds every τ ≤ 1/2. r=2, Q=3 needs 2(b+w)D² = 16384 > M | `CHECKED` |
| 16384 alternating units | r=2, Q=3, R=2, τ=5/16, interval {0, 1}, tail 5/32 = τ/(rb). D=64. Not R ≥ 32 and not 6R ≤ Q | `CHECKED` |
| same n, trial budget 1 | the same tuple is admitted and then refused. Alphabet 8 and 262144 candidate trials are not started | `REFUSED` |

The hand-picked point at this n used τ=1/2. The tail is unchanged. The measurement allowance 2τn falls from 16384 to 10240. Sparse placement at the chosen tuple reproduces the Step 12 record, because the schedules do not read τ:

| quantity | value |
| --- | --- |
| label copies | 256 |
| largest group | 8192 |
| group ranks at each of three nodes | (2048, 2048) |
| independent measurements | 4096, within the allowance 10240 |
| sparse polynomial entries | 59392 |
| dense polynomial slots | 536870912 |
| sparse measurement entries | 14336 |
| mode-coordinate queries | 163840 |
| pullback candidates | 34816 |

Counting uses the existing backends. A refusal returns no count and does not start enumeration.

| backend | result | label | operation record |
| --- | --- | --- | --- |
| common | no count | `REFUSED` | q=5461, J=5462, Boolean-vector bound of 36697 bits, 0 vectors scanned |
| support, memoization on | count 1 | `COUNTED` | 10240 components, largest support 4, degree at most 4, 13 evaluations, 10227 memo hits, 96 local Boolean vectors, plan 49152 |
| tensor | count 1 | `COUNTED` | same 10240 components, bucket width 4, 176128 local Boolean vectors, 98304 factor evaluations, 262144 jet products, peak table size 72, 360 peak coefficients, 268288 planning steps |

The count 1 matches the single assignment of alternating-polarity unit clauses. The support count and the tensor count are disjoint-support evaluations of this one placement. They are not a general connected-instance bound. The tensor backend is still Step 15's finite Taylor elimination, not Appendix D. `build_auxiliary_slice(2, 192, 32)` still refuses before allocation; the reported side length remains 2³⁸⁴.

All 214 pytest tests and all 31 reference checks pass at version 0.17.0. `--reference` writes Steps 9–16 to `benches/reference_alignment.json`. `--step16` selects these rows, also saved in `artifacts/step16-verification.json`. Both reports keep `theorem_1_1_complete=false` and `lean_proof_complete=false`. No wall times are recorded.

## Step 17 — sparse systems, Lemma 5.1, and the truncated product

Shipped in `sparse_systems.py`. A polynomial is a list of variable-index sets. Duplicate monomials cancel modulo 2, a repeated index inside one monomial is stored once, and the empty monomial is the constant 1. The supplied list must have degree at most \(k\) and length at most \((\log_2 n)^c\), with that inequality decided by exact integer arithmetic. The integer lift sums the monomial products. On a Boolean point it is even exactly when the \(\mathbb F_2\) equation is zero.

Lemma 5.1 builds
\[
p_\ell(Y)=(1-Y)^\ell\sum_{j=0}^{\ell-1}\binom{\ell+j-1}{j}Y^j
\]
and divides by the monic \(h(Y)=Y^2-Y\). The construction uses only integers. For \(\ell=2\), \(p_2(Y)=1-3Y^2+2Y^3=(1-Y)+(2Y-1)h(Y)\). For \(\ell=1,2,3,4\) and every integer \(y\) from \(-6\) through \(12\), \(p_\ell(y)\equiv 1\{2\mid y\}\pmod{2^\ell}\). Label `CHECKED`.

The product keeps a running count of \(h\) factors and drops any term whose count reaches \(\ell\). At \(\ell=2\) the retained powers are \((0,0)\), \((0,1)\) and \((1,0)\). The omitted term is divisible by 4 at integer inputs. Each retained term has at most \(m_{test}\) linear factors and at most \(\ell-1\) quadratic factors. Label `CHECKED`.

| system | result | label | record |
| --- | --- | --- | --- |
| \(n=4\), \(x_0x_1+x_2=0\), \(x_3+1=0\), \(\ell=5\) | count 4 | `COUNTED` | 16 Boolean points, 15 retained terms, coefficient width 7 bits, input-degree bound 20 |
| empty equations, \(n=5\) | count 32 | `COUNTED` | \(2^5\), no scan and no term table |
| constant equation \(1=0\), \(n=4\) | count 0 | `COUNTED` | no scan |
| same \(n=4\) system, term budget 2 | no count | `REFUSED` | term table not allocated, coefficients not constructed, 0 vectors |
| same system, coefficient budget 1 bit | no count | `REFUSED` | construction stops, term table not allocated, 0 vectors |
| same system, cube budget 1 | no count | `REFUSED` | no scan |

The count 4 is the least nonnegative residue of the truncated indicator modulo \(2^5\). It matches an independent \(\mathbb F_2\) scan, and the residue equals the common-zero indicator at every one of the 16 points. Sections 5.3–5.5 are not run: there are no replacement rows and no sum in the Section 3 algebra \(A_\ell\). Section 6 and Appendix D are not run. The cube scan is not Theorem 2.2's bit-complexity bound.

All 227 pytest tests and all 40 reference checks pass at version 0.18.0. `--reference` writes Steps 9–17 to `benches/reference_alignment.json`. `--step17` selects these nine rows, also saved in `artifacts/step17-verification.json`. Both reports keep `theorem_1_1_complete=false` and `lean_proof_complete=false`. No wall times are recorded. Step 17 does not implement Theorem 1.1.

## Step 18 — replacement forms and the sum in \(A_\ell\)

Shipped in `algebra_count.py`. Each non-constant equation supplies the factor
\[
B_i=\sum_{j=0}^{\ell-1}e_j(L_i)h(L_i)^j V^j\in A_\ell,\qquad A_\ell=\mathbb Q[V]/(V^\ell),
\]
with \(e_j\) from Lemma 5.1 and \(h(Y)=Y^2-Y\). A neutral factor equal to 1, with scope \(\{v\}\), is added for every declared variable. The product degree bounds are the ones in Lemma 5.2 and (5.5): \(d=k(m_{test}+2(\ell-1))\) and \(d_H=kR(m_{test}+2(\ell-1))\). For the \(n=4\), \(k=2\), two-equation system these are \(d=d_H=20\). Label `CHECKED`.

The sum is Proposition 3.4. `separated_group_sum` is the same Lemma 3.2 scan the CNF counter uses; the coefficient ring is \(A_\ell[z]/(z^J)\). Each coefficient of \(V\) is interpolated separately. The rows are a finite design that restores \(\delta_v\) at \(T=0\): one block holds every variable, the other is empty, and \(P_{f,v}(T)=(1+T/Q)e_v\) with \(Q=d_H+1=21\) and \(R=1\). This meets Definition 3.1 for these scopes. It is not the Section 4 auxiliary array and it is not a Lemma 4.8 row set. On this system Lemma 4.8 returns \(r=1\), \(Q=0\), \(R=10\), so \(\tau\) and \(D\) are not formed.

| system | result | label | record |
| --- | --- | --- | --- |
| \(n=4\), \(x_0x_1+x_2=0\), \(x_3+1=0\), \(\ell=5\) | count 4 | `COUNTED` | \(V\)-coefficients \((-2,6,0,-192,2496)\), unreduced sum 2308, residue 4 modulo \(2^5\); 399 Boolean vectors, bound 714, 2016 factor evaluations; \(H(1)\) is not 4 |
| same system after the Section 5.5 split is unnecessary | heavy set empty | `CHECKED` | \(M=8\), \(b_0=4\), four neutral scopes |
| five copies of \(x_0=0\) and \(x_0+x_1=0\), \(n=4\) | count 4 | `COUNTED` | heavy variable \(x_0\), \(b_0=6\), branch counts \((4,0)\), 110 Boolean vectors |
| that branch, \(n'=3\) | guard fails | `CHECKED` | Lemma 4.8 gives \(r=1\), \(Q=0\), \(R=8\); \(D\) is not formed |
| vector budget 1 on the two-equation system | no count | `REFUSED` | bound 714, 0 vectors, jets not allocated, coefficients not constructed |
| branch budget 1 on the six-equation system | no count | `REFUSED` | 2 branches refused before the branch loop |
| coefficient budget 1 bit | no count | `REFUSED` | parity coefficients not constructed, jets not allocated |

The residue 4 matches both an independent \(\mathbb F_2\) scan and the Boolean cube sum of the truncated factors in \(\mathbb Z[V]/(V^5)\). A length-2 jet at \(Q=3\) and a two-block unit design each recover the count 2 of \(x_0=0\) on two variables. Those checks are in pytest. Section 6 and Appendix D are not run. The finite sum is not Theorem 2.2's bit-complexity bound.

All 237 pytest tests and all 47 reference checks pass at version 0.19.0. `--reference` writes Steps 9–18 to `benches/reference_alignment.json`. `--step18` selects these seven rows, also saved in `artifacts/step18-verification.json`. Both reports keep `theorem_1_1_complete=false` and `lean_proof_complete=false`. No wall times are recorded. Step 18 does not implement Theorem 1.1.

## Step 19 — equation selection and disjoint sparsification

Shipped in `equation_selection.py`. A clause violation is a square-free polynomial over \(\mathbb F_2\): a positive literal contributes the factor \(1+x\), a negative literal contributes \(x\), a tautology is 0, and the empty clause is the constant 1. A binary row is a list of \((e,c)\) pairs. Coefficients of 0 contribute nothing, and a repeated clause index cancels modulo 2. Label `CHECKED`.

Lemma 6.1's integers are \(\mu=\lfloor\log_2 m\rfloor+1\), \(P=20000 B^2\mu^2\), and \(D_0=320000 B^2\mu^3\). On \(n=m=3\), \(B=1\) these are \(\mu=2\), \(P=80000\), \(D_0=2560000\). Here \(q_G=1\) and \(h\le 0\), so \(A(0)\) is the two even-parity checks. Every nonzero kernel vector has weight 2, which is at least \(m/P\). There are 2 rows, each of weight 3, and \(2\le 3n\). Label `CHECKED`. The local-code search at \(n=m=2048\), \(B=1\) is 524288 bits and is refused before a matrix is allocated. Label `REFUSED`. No kernel vector is scanned.

Lemma 6.2 appends rows only when \((2m)^P\) fits the candidate budget. Ties keep the first candidate in enumeration order. On two unit clauses in two variables, list length 1 uses 6 rounds and 24 counts. The final \(N(A)\) is 1, equal to \(\#F\), and `count_common_zeros` on the same system returns 1. Label `COUNTED`. The paper-scale \(P=80000\), a candidate budget of 1, and a cube budget of 1 each return no count and enumerate nothing. Label `REFUSED`.

Lemma 6.3 keeps the declared \(n\). Normalization deletes repeated literals, tautologies, duplicate clauses, and clauses that contain another clause. It does not propagate units and it does not copy variables. Children of a branch partition the parent's satisfying assignments. On the 3-variable 2-CNF the leaf counts are \((3,1)\), summing to 4, with occurrence at most \(B=6\). Label `COUNTED`. A node budget of 1 returns no leaves and no count. Label `REFUSED`. Step 5 implication-cycle copying increases \(n\); this tree does not.

Proposition 6.4 at \(n=20\) gives \(\lambda_n=5\), \(\nu_n=3\), \(s=0\), so the schedule is the exhaustive-search case. The integer witness \(s=2\) is only \(\lambda=521284\), \(\nu=19\). For \(k=2\) and \(\sigma=1\) the Appendix A power of two is \(M=64\). Label `CHECKED`. The \(O_k(\rho(n))\) bound is not claimed.

| system | result | label | record |
| --- | --- | --- | --- |
| \(n=m=3\), \(B=1\), even-parity \(A(0)\) | min weight 2 | `CHECKED` | \(\mu=2\), \(P=80000\), \(D_0=2560000\), 2 rows of weight 3, 8 vectors |
| \(n=m=2048\), local-code search | no matrix | `REFUSED` | 524288 bits, 0 vectors, no allocation |
| two units, \(n=2\), list length 1 | \(N(A)=1\) | `COUNTED` | 6 rounds, 24 counts, 96 cube vectors; `count_common_zeros` returns 1 |
| paper-scale \(P=80000\), candidate budget 1, cube budget 1 | no count | `REFUSED` | nothing enumerated |
| 3-variable 2-CNF, \(M=2\) | leaf counts \((3,1)\) | `COUNTED` | sum 4, same \(n\), occurrence at most 6 |
| tree node budget 1 | no count | `REFUSED` | no leaves, 0 vectors |
| \(n=20\) schedule | \(s=0\) | `CHECKED` | \(\lambda_n=5\), \(\nu_n=3\), exhaustive case; witness \(s=2\) only at \(\lambda=521284\) |

Appendix D is not run. These finite checks are not a proof of Proposition 6.4's running time, and they are not Theorem 1.1.

All 249 pytest tests and all 56 reference checks pass at version 0.20.0. `--reference` writes Steps 9–19 to `benches/reference_alignment.json`. `--step19` selects these nine rows, also saved in `artifacts/step19-verification.json`. Both reports keep `theorem_1_1_complete=false` and `lean_proof_complete=false`. No wall times are recorded. Step 19 does not implement Theorem 1.1.

## Step 20 — Appendix D finite modular networks and local recovery

`modular_ring.py` implements the base ring B_{p,ell} and its unramified extension A_p. Coordinates are integer residues, in V-major order. The extension degree is e=4 ceil(log2(s)) when selected by `extension_ring_for_size`; s=2 chooses f=U^4+U+1. Multiplication reduces by the integer polynomial f, and discards powers V^ell and above. Units are inverted through the residue field, Newton lifting and truncated V-series inversion. Binary-coordinate residue lifts provide points with unit differences for interpolation.

`modular_network.py` implements ordered ports in a bipartite multigraph, including distinct parallel edges and isolated scalar vertices. Each implicit tensor carries its creation ring. Entries are evaluated at that precision and canonically lifted separately; the lifting operation is not a ring homomorphism. The direct leaf contraction streams every Boolean edge assignment and refuses an oversized scan before evaluating entries.

D.2's CNF reduction uses trivalent equality trees followed by a subdivision of every edge. Repeated literals and tautologies are normalized, empty clauses give zero, and declared free variables contribute their power of two. The example has 12 edges, 14 vertices and maximum degree 3. Its residue modulo 2^9 is **196**, with 4096 edge assignments and 57344 tensor evaluations. Four identical positive units on three declared variables give **4** through a 10-edge, 11-vertex equality-tree network. Label `COUNTED`.

`modular_recovery.py` takes an explicit finite packet layout. It cancels known powers of two in W_xy(t+alpha z), and inverts only odd denominators. For a=1, t=1, x=0, y=1, the exact local coefficients are (-3,-22,-48,-32). W(0)=I is checked. The layout does not supply symbol groups and does not establish Lemma D.3.

The D.8 evaluator streams independent row/column port pairs. At each active node, interpolation at degree_bound+1 residue-field points recovers the whole local polynomial; D.5 keeps the specified prefix. D.9 uses modular quotient coefficients and rational final weights, then reduces at the root precision. Every node and its precision, coefficient slots and pair scans are planned before numerical arrays or entry calls.

| instance | result | record |
| --- | --- | --- |
| two-edge network, p=5, ell=2, a=u=1 | 9+4V, `CHECKED` | nodes -1/1; multiplicity 6; precision 133; 224 port pairs, 448 entries |
| same network, a=u=2, two unused packet slots | 9+4V, `CHECKED` | nodes -3,-2,-1,1,2,3; multiplicities 6,3,6,6,3,6; precisions 149,101,149,149,101,149; 1440 port pairs, 2880 entries |
| nonlinear row, degree-6 deformation | 8+5V, `CHECKED` | independent rational polynomial oracle; the leading coefficients are 11/4 and 155/4 |
| network vector budget 4095 on the example | `REFUSED` | no network allocated; no count; zero scans |
| recovery precision 132, pair budget 223, or jet-slot budget 751 | `REFUSED` | no jets or kernels; no value; zero scans and entry calls |
| base residue field of size 2 with 7 required points | `REFUSED` | interpolation cannot start; no value or scan |

The coefficients above are modulo 32 and V^2=0. At T=1 the first network has coefficients (26,13), so that value is not its contraction. Creation-precision tests lift entries -1 and 6 from precision 2 to representatives 3 and 2 at precision 8; two edge terms sum to 12.

Step 20 has fourteen reference checks saved by `verify --step20 --output artifacts/step20-verification.json`. The independent pytest checks include irreducibility by trial division, integer extension arithmetic, all small Cauchy coefficients against rational products, complete-degree polynomial recovery, multiple valuation levels, CNF normalization, tree degrees, fixed creation precision, and refusal before entry calls.

All 289 pytest tests and all 70 reference checks pass at version 0.21.0. `verify --quick` and both CI report validators pass. The full report records Steps 9–20; Step 20 has 40 pytest tests and 14 reference rows. The theorem and Lean completion flags remain false. No wall times are recorded in the reference reports.

At the Step 20 checkpoint, D.4.2's collision packing, D.4.5–D.4.7's modes and recursive children, and D.3's incidence encoding remained open. Step 21 below supplies finite packing and child port topology. Direct leaf enumeration and D.8 port-pair enumeration are exponential. The resource counters and conservative recovery bounds are not Theorem D.1's time/space proof. Step 20 does not implement Theorem 1.1. See `APPENDIX_D.md` for API examples and exact scope.

## Step 21 — collision labels, packet layouts and child topology

`modular_packing.py` implements D.4.2's integer schedule and deterministic finite collision packing. D.4 uses r=64 and epsilon=1/131072. H and Delta use exact integer roots; u and cumulative a use a popcount identity, without a loop through large depths. On s=16, M=20, lambda=4, H=8, Delta=13 and u=a=0, so the active guard stops before L or D_arr is formed. The logarithmic witness lambda=131072 gives u=1, a=2, previous_a=0; lambda=262144 at depth 1 gives u=1, a=6, previous_a=4. Neither witness constructs the corresponding s.

The first D.6 condition is **M >= 2*D_arr^2*Delta**, with D_arr squared. The other conditions are M >= L*D_arr and 2*L <= Delta. `pack_network` requires all three by default. Explicit `require_hypotheses=False` checks the actual finite conclusions after packing; it does not change the numerical flags or establish membership in the fixed D.4 schedule.

Exact conditional expectations use completed pair loads and assigned-neighbor incidence histograms. Shared endpoints and parallel edges are correlated. If K=r^a and P is the sum of squared parallel multiplicities, then K^2 E[Phi]=(K^2-1)P+(K-1)(sum d_v^2-2P). Vertex order is preserved and ties take the lowest label. Pytest independently enumerates every remaining vertex label to check each candidate expectation, including interleaved sides and parallel edges. The shared/parallel six-edge example has initial potential 23/2, final potential 1, labels (0,1,0,1), and pair loads ((2,1),(1,2)).

Every pair-label class is padded to the same H_pad, a multiple of L. Slot bits select orientations for unequal symbols and distinct copies for equal symbols. Empty slots have no ports. Certificates recompute the loads, potentials, D.6 flags, slot orientations, packet counts and real group port bounds. Label trials, score updates, optional trace and conservative padded storage are all budgeted before label search. A finite conclusion failure returns no packing but reports completed work accurately.

| instance | result | record |
| --- | --- | --- |
| 128-edge matching, r=2, a=1, Delta=4 | `CHECKED` | D.6 admitted; H_pad=32, 64 packets, potential 96 to 0, 64 ports per side/group, largest child topology 128 edges and degree 4; no M/8 shrink |
| 32-edge matching, r=2, a=2 | `CHECKED` | explicit finite mode; H_pad=4, 16 packets, orbits with 0, 1 and 2 equal coordinates; D.6 not admitted |
| 8192-edge matching, r=64, a=1, Delta=4, default mode | `REFUSED` | first D.6 threshold is 134217728; no labels or packets allocated, zero trials |
| same 8192-edge matching, explicit finite mode | `CHECKED` | potential 8190 to 0; all 4096 pair classes have load 2; H_pad=2, 4096 packets, 128 ports per side/group; largest child topology has 256 edges and degree 4, within M/8=1024 |
| two repeated positive units on two declared variables | count 2, `COUNTED` | automatic packet layout plus Step 20 recovery; 4 network edges, residue 1 and one free variable; 6656 port pairs, 33280 tensor evaluations |
| label, trace, padding or score budget too small | `REFUSED` | no packing, labels or packet allocation; zero trials and entry evaluations |
| 8192-edge matching, child edge budget 255 | `REFUSED` | existing packing retained; 256 required child edges; no child topology allocated |

The r=2 and r=64 matching tuples are not D.4 schedule members. The latter also fails D.6 despite satisfying the actual finite conclusions. Its label search uses 1048576 trials and 1572864 score updates. Packing and child topology perform zero tensor evaluations.

`child_topology` implements only D.4.6 port geometry: selected source vertices on one side, packet vertices on the other, and one edge for each selected original port. It preserves source port order, isolated original vertices and zero-degree packet vertices. Edge and degree bounds are checked on the resulting finite topology; no packet tensor, mode coefficient or recursive contraction is evaluated. The automatic CNF recovery still uses Step 20's port-pair enumeration.

All **322 pytest tests** and **84 reference checks** pass at version **0.22.0**. `verify --quick` and both CI report validators pass. The full report records Steps 9–21. Step 21 has 33 pytest tests and 14 reference rows; `verify --step21 --output artifacts/step21-verification.json` saves those rows separately. The theorem, Lean and recursive-contraction flags remain false. No wall times are recorded in reference reports. API examples, resource preflight formulas and scope are in [PACKET_PACKING.md](PACKET_PACKING.md).

## Step 22 — finite weighted modes

Version 0.23.0 checks one packet's D.4.5 mode table. A crossing slot at \(a=2\), \(t=1\), precision 4 recovers Cauchy residue 8. The two-edge network still recovers \(9+4V\) through Step 20 after the mode budget passes, with 224 port pairs and 448 entry evaluations. Mode refusal returns no count and does not call recovery. Packet tensors, recursive contraction, and the theorem flags stay untouched. See [WEIGHTED_MODES.md](WEIGHTED_MODES.md).

## Step 23 — Lean syntax and a streamed disjoint frontier

`lean_encoding.py` preserves the declared n, Lean's literal polarity and framed binary lists; the four-variable star encodes in 138 bits. Count 196 has the little-endian word `(0,0,1,0,0,0,1,1)`. `lean_parameters.py` computes the first passing dyadic choice with exact ceiling logs: k=3 and sigma=1 choose M=512, exponent 9; every smaller exponent fails.

`lean_sparsify.py` normalizes clauses, selects the priority core and branches over all clauses of that width containing it. The original variable set stays fixed. A structural preflight is replayed for counting; only the depth-first frontier is retained. The four-variable star at M=2 has two leaves and count **9**. Three independent stars have eight leaves, depth 3, at most four pending nodes and count **125**. Completed leaves and their histories are not stored. Structural refusals precede leaf evaluation; a later leaf refusal reports completed work and no partial total.

Fourteen reference rows are exposed by `verify --step23`. These are finite syntax, semantics, parameter and resource checks, rather than the complete Lean leaf algorithm or its Python asymptotic transfer.

## Step 24 — constructed ternary queries and native leaves

`lean_ring.py` constructs Newton coefficients modulo `3^t` and divides the normalized polynomial into t three-coefficient digits. The dense ring has independent U and Z cutoffs and reads out U only after the product. `lean_squarefree.py` implements support-union multiplication and Boolean monomial weights. `lean_query.py` counts history acceptance or acceptance intersected with a marked original clause violation; it preserves raw clause IDs, duplicates, original n and fixed ambient coordinates.

For three positive units on five declared variables, one all-one history row gives accepted count **8** and marked count **4**. At t=4 the accepted marker sums are `(8,36,45,0)` modulo 81, so higher coefficients participate. Unit history rows count the original CNF. `--backend lean --leaf-backend native` returns **196** on the supplied example and **9** on the four-variable star; the star uses 51 monomial products and 1741 scalar operations across two leaves. The nine-variable three-star formula returns **125**. All these native runs scan zero Boolean assignments.

Every native query accounts for its scalar and monomial work; the entire frontier shares those budgets. A 1000-scalar-operation limit finishes the star's first leaf and refuses the second with `count=None`. Live immutable tables, scratch and output coefficients are charged in logical encoded bits; they are not Python heap or visited-tape measurements. The complete squarefree support table may still grow exponentially.

Seventeen reference rows are exposed by `verify --step24`. Independent unit tests check small-precision Hermite/Hasse conditions, dense convolution, raw accepted/marked semantics, fixed coordinates and global budget refusal. The source mapping and API are in [LEAN_NATIVE.md](LEAN_NATIVE.md). The separate static Lean audit covers 3969 modules without missing local imports, code-level placeholders or local axiom declarations; `kernel_verified=false`.

Version 0.30.0 passes **582 pytest tests** and **194/194 reference checks** across Steps 9–30. The full report and every Step 15–30 artifact were regenerated by executing their suites; each records version 0.30.0 and false completion flags.

## Step 25 — scheduled sparse codes and violation-mass greedy selection

`lean_code.py` builds the scheduled ternary Tanner code from the exact `SparseCodeParameters` integers (floor-log `logSize`, `gapDenom = 200*B*logSize`, `degreeBound = 8*gapDenom^2*logSize`, `degree = min(m, degreeBound)`, `localRows = degree/(8B)+1`, local distance `degree/(100*B*logSize)+1`, `distanceDenom = 20000*B^2*logSize^2`, least cube dimension `r` with `degree*2^r >= m`). Both promised searches are complete increasing odometer enumerations: Walsh shifts meeting `q*|sum +-1| <= D`, then a local matrix whose every nonzero codeword has weight at least the distance. The saved witness is the first accepted candidate. Independent checks cover the exact rational identity `(1/(100 B l))*((1/(100 B l)) - 1/gapDenom) = 1/distanceDenom`, the kernel lower bound `distanceDenom*minimum_weight >= m`, the depth bound, and hand-built accept/reject certificates. Instances (4,1) through (256,2) are `checked`; m=16384 and 2^20 refuse before allocation on the matrix-size budget, and candidate or word budget refusals never keep a partial code.

Rows regenerate from slots in the reverse `serializedRows` order: row i has chronological index `row_count - 1 - i`, column `(side ? vertex ^ shifts[port] : vertex)*degree + port`, kept only below m, coefficient from the local matrix. No candidate library is stored; the greedy re-enumerates the canonical ternary family rank by rank.

`lean_greedy.py` starts from the code rows and scores each candidate with real marked queries accumulating every clause's violation count: the mass of `candidate :: history`. A strict minimum keeps the earliest candidate, and the fixed round count is `loopCap = 2*(n + size m + 1)`. Once the mass is zero, further rounds only prepend zero-mass rows and cannot shrink the accepted set, so the port stops at the first zero-mass history; the final accepted query then runs separately. The star's full trajectory is reproduced by brute force over assignments (mass 12 to 4 to 1 to 0 in 3 rounds, 84 marked queries, count 9); the six-variable double star counts 25 and the supplied example 196. `--backend lean --leaf-backend hash` runs the pipeline on every streamed leaf (star: 2 leaves, 200 marked queries); marked-query budgets are shared across the frontier, and a 50-query cap finishes the first leaf but refuses the second with no partial total.

Fifteen reference rows are exposed by `verify --step25`. The Lean schedule's candidate depth is `distanceDenom` and lies above the Python candidate budgets, so the depth is an explicit finite parameter. These are finite schedule, regeneration and scoring checks; they are not Lean primitive-machine step or tape bounds and do not establish a general-input asymptotic claim.

## Step 26 — high-variable frontend and the streamed rank loop

`lean_high.py` computes `EarlyThreshold`'s paid integers without floating point: `Asymptotic.L n = floor(log2 (n + 16)) + 1`, threshold `isqrt(L) + 1`, cutoff `ceil(totalScope * threshold / ambient)`, and the high set of coordinates whose degree exceeds the cutoff, with an ambient-order mask. Independent checks cover both integer bounds (`cardinality*(cutoff+1) <= totalScope`, `cardinality*threshold <= ambient`) and the equal-degree boundary exclusion. The hub instance (8 variables, 7 clauses) has cutoff 6, one high coordinate and two ranks; the star has none.

The rank loop walks all `2^|H|` assignments of the sorted high coordinates with the binary odometer. Every rank runs one real conditioned query with the high bits fixed and every remaining coordinate retained, and the exact totals accumulate. The rank-sum identities are checked against direct queries for unit rows (hub 129, equal to the brute model count), empty histories (hub 256, star 16), explicit histories, marked clauses, base-fixed coordinates and a two-high-coordinate instance with four ranks. `--high-loop` routes the hash leaf's marked and accepted queries through the loop; the hub, star, forest and streamed examples keep identical counts, the routed hub greedy records two ranks per loop query, and rank, query and algebra budgets refuse without a partial count.

Twelve reference rows are exposed by `verify --step26`. The full parameter frontend, degree and occurrence guards, and the counting extension remain open; none of these finite checks is a Lean primitive-machine bound.

## Step 28 — streamed correction pairs and Taylor readback

`lean_taylor.py` streams the correction pairs `(E, U)` with `|E| < J` and `U ⊆ E` in the exact `DifferencePairs.fold` order: coordinate by coordinate, none first, then out and in, pruning the moment the strict budget is exhausted, so no completed-pair list is ever retained. The counts are `pairCount n J = sum_{j<J} 2^j * C(n,j)` (0 at J = 0, 1 at n = 0), and the Boolean finite-difference identity rewrites every squarefree sum as `sum_{(E,U)} (-1)^(|E|-|U|) * 2^(n-|E|) * p(U)`: the alternating sums reproduce the support-table coefficients exactly, and with `J` above the maximum support size the streamed total equals both the direct weight sum and the brute-force point sum (a signed probe gives 116; the model indicators give the star 9, the double star 25 and the hub 129).

The readback streams per node the `J` exact Taylor coefficients of the node-shifted target as `d` little-endian 2t-bit words: `Q*J*d` flat words including zero words, with a linear prefix charge. `decode_stream` rebuilds the blocks from the stream alone, and Newton-Hermite with `Q*J` above the degree returns the target's value at zero, the unique count word below `3^t`. Pair, assignment, term, jet-order and stream-word budgets refuse without a partial sum.

Nine reference rows are exposed by `verify --step28`. The production native leaf still evaluates the full squarefree support table; switching it to this streamed path, the native placement and the whole-machine resource transfer remain open.

## Step 29 — streamed pointwise correction-pair leaf

`lean_stream.py` evaluates the native ternary queries through the pointwise pair identity: for the marker-sum polynomial `p` of the query, `sum_x p(x) = sum_{(E,U), |E|<J} (-1)^(|E|-|U|) * 2^(free-|E|) * p(U)`. The order bound is `J = (sum of row degrees) * (3*(t-1) + 2) + mark_degree + 1`, where a row degree is the largest post-substitution distinct-variable count and a fully forced-violated clause (fixed literal false) contributes degree zero as a constant-zero polynomial. Every pair performs exactly one point evaluation built from the scalar row values; no squarefree support table is materialized, and `polynomial_terms` stays zero. The precision rule mirrors the native path (`3^t > 2^n` on the declared n), so marker sums agree coefficient by coefficient: the fixture gives accepted 8 and marked 4 with marker sums `(8, 36, 45, 0)`, and the unit-row leaf counts the star 9, the double star 25, the hub 129 and the supplied example 196, equal to the native leaf and to brute force. A star instance whose native run refuses on the polynomial-term budget still counts 9 here.

`--leaf-backend stream` streams the frontier with shared pair and point-evaluation budgets (`--max-stream-pairs`, `--max-stream-evaluations`); a refusal after earlier leaves keeps the completed work but publishes no total. Randomised cross-checks against the native algebra (seeded CNFs with random rows, marks and fixed coordinates) match counts, marker sums and precision exactly.

Nine reference rows are exposed by `verify --step29`. The full machine's degree reduction, the native placement and the whole-machine resource refinement remain open; the streamed cost scales with the pair count, so the default leaf remains `common`.

## Step 30 — leaf assembly and placement

`lean_assembly.py` places one leaf's state at the recorded `QueryReady` port map: modulus bits (length at most `2t+1`) at 194, the `d^3` structure words of `X^(i+j) mod F` at 201, the capacity width word at 202, scalar codes of zero, one, minus one and one half at 203/207/211/253, the unary degree, precision, jet and `3^d-2` exponent words at 204/212/243/209, the `3t`-coefficient idempotent code (`6t^2` bits, identical to the constructed indicator encoding) at 281, flat and framed node tables at 256/300, the interpolation caches (node inverses, negatives, pairwise unit differences) at 713/714/715 and the unary widths at 716/717. Every placed length matches the recorded formula; the structure words reproduce the power oracle, the caches invert, negate and invert differences independently, and the framed node table decodes back to the flat words.

The row stream is framed exactly like `Encoding.encodeList`, decodes back to the raw payloads and to the regenerated slot lists, satisfies the `HashGreedyCode.initial_encoding` identity, and the append encoding is bit-identical to `true :: frame (payload.reverse) ++ history` with the newest row first. `run_assembled_leaf` threads one run from the scheduled code through placement, the framed stream, the decoded history and the greedy rounds to the little-endian count word, verifying the decoded-history marked score against the direct rows; the star, hub, forest and supplied example give 9/129/25/196, equal to brute force and to the direct greedy. Budget refusals (nodes, degree, exponent unary, table, stream) publish no count.

Nine reference rows are exposed by `verify --step30`. The physical 896-stack machine, building the ring target from the query heads, grouping replacement and the whole-program prefix-space refinement remain open.

## Current distance to the reference conclusion

The next dependencies are the counting extension ring and real placement, then streamed correction pairs, node jets and Hermite reconstruction. The current native CNF leaf is a J=1 full-support algebraic sum. It supplies an exact finite query, but its support table must be replaced to obtain the formal program's resource behavior.

Whole-program bit complexity, every-prefix working space, Lean kernel verification and Python refinement remain open. The new library's theorem source is distinct from the old single-file scaffold; the pending work does not mean its final theorem still contains `sorry`. Appendix D packet coefficients, recursive children and incidence encoding remain a separate later route. The recorded finite placement, component, support, Taylor and modular recovery results do not establish general asymptotic bounds. See [ROADMAP.md](../ROADMAP.md) for the dependency order and [REFERENCE_ALIGNMENT.md](REFERENCE_ALIGNMENT.md) for the distinct exponents.
