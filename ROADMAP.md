# Roadmap

Build the pieces of the counting algorithm on instances small enough to compare with brute force. A step is done when that comparison passes and the cost matches the step's own analysis. Refuting ETH is not a success criterion for this repository.

## Done

### Step 0 — counting polynomial and interpolation

- [x] Clause polynomial Φ_e and brute-force #F
- [x] Product identity on disjoint variable groups
- [x] Random replacements of degree ≤ R with value X_v at T = 0; Lagrange recovery of H(0)
- [x] Self-reduction demo that searches for one satisfying assignment
- Cost: each H(t), t ≠ 0, still enumerates 2^n assignments

### Step 1 — bounded-occurrence testbed

- [x] Generator and checker for 3-CNF with at most 3 occurrences of each variable
- [x] Degree-cap experiment at D = m·k·R and D = 3·n·R, compared with brute-force #F
- [x] Fixed-seed formulas under `benches/` and pytest coverage
- See `docs/PROGRESS.md` for the numbers from the fixed seed

### Step 2 — simplified (r, τ, D)-design and blocked sum

- [x] On a hand-built replacement family, group the variables and sum inside each group after the shared measurements are fixed
- [x] Check blocked sum = full-space sum H(t)
- [x] Enumeration of measurements may stay exponential; record how time grows with the number of groups
- Success criterion met: agreement on small n. The popcount tuple count is (r+1)^b, exponential in the number of groups. A subexponential runtime was not required and is not claimed
- The design here is a block partition plus τ explicit linear measurements, with D carried and unused. It is not the matrix M0

### Step 3 — small-parameter replacement matrix

- [x] One machine-built slice of M0 = V0 f_I(W0) and the rows P_u(T), at side lengths 4, 9, 16, and 64. The 4×4 exchange and the row ζ_(00,11) − (T/2) ζ_(10,01) − T ζ_(01,10) are also written out by hand
- [x] T = 0 restores the auxiliary coordinate and, on four explicit placements, restores X_v. H(0) matches brute force. The measured degree sits inside M·(R−1) and inside m·k·R
- [x] A legal design (partition, absolute row sum at most the side length, kernel condition) holds at every t ∈ {1, …, Q} on those placements, including one 3-clause that meets three groups
- [x] Full M0 is absent. The parameter list with R ≥ 32 and 6R ≤ Q forces side length 2^384 at r = 2, above the cap of 64. Schedules are hand-written. The point R = Q = 2 is an exact identity and sits outside the hypothesis R < Q
- Not claimed: Theorem 1.1, a subexponential runtime, or 2^{o(n)}

### Step 4 — Hermite / multi-coefficient interpolation

- [x] Exact rational Hermite interpolation from Hasse coefficients of H. The coefficient of index k is the k-th Taylor coefficient divided by k!. Both were checked on three closed forms
- [x] The jet at a sample t is the formal derivative of the clause product, summed over {0,1}^n in one pass. Derivative orders are not separate passes. The condition count stays D+1, packed into fewer distinct points than the Lagrange path: at most two coefficients per point, at most three, or every coefficient at t = 1
- [x] H(0) matches brute-force #F and `interpolate_H0` on the closed forms, on Step 0 random formulas, and on the Step 1 bench row n6-m4-r1 (7 paired points against 13, #F = 44). On these draws the measured degree equalled the cap
- Not claimed: fewer than D+1 conditions, a subexponential evaluation, Theorem 1.1, or 2^{o(n)}. A full jet of order D does quadratic arithmetic per assignment. Each sample point still costs 2^n

### Step 5 — sparsification and a reduction from general small k-CNF

- [x] Implication-cycle copying at B = 3. A variable of occurrence d > 3 becomes d copies, one per occurrence, plus (¬y_1 ∨ y_2) ∧ … ∧ (¬y_d ∨ y_1). Each copy occurs three times. A variable that already occurs at most three times stays one variable. #F' = #F, with no extra factor of 2
- [x] Brute force agrees on the hand formulas inside n' ≤ 20, including a free bit, an empty clause, one heavy width-3 variable, and a width-2 formula whose clauses sit in two blocks. The two-block row goes through the Step 2 blocked sum. The rows with n' ≤ 7 go through `interpolate_H0`, and on those draws the measured degree equalled M·R. Sixteen units are counted; their interpolation exceeds the assignment budget and is refused
- [x] Recorded, and not counted, when n' > 20. A seeded 3-CNF with n = 5 and m = 12 has n' = 36; the cap-3 bound n' ≥ 3m − 2n is already 26. Sixteen units fit in the brute window and miss the interpolation budget
- The output is mixed width. The exact-width-3 checker rejects the cycle clauses. Step 3 is not applied. B = 3 is the smallest cap this gadget meets
- Not claimed: Theorem 1.1, a subexponential runtime, 2^{o(n)}, or a reduction that keeps n' = n + o(n)

### Step 6 — parameter scan

- [x] For n in a machine-limited range, time brute force against the popcount blocked sum. The grid runs blocks of size 2 for b = 1..7 (n = 2..14), plus rows at r = 1, r = 3, and R = 2. On every row that runs, blocked H(0) matches brute force and the full-cube sum
- [x] Lagrange, paired Hermite, and the Step 1 degree check run on that family where the assignment budget allows, and are refused where it does not. At n = 14 the blocked sum still runs; Lagrange and paired Hermite do not
- [x] Occurrence copying is recorded on the popcount row (the map is the identity), on the seeded 3-CNF with n = 5 and m = 12 (n' = 36, cube not enumerated), on n' = 21, and on a clause that crosses blocks. Step 3 is not applied to that clause. Full M0 at R = 32 and Q = 192 is refused before allocation
- [x] The reproducible record is the operation counts: 2^n, (r+1)^b, the local-evaluation count, and the assignment budget. Wall times are one machine's clock. A refused pipeline has no timing
- Not claimed: Theorem 1.1, a subexponential runtime, 2^{o(n)}, or that the scan closes the distance to a full counting algorithm. That distance is the list in `docs/PROGRESS.md`

### Step 7 — larger slice and a balance-check scaffold

- [x] The slice cap is 256. `auxiliary_problems` passes on two points the cap of 64 refused: r = 3, Q = 2, R = 1, D = 81, and r = 2, Q = 4, R = 3, D = 256. Both have R < Q. T = 0 restoration, degree ≤ R−1, support, the measurement rows, and the kernel condition were checked, and on one explicit placement each, H(0) matches brute force inside M·(R−1) and m·k·R
- [x] `schedule_placement` checks 2(b+w)D² ≤ M and 2Qr² ≤ n, with n ≤ M ≤ bn. Those inequalities need n on the order of D². They fail for every formula this repository compares with brute force, and the call then refuses. It does not draw a placement
- [x] When the inequalities hold (64 unit clauses, r = 2, Q = 1, D = 4), the same call still refuses to choose schedules. The conditional-expectation balancer is not implemented. Caller-supplied schedules, including the four Step 3 placements and the two new points, are placed and `design_problems_at` is run
- [x] `build_auxiliary_slice(2, 192, 32)` and `schedule_placement` at those parameters both refuse before a matrix is allocated. The next pure powers, 729 and 1024, stay above the cap. `verify --quick` runs the D = 81 point and the refusals. `verify --step7` also runs D = 256
- Not claimed: Theorem 1.1, a refutation of ETH, a subexponential runtime, or 2^{o(n)}. A passing balance check is not a schedule

## Not started

No further step is listed. The distance to a full counting algorithm is still the list in `docs/PROGRESS.md`: full M0, the conditional-expectation balancer, sparsification with n' = n + o(n), a measurement space that is not exponential in the number of blocks, and a Hermite evaluation that does not walk the cube. Those constructions are not implemented here.

## Out of scope

- A claim that this repository runs the full counting algorithm, or that ETH is false
- A race against CDCL solvers
- Unchecked runs at n too large to compare with brute force

## File map

```
src/counting_sat_algebra/
  cnf.py, replacements.py, independent.py, self_reduce.py, verify.py   # Step 0
  bounded_occurrence.py                                                # Step 1
  design.py, blocked_sum.py                                            # Step 2
  replacement_matrix.py                                                # Step 3
  hermite.py                                                           # Step 4
  sparsify.py                                                          # Step 5
  parameter_scan.py                                                    # Step 6
  scheduler.py                                                         # Step 7
benches/
docs/PROGRESS.md
```
