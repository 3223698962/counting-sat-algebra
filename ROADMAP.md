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

## Not started

### Step 5 — sparsification and a reduction from general small k-CNF

- [ ] Counting-preserving sparsification or variable copying at small n
- [ ] Pipe a general small k-CNF into a bounded-occurrence instance and then through Steps 2 and 3

### Step 6 — parameter scan

- [ ] For n in a machine-limited range, time brute force against a blocked sum
- [ ] Record which steps match brute force, which stopped on an unfinished construction, and the remaining distance to a full counting algorithm

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
benches/
docs/PROGRESS.md
```
