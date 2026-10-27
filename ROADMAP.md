# Roadmap

Implement the deterministic exact counter described in references, including the stated time and space bounds. Finite checks are milestones; complete alignment with Theorem 1.1 requires the general-input pipeline and bit-complexity analysis. See docs/REFERENCE_ALIGNMENT.md for the stronger Appendix D bound and the weaker Sections 3–6 / Lean target.

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
- [x] Caller-supplied schedules, including the four Step 3 placements and the two new points, are placed and `design_problems_at` is run. Choosing schedules when the inequalities hold is Step 8
- [x] `build_auxiliary_slice(2, 192, 32)` and `schedule_placement` at those parameters both refuse before a matrix is allocated. The next pure powers, 729 and 1024, stay above the cap. `verify --quick` runs the D = 81 point and the refusals. `verify --step7` also runs D = 256
- Not claimed: Theorem 1.1, a refutation of ETH, a subexponential runtime, or 2^{o(n)}. A passing balance check is not a schedule

### Step 8 — original balancer and legacy Boolean budget

Historical scaffold from 0.9.0. The explicit Boolean mode reproduces its budget stops; Step 9 implements the exact reference scheduler and default rational check.

- [x] Successive choice runs only when the Lemma 4.5 inequalities hold and the slice is under the cap. At each variable and each clause it tries every symbol in {0, …, r−1}^Q and keeps the symbol that minimizes a quadratic load potential. Ties break toward the lexicographically first symbol. Nothing is sampled. The schedule cube is not enumerated
- [x] In the explicit `design_method="boolean"` mode, on 64 disjoint unit clauses, r = 2, Q = 1, D = 4, the search meets label multiplicity 16 ≤ 2M/D = 32 and largest group 32, with 32·r ≤ 2n. `design_problems_at` is not called: an outside group has size 32, and that walk enumerates 2^32 Boolean points, above the budget 4096. `place_formula` is not called. The same kind of refusal, with 1 ≤ R < Q holding, is what 1024 disjoint units give at r = 2, Q = 2, R = 1, D = 16 (outside group 512, multiplicity 64 ≤ 128)
- [x] Formulas that fail the inequalities, including the brute-comparable ones, still refuse with no placement and without running the search. Full M0 and side length 1024 still refuse before allocation. A direct call on eight units at r = 2, Q = 2, which is not admitted, passes `design_problems_at` and the T = 0 restoration. `schedule_placement` does not make that call
- [x] `verify --step8` and `verify --quick` record these rows. Refusals are not timed. The direct check on eight units took about 0.01 seconds on one run. That clock is not a complexity result
- Not claimed: Theorem 1.1, a refutation of ETH, a subexponential runtime, or 2^{o(n)}. The search is not a proof of Lemma 4.5. The binomial tail is not checked. A budget refusal is not a placement

### Step 9 — reference conditional expectation and rational design

- [x] Minimize the exact reference potential Phi/(M(b+w)) + Xi/(Qn), including shared-incidence dependencies
- [x] Check each prefix expectation against all completions on tiny inputs; independently check the final potential
- [x] Eliminate each group's measurement matrix once and check the entire rational kernel through row-span membership
- [x] Complete a real admitted n=1024, Q=2, R=1 placement without a Boolean cube check
- [x] Check Lemma 4.7's exact binomial tail. The n=1024 point fails it; n=16384, Q=3, R=2, tau=1/2 meets every numerical hypothesis

### Step 10 — Section 3 common counting routine

- [x] Stream integer measurement tuples and scan each local Boolean cube twice, with no attained-measurement table
- [x] Recover truncated Taylor jets using U subset E and fixed-bit restrictions, then restore H(0) by Hermite interpolation
- [x] Compare all requested coefficients and counts on the four Step 3 examples with independent cube oracles
- [x] Check second-order recovery, negative literals, empty clauses, free variables, invalid designs and pre-enumeration budgets
- This counts supplied finite designs. Their measurement space can still be exponential; no general-formula complexity bound follows from these runs

### Step 11 — structured auxiliary operator

- [x] Query M0 and replacement-polynomial entries by Fourier/Krawtchouk sums over swap orbits
- [x] Compare every entry with dense slices of side 4, 9, 16, 64 and 81
- [x] Query coordinates at r=2, Q=192, R=32, D=2^384 with no dense allocation
- [x] Describe individual integer measurement modes and refuse oversized enumeration
- Step 12 connects this descriptor to finite sparse placements and the common counting entry point

### Step 12 — sparse admitted placement

- [x] Store occupied polynomial and measurement coordinates without M×n dense arrays
- [x] Bound replacement support, mode-coordinate queries, pullback candidates and sparse elimination before resource overuse
- [x] Check restoration, common denominator, coefficient norms and exact rational kernels
- [x] Complete n=16384, r=2, Q=3, R=2 placement with all Lemma 4.7 numerical hypotheses satisfied
- [x] Record 59392 polynomial entries and 14336 measurement entries versus 536870912 dense polynomial slots
- [x] Plan the common count and refuse before enumeration: q=5461, J=5462, zero Boolean vectors examined

### Step 13 — DIMACS and independent components

- [x] Parse strict DIMACS with declared free variables, multiline clauses and empty clauses
- [x] Normalize clauses and identify independent components; place each component afresh
- [x] Count each nontrivial component through the sparse shared-measurement backend; cap optional memoization
- [x] Add actual count CLI, JSON reporting, explicit budgets and exact hex output for huge counts
- [x] Count the eight-variable example as 196 and the 32769-variable disjoint 3-CNF as 7^10923
- [x] Verify no cube-oracle backend, free-bit factors, disabled memoization and refusals without partial counts
- This completes a finite supplied-CNF counter within component budgets. General large connected instances and the reference's asymptotic bounds remain open

## Remaining work toward the reference conclusion

- [x] Connect structured auxiliary entries and measurement descriptors to finite sparse placed formulas
- [ ] Choose and run the full asymptotic parameter family
- [ ] Implement Section 5's sparse polynomial systems and coefficient algebra
- [ ] Implement Section 6's equation selection and disjoint sparsification with the required size control
- [ ] Implement Appendix D's tensor contraction to reach Theorem 1.1's sharper bound
- [ ] Establish exact bit-operation and space bounds; complete the Lean machine and proof without the existing sorry

Unsupported claims of theorem completion are not an acceptance criterion. Results must distinguish counted formulas, checked hypotheses or coordinates, and budget refusals.

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
  scheduler.py, conditional_balance.py, linear_algebra.py                # Steps 7–9
  reference_parameters.py                                               # reference hypothesis checks
  common_count.py                                                       # Step 10
  structured_auxiliary.py, reference_verify.py                           # Step 11 and alignment checks
  sparse_placement.py                                                   # Step 12
  component_count.py, dimacs.py, cli.py                                  # Step 13
benches/
examples/count_demo.cnf
docs/PROGRESS.md
```
