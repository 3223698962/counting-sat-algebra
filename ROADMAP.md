# Roadmap

Build the pieces of the counting algorithm on instances small enough to compare with brute force. A step is done when that comparison passes and the cost matches the step's own analysis. Refuting ETH is not a success criterion for this repository.

## Done

### Step 0 — counting polynomial and interpolation

- [x] Clause polynomial Φ_e and brute-force #F
- [x] Product identity on disjoint variable groups
- [x] Random replacements of degree ≤ R with value X_v at T = 0; Lagrange recovery of H(0)
- [x] Self-reduction demo that searches for one satisfying assignment
- Cost: each H(t), t ≠ 0, still enumerates 2^n assignments

## Not started

### Step 1 — bounded-occurrence testbed

- [ ] Generator and checker for 3-CNF with at most 3 occurrences of each variable
- [ ] Degree-cap experiment at D = m·k·R and D = 3·n·R, compared with brute-force #F
- [ ] Fixed-seed formulas under `benches/` and pytest coverage

### Step 2 — (r, τ, D)-design and blocked sum

- [ ] On a hand-built replacement family, group the variables and sum inside each group after the shared measurements are fixed
- [ ] Check blocked sum = full-space sum H(t)
- [ ] Enumeration of measurements may stay exponential; record how time grows with the number of groups
- Success criterion: agreement on small n. A subexponential runtime is not required

### Step 3 — a runnable slice of the low-degree replacement / matrix construction

- [ ] One small-parameter slice that can be computed by hand or by the machine
- [ ] Check the T = 0 restoration, the degree cap, and that a legal design exists at a few t
- [ ] If the full matrix M0 is too heavy, write down the obstacle and the simplification

### Step 4 — Hermite / multi-coefficient interpolation

- [ ] Several Taylor or Hasse coefficients at each sample point
- [ ] Fewer sample points, compared with the degree cap

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
  # bounded_occurrence.py                                                # Step 1, absent
benches/          # Step 1
docs/PROGRESS.md  # Step 1
```
