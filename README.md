# counting-sat-algebra

Small Python package that checks the counting-polynomial skeleton for #3-SAT on tiny instances.

> **Honesty banner**
>
> This repository checks five things on instances small enough to brute-force:
>
> 1. **Step 0.** Write #F as a sum of clause polynomials, replace each literal input by a rational linear form that restores the original variable at T = 0, and recover H(0) = #F by Lagrange interpolation.
> 2. **Step 1.** Generate and recognize 3-CNFs in which every variable occurs in at most three clauses, and check the degree caps D = m·k·R and D = 3·n·R on those formulas.
> 3. **Step 2.** On a hand-built replacement family whose higher coefficients lie in τ shared linear measurements, partition the variables into blocks of size at most r and check that the blocked sum equals the full sum H(t). The integer D is carried as a degree cap and is not used to build a matrix.
> 4. **Step 3.** Build the replacement matrix M0 = V0 f_I(W0) and the rows that come from it, only when the side length r**(2Q) is at most 64. Check that T = 0 restores the variable, that the degree stays inside M·(R−1) and m·k·R, and that a legal design holds at every sample t on four explicit placements.
> 5. **Step 4.** Read several Hasse coefficients of H at each sample point and recover H(0) = #F from fewer distinct points than the D+1 value samples of the Lagrange path. The number of conditions stays D+1.
>
> H(t) for a general replacement is still a sum over all 2^n assignments. A Hermite node is one such pass: the jet is produced inside it, and a full jet of order D does quadratic arithmetic per assignment. The Step 2 measurement enumeration is exponential in the number of blocks: for one popcount measurement and blocks of size r it walks (r+1)^b tuples. The matrix at the parameter list R ≥ 32 and 6R ≤ Q is not built: at r = 2 that side length is 2^384. Sparsification is not in this tree. Theorem 1.1 is not implemented. Nothing here refutes the exponential-time hypothesis or runs in 2^{o(n)} time.

## Layout

```
src/counting_sat_algebra/   # import name counting_sat_algebra
tests/                      # pytest
benches/                    # fixed-seed bounded-occurrence formulas
docs/PROGRESS.md            # what was actually checked
```

Core arithmetic uses `fractions.Fraction` from the standard library.

## Run

```bash
python -m pip install -e ".[dev]"
pytest -q
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra.verify --step2
python -m counting_sat_algebra.verify --step3
python -m counting_sat_algebra.verify --step4
```

`python -m counting_sat_algebra.verify` without `--quick` draws more random 3-CNFs (n up to 10) and extends the Step 2 group scan by one block. That path is exact but slow, because each Step 0 sample still costs 2^n. `--step2` runs only the blocked-sum checks. `--step3` runs only the tiny replacement-matrix checks. `--step4` runs only the Hermite checks. `--quick` includes Steps 3 and 4; those slices are small.

Exit code 0 and the line `ALL TESTS PASSED` mean the checks in that run agreed with brute force.

## What the checks correspond to

| Code | Idea being checked |
|------|--------------------|
| `cnf.py` | Clause polynomial Φ_e and brute-force #F |
| `independent.py` | If no clause crosses variable groups, #F = ∏ N_g |
| `replacements.py` | Degree-≤ R replacements with value X_v at T = 0; H(0) by Lagrange |
| `bounded_occurrence.py` | 3-CNF, ≤ 3 occurrences per variable; degree caps m·k·R and 3·n·R |
| `design.py` | Simplified (r, τ, D) block partition and τ shared linear measurements |
| `blocked_sum.py` | Hand-built replacements in that span; blocked sum versus full H(t) |
| `replacement_matrix.py` | Tiny M0 = V0 f_I(W0) and the rows built from it; not the full matrix |
| `hermite.py` | Hasse jets of H and Hermite recovery of H(0); still one 2^n pass per node |
| `self_reduce.py` | Optional search for one satisfying assignment by fixing bits |

The degree cap used for a width-k formula is `m * k * R`. On the bounded-occurrence testbed, `M = Σ |I_e| ≤ 3n`, so `M R ≤ 3 n R` is a second, sometimes looser, cap. Step 1 records both.

Step 2's design is a simplification. Blocks partition the variables, each block has size at most r, and τ rational linear measurements are shared across blocks. Clauses may not cross blocks. The replacement is hand-built so that

`⟨P, X⟩ = X_v + Σ_{j≥1} T^j ⟨α_j, μ(X)⟩`.

After the local measurements are fixed, the sum factors into a product of block sums. That identity is checked against the full sum on small n. D is stored and not read. A measurement space small enough for a subexponential count is not part of this step.

Step 3 builds M0 = V0 f_I(W0) when D = r**(2Q) ≤ 64. The rows restore the auxiliary coordinate at T = 0 and have degree at most R−1. On four hand-written placements they restore X_v, H(0) matches brute force, and the measured degree sits inside both caps. One of those placements is a single 3-clause whose variables sit in three different groups; Step 2 rejects that shape, and the Step 3 kernel condition accepts it. The point Q = R = 2 is checked as an identity. It does not satisfy R < Q. The matrix at R ≥ 32 and 6R ≤ Q is not built.

Step 4 differentiates the clause product in the jet ring and sums over the cube. The Hasse coefficient of index k is the Taylor coefficient divided by k!. Newton interpolation with repeated nodes then evaluates H(0). Three patterns are checked: at most two coefficients per point, at most three, and every coefficient at the single point t = 1. Each pattern uses D+1 conditions, and each uses fewer distinct points than Lagrange on t = 1, …, D+1. On the recorded draws the measured degree equalled the cap, so no condition was dropped for a slack bound. The single-point jet replaces D+1 passes by one pass whose arithmetic per assignment grows quadratically with D.

## Next step

Step 5 in `ROADMAP.md`: sparsification, and a reduction from a general small k-CNF into the bounded-occurrence testbed. That step is not in this tree.
