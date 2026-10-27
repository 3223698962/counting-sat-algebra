# counting-sat-algebra

Small Python package that checks the counting-polynomial skeleton for #3-SAT on tiny instances.

> **Honesty banner**
>
> This repository checks seven things on instances small enough to brute-force:
>
> 1. **Step 0.** Write #F as a sum of clause polynomials, replace each literal input by a rational linear form that restores the original variable at T = 0, and recover H(0) = #F by Lagrange interpolation.
> 2. **Step 1.** Generate and recognize 3-CNFs in which every variable occurs in at most three clauses, and check the degree caps D = m·k·R and D = 3·n·R on those formulas.
> 3. **Step 2.** On a hand-built replacement family whose higher coefficients lie in τ shared linear measurements, partition the variables into blocks of size at most r and check that the blocked sum equals the full sum H(t). The integer D is carried as a degree cap and is not used to build a matrix.
> 4. **Step 3.** Build the replacement matrix M0 = V0 f_I(W0) and the rows that come from it, only when the side length r**(2Q) is at most 64. Check that T = 0 restores the variable, that the degree stays inside M·(R−1) and m·k·R, and that a legal design holds at every sample t on four explicit placements.
> 5. **Step 4.** Read several Hasse coefficients of H at each sample point and recover H(0) = #F from fewer distinct points than the D+1 value samples of the Lagrange path. The number of conditions stays D+1.
> 6. **Step 5.** Replace a variable that occurs more than three times by one copy per occurrence, tied together by an implication cycle, and check #F' = #F where both formulas fit in the brute-force window. Pass the reduced formula through `interpolate_H0`, and through the Step 2 blocked sum when the source clauses already sit inside blocks.
> 7. **Step 6.** On a fixed machine-limited grid of n, time brute force against the popcount blocked sum, and record the other steps as a match or as a refusal. A refusal is what the grid does when the cube passes the brute-force wall, the assignment budget, or the measurement-tuple cap, and when the construction is unfinished: full M0, a crossing clause with no Step 3 placement, copying that pushes n' past 20, a Hermite jet that would still walk the cube. Wall times are one machine's clock. The reproducible record is the operation counts.
>
> H(t) for a general replacement is still a sum over all 2^n assignments. A Hermite node is one such pass: the jet is produced inside it, and a full jet of order D does quadratic arithmetic per assignment. The Step 2 measurement enumeration is exponential in the number of blocks: for one popcount measurement and blocks of size r it walks (r+1)^b tuples. The matrix at the parameter list R ≥ 32 and 6R ≤ Q is not built: at r = 2 that side length is 2^384. Step 5's copying can make the reduced length larger than the source length; past n' = 20 the model counts are recorded as unchecked. The cycle clauses have width 2, so the exact-width-3 checker does not accept a formula that used the cycle. The Step 6 scan writes those stops down. It does not remove them, and it does not time a construction it refused. Theorem 1.1 is not implemented. Nothing here refutes the exponential-time hypothesis or runs in 2^{o(n)} time.

## Layout

```
src/counting_sat_algebra/   # import name counting_sat_algebra
tests/                      # pytest
benches/                    # fixed-seed formulas and the Step 6 scan record
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
python -m counting_sat_algebra.verify --step5
python -m counting_sat_algebra.verify --step6
```

`python -m counting_sat_algebra.verify` without `--quick` draws more random 3-CNFs (n up to 10), extends the Step 2 group scan by one block, and runs the full Step 6 grid. That path is exact but slow, because each Step 0 sample still costs 2^n. `--step2` runs only the blocked-sum checks. `--step3` runs only the tiny replacement-matrix checks. `--step4` runs only the Hermite checks. `--step5` runs only the occurrence-copying checks. `--step6` runs only the parameter scan, on the full grid. `--quick` includes Steps 3, 4, 5, and the small Step 6 grid. Step 5 records a row instead of counting it when the reduced n is above 20. Step 6 records a cell as refused, and does not time it, when the construction is unfinished or the enumeration is past a rail.

Exit code 0 and the line `ALL TESTS PASSED` mean the checks in that run agreed with brute force, or were the refusals the scan records. A Step 5 row whose reduced n is above 20 is recorded as unchecked and still passes under that line. A Step 6 refusal is the same kind of recorded stop.

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
| `sparsify.py` | Implication-cycle copying to occurrence ≤ 3; #F' = #F inside the brute window |
| `parameter_scan.py` | Machine-limited grid: matches, refusals, operation counts; not a complexity result |
| `self_reduce.py` | Optional search for one satisfying assignment by fixing bits |

The degree cap used for a width-k formula is `m * k * R`. On the bounded-occurrence testbed, `M = Σ |I_e| ≤ 3n`, so `M R ≤ 3 n R` is a second, sometimes looser, cap. Step 1 records both.

Step 2's design is a simplification. Blocks partition the variables, each block has size at most r, and τ rational linear measurements are shared across blocks. Clauses may not cross blocks. The replacement is hand-built so that

`⟨P, X⟩ = X_v + Σ_{j≥1} T^j ⟨α_j, μ(X)⟩`.

After the local measurements are fixed, the sum factors into a product of block sums. That identity is checked against the full sum on small n. D is stored and not read. A measurement space small enough for a subexponential count is not part of this step.

Step 3 builds M0 = V0 f_I(W0) when D = r**(2Q) ≤ 64. The rows restore the auxiliary coordinate at T = 0 and have degree at most R−1. On four hand-written placements they restore X_v, H(0) matches brute force, and the measured degree sits inside both caps. One of those placements is a single 3-clause whose variables sit in three different groups; Step 2 rejects that shape, and the Step 3 kernel condition accepts it. The point Q = R = 2 is checked as an identity. It does not satisfy R < Q. The matrix at R ≥ 32 and 6R ≤ Q is not built.

Step 4 differentiates the clause product in the jet ring and sums over the cube. The Hasse coefficient of index k is the Taylor coefficient divided by k!. Newton interpolation with repeated nodes then evaluates H(0). Three patterns are checked: at most two coefficients per point, at most three, and every coefficient at the single point t = 1. Each pattern uses D+1 conditions, and each uses fewer distinct points than Lagrange on t = 1, …, D+1. On the recorded draws the measured degree equalled the cap, so no condition was dropped for a slack bound. The single-point jet replaces D+1 passes by one pass whose arithmetic per assignment grows quadratically with D.

Step 5 copies variables. For the default cap B = 3, a variable in d > 3 clauses becomes d variables, and the copies are forced equal by a cycle of d binary clauses. The recovery on a counted instance is #F' = #F. `interpolate_H0` is run when the sample window fits a fixed assignment budget. The Step 2 blocked sum is run when the caller supplies a partition of the source variables that already contains every clause; copies stay in the owner's block, and a clause that already crossed blocks still crosses. A 3-CNF with 3m − 2n > 20 has reduced length above the brute-force wall n ≤ 20, and that equality is not checked. The numbers are in `docs/PROGRESS.md`.

Step 6 does not add an algorithm. It runs the pieces above on one fixed grid and writes down the result of each cell. Where the popcount blocks fit, blocked H(0) equals the brute-force count. Where a sample window, a measurement tuple count, or n itself is past a rail already used in this repository, the cell is refused before that loop. Full M0 is one of those refusals: the call is the parameter list R ≥ 32 and 6R ≤ Q, and it raises before the matrix is allocated. A copied clause that still meets two blocks is another: Step 2 refuses it, and the scan does not invent a Step 3 placement. The single-point Hermite jet is not run. It would still be one pass over the cube. Pytest locks the small grid. `benches/parameter_scan.json` is the full grid, including wall times from the run that wrote the file.

## Next step

There is no further step in `ROADMAP.md`. The remaining distance to a full counting algorithm is the unfinished list in `docs/PROGRESS.md`: full M0, a scheduler, sparsification that keeps n' close to n, a measurement space that is not exponential in the number of blocks, and a Hermite evaluation that costs less than a pass over the cube. None of those is in this tree. Theorem 1.1 is not implemented.
