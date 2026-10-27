# counting-sat-algebra

Small Python package that checks the counting-polynomial skeleton for #3-SAT on tiny instances.

> **Honesty banner**
>
> This repository checks two things on instances small enough to brute-force:
>
> 1. **Step 0.** Write #F as a sum of clause polynomials, replace each literal input by a rational linear form that restores the original variable at T = 0, and recover H(0) = #F by Lagrange interpolation.
> 2. **Step 1.** Generate and recognize 3-CNFs in which every variable occurs in at most three clauses, and check the degree caps D = m·k·R and D = 3·n·R on those formulas.
>
> Every value H(t) for t ≠ 0 is still a sum over all 2^n assignments. The package does not contain an (r, τ, D)-design, a blocked sum, the matrix M0, or sparsification. Theorem 1.1 is not implemented. Nothing here refutes the exponential-time hypothesis or runs in 2^{o(n)} time.

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
```

`python -m counting_sat_algebra.verify` without `--quick` draws more random 3-CNFs (n up to 10). That path is exact but slow, because each sample still costs 2^n.

Exit code 0 and the line `ALL TESTS PASSED` mean the checks in that run agreed with brute force.

## What the checks correspond to

| Code | Idea being checked |
|------|--------------------|
| `cnf.py` | Clause polynomial Φ_e and brute-force #F |
| `independent.py` | If no clause crosses variable groups, #F = ∏ N_g |
| `replacements.py` | Degree-≤ R replacements with value X_v at T = 0; H(0) by Lagrange |
| `bounded_occurrence.py` | 3-CNF, ≤ 3 occurrences per variable; degree caps m·k·R and 3·n·R |
| `self_reduce.py` | Optional search for one satisfying assignment by fixing bits |

The degree cap used for a width-k formula is `m * k * R`. On the bounded-occurrence testbed, `M = Σ |I_e| ≤ 3n`, so `M R ≤ 3 n R` is a second, sometimes looser, cap. Step 1 records both.

## Next step

Step 2 in `ROADMAP.md`: an (r, τ, D)-design and a blocked sum compared with the full sum on a hand-built replacement family. That step is not in this tree.
