# counting-sat-algebra

Small Python package that checks the counting-polynomial skeleton for #3-SAT on tiny instances.

> **Honesty banner**
>
> Step 0 checks the following on instances small enough to brute-force. Write #F as a sum of clause polynomials, replace each literal input by a rational linear form that restores the original variable at T = 0, and recover H(0) = #F by Lagrange interpolation.
>
> Every value H(t) for t ≠ 0 is still a sum over all 2^n assignments. The package does not contain an (r, τ, D)-design, a blocked sum, the matrix M0, or sparsification. It does not refute the exponential-time hypothesis, and it does not claim a 2^{o(n)} algorithm.

## Layout

```
src/counting_sat_algebra/   # import name counting_sat_algebra
tests/                      # pytest
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
| `self_reduce.py` | Optional search for one satisfying assignment by fixing bits |

The degree cap used for a width-k formula is `m * k * R`. Interpolation uses the nodes T = 1, …, D+1.

## Next step

Step 1 in `ROADMAP.md`: a generator and checker for 3-CNFs with at most three occurrences of each variable, and a degree-cap experiment. That step is not in this tree yet.
