# Lean reference materials

Python package version **0.26.0** uses the newly supplied Lean library as its implementation reference. Source availability, Lean kernel verification and completion of the Python algorithm are separate facts.

## Supplied library

The current source root is [`math/lean/`](math/lean/), with the algorithm under [`OAI/Algorithms/SubexponentialSAT/`](math/lean/OAI/Algorithms/SubexponentialSAT/). Its toolchain is `leanprover/lean4:v4.34.1`; Lake dependencies are pinned in `lake-manifest.json`.

- [`Main.lean`](math/lean/OAI/Algorithms/SubexponentialSAT/Main.lean) imports the current algorithm endpoints.
- [`SubexpOrderedRawMainTheorem.lean`](math/lean/OAI/Algorithms/SubexponentialSAT/SubexpOrderedRawMainTheorem.lean), `mainTheorem_explicit`, gives the original externally encoded SAT signature. Its proof is a theorem reference, rather than a `by` block.
- [`SqrtLogSpace/Final.lean`](math/lean/OAI/Algorithms/SubexponentialSAT/SqrtLogSpace/Final.lean), `compiled_exact_count_and_sat`, gives one fixed compiled program, exact count and SAT tapes, a canonical binary count word of length at most n+1, and a polynomial visited-tape bound at every prefix. `exact_count_and_sat` exposes the primitive-machine signature; `exact_count_eth` transports the bounds to n and the clause count.
- `formalization.yaml` lists the SqrtLogSpace exact-count endpoint.

The exact count's exponent uses the square of log log n, as explicitly defined in `SqrtLog/StatementScale.lean` and `SqrtLogSpace/Asymptotics.lean`. It is the weaker time target compared with the local manuscript's Appendix D bound. The roadmap follows this supplied counting program first.

## Local audit and kernel verification

`python -m counting_sat_algebra.lean_reference --output artifacts/lean-reference-audit.json` follows the current endpoints' local import closure and records a content hash. The current static audit covers **3969 local modules**, with no missing local imports, code-level `sorry`/`admit`, or local `axiom` declarations. Nested comments and strings are excluded from these searches.

The audit records `kernel_verified=false`. This environment has not compiled the supplied library with Lean; Mathlib and the toolchain still need to be resolved for that check. In a prepared environment, run these commands from `math/lean/`:

```bash
lake build OAI.Algorithms.SubexponentialSAT.Main
lake build OAI.Algorithms.SubexponentialSAT.SubexpOrderedRawMainTheorem
lake env lean OAI/Algorithms/SubexponentialSAT/SqrtLogSpace/Final.lean
```

The final file includes `#print axioms` for its endpoints. Save that actual output when carrying out kernel verification. A source scan does not establish the theorem's transitive axiom set.

## Historical scaffold

[`SubexponentialSAT.lean`](SubexponentialSAT.lean) and [`378.md`](378.md) remain historical material. The single-file scaffold's `exact_count_and_sat` still ends in `sorry`; it is not the current library endpoint. The statement that all supplied Lean material is merely that scaffold is obsolete.

## Python transfer

Step 23 implements the binary syntax and streaming disjoint frontier. Step 24 implements ternary Newton digits, shared-marker truncation and finite accepted/marked query counts, with a native algebraic leaf available in the DIMACS command. Step 25 implements the scheduled ternary sparse code with regenerated rows and the violation-mass greedy selector, available as the hash leaf. Step 26 implements the early high-incidence threshold, the high-variable frontend and the streamed conditioned rank loop, which the hash leaf can use for every query. The counting extension, packing and streamed Taylor leaf remain open. Python scalar operation counts and logical encoding charges are not primitive Lean-machine step or visited-tape measurements.

`theorem_1_1_complete`, `lean_proof_complete` and `recursive_contraction` stay false in Python reports. In particular, `lean_proof_complete=false` records the pending local verification and Python refinement; it does not assert that the newly supplied source has a missing theorem proof.
