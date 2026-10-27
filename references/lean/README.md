# Lean materials (separate effort)

These files are a local copy of a Lean formalization effort. They are not part of the Python package `counting_sat_algebra`. This repository does not build them, import them, or check them with Lean.

- [`378.md`](378.md) — scope note from `openai/math`, `lean/docs/378.md`. It states the aim of that formalization: for fixed \(k \ge 2\), a deterministic machine that exactly counts satisfying assignments of a canonically encoded \(k\)-CNF, within the time and visited-tape bounds written in the note. The October 25 manuscript’s sharper counting bound is quoted there as scope for that effort.
- [`SubexponentialSAT.lean`](SubexponentialSAT.lean) — scaffold linked from that note as `ComparatorChallenges/SubexponentialSAT.lean` (namespace `OAI.Algorithms.SubexponentialSAT`). It defines CNF evaluation, a canonical encoding, a primitive-machine trace, and visited tape cells. The main statement, `exact_count_and_sat`, ends in `sorry`.

In Lean, `sorry` is a placeholder for an unfinished proof. The declaration can be type-checked while the proof is absent. This file does not prove `exact_count_and_sat`, and it does not prove Theorem 1.1 of the preprint.

The Python package version 0.14.0 checks finite instances of interpolation, a blocked sum, a matrix slice of side length at most 256, Hermite recovery of \(H(0)\), occurrence copying, a parameter scan, exact reference conditional-expectation schedules, rational design checks, streamed shared-measurement Taylor recovery, and structured full-parameter coordinate queries. It also constructs a sparse admitted placement at n=16384, Q=3, R=2, and counts supplied DIMACS CNFs through independently placed components. A disconnected 32769-variable instance has exact count 7^10923; the global sum on the 16384-variable placement is refused before enumeration. The complete general-input and asymptotic pipeline remains open. These finite results are separate from a proof of this Lean statement.
