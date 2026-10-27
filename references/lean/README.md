# Lean materials (separate effort)

These files are a local copy of a Lean formalization effort. They are not part of the Python package `counting_sat_algebra`. This repository does not build them, import them, or check them with Lean.

- [`378.md`](378.md) — scope note from `openai/math`, `lean/docs/378.md`. It states the aim of that formalization: for fixed \(k \ge 2\), a deterministic machine that exactly counts satisfying assignments of a canonically encoded \(k\)-CNF, within the time and visited-tape bounds written in the note. The October 25 manuscript’s sharper counting bound is quoted there as scope for that effort.
- [`SubexponentialSAT.lean`](SubexponentialSAT.lean) — scaffold linked from that note as `ComparatorChallenges/SubexponentialSAT.lean` (namespace `OAI.Algorithms.SubexponentialSAT`). It defines CNF evaluation, a canonical encoding, a primitive-machine trace, and visited tape cells. The main statement, `exact_count_and_sat`, ends in `sorry`.

In Lean, `sorry` is a placeholder for an unfinished proof. The declaration can be type-checked while the proof is absent. This file does not prove `exact_count_and_sat`, and it does not prove Theorem 1.1 of the preprint.

The Python package checks small instances of isolated pieces: interpolation, a blocked sum, a matrix slice of side length at most 256, Hermite recovery of \(H(0)\), occurrence copying, a parameter scan, and a successive-choice balancer that refuses when its Boolean-point budget or the balance inequalities fail. Agreement with brute force on those instances is a different artifact from this Lean statement.
