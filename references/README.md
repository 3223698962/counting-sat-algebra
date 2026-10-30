# References

Offline copies for reading this repository. `pytest` and `python -m counting_sat_algebra.verify` do not execute these files.

| Path | What it is |
|------|------------|
| [`paper/`](paper/README.md) | Preprint *The Exponential Time Hypothesis is False* (OpenAI; 17 October 2026, revised 25 October 2026): PDF and plain-text extract |
| [`lean/`](lean/README.md) | Supplied SubexponentialSAT library with final exact-count and space proof source, plus the separately retained historical single-file scaffold |
| [`others/`](others/) | Deterministic factoring application note: `DeterministicFactoring.tex`, its arXiv 2610.00600v1 PDF and source zip; the application mapping is recorded in [`docs/APPLICATION_FACTORING.md`](../docs/APPLICATION_FACTORING.md) |

The supplied library is the active Python implementation reference in version 0.30.0. Its static local-import audit passes; local Lean kernel compilation and Python refinement remain pending. The old scaffold's `sorry` does not describe the new library. Source availability, kernel verification and completion of the Python algorithm are recorded separately in [`lean/README.md`](lean/README.md) and [`docs/STATUS.md`](../docs/STATUS.md).
