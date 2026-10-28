# Roadmap

Implement the deterministic exact counter described in references, including the stated time and space bounds. Finite checks are milestones; complete alignment with Theorem 1.1 requires the general-input pipeline and bit-complexity analysis. See [docs/REFERENCE_ALIGNMENT.md](docs/REFERENCE_ALIGNMENT.md) for the stronger Appendix D bound and the weaker Sections 3–6 / Lean target.

This file lists remaining work only. The Step 0–21 checklist lives in [docs/steps/completed.md](docs/steps/completed.md). The layering of core, backends, and witnesses is in [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md). Package 0.22.1 does not add an algorithm. `theorem_1_1_complete`, `lean_proof_complete`, and `recursive_contraction` stay false. There is no general \(2^{o(n)}\) algorithm in this repository, and this repository does not refute ETH.

## Remaining work

- [ ] Implement D.4.5's weighted modes
- [ ] Implement D.4.6's packet tensor coefficients and immutable creation-time descriptions
- [ ] Implement D.4.7's recursive child contractions and establish edge, degree and precision invariants throughout recursion. Step 15's finite Taylor elimination is not this algorithm
- [ ] Implement D.3's incidence-monomial network encoding and connect the general-density Section 6 queries to the completed modular recursion
- [ ] Establish exact bit-operation and space bounds; complete the Lean machine and proof without the existing sorry

Unsupported claims of theorem completion are not an acceptance criterion. Results must distinguish counted formulas, checked hypotheses or coordinates, and budget refusals.
