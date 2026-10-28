# Architecture

当前包版本 **0.22.1**。这一页区分仍在推进的 Appendix D 路线和已经落地的有限见证。Theorem 1.1 没有实现，本仓库不否定 ETH，也没有一般的 \(2^{o(n)}\) 算法。`theorem_1_1_complete`、`lean_proof_complete` 与 `recursive_contraction` 保持 false。Lean 主定理 `exact_count_and_sat` 仍以 `sorry` 结束。

Package 0.22.1 is a split of the reference-check code plus repository gates. The algebraic routines are unchanged.

## Layers and step numbers

Step numbers are historical milestones. The code is easier to read as three layers.

**Core.** Formulas, exact rationals, designs, and rings. These modules do not claim a runtime bound: `cnf.py`, `dimacs.py`, `replacements.py`, `linear_algebra.py`, `design.py`, `modular_ring.py`, `reference_parameters.py`.

**Backends.** Finite procedures that return a count, a checked object, or a refusal before enumeration:

| Backend | Module | What it actually does |
| --- | --- | --- |
| Shared measurement | `common_count.py` | Section 3 local sums and Taylor recovery |
| Disjoint supports | `support_count.py` | Counts one placed formula by replacement support |
| Finite Taylor elimination | `tensor_count.py` | Step 15 tables inside an explicit budget. Not Appendix D recursion |
| DIMACS components | `component_count.py`, `cli.py` | Re-places each component. Backends are `common` or `tensor` only |
| Section 5–6 | `sparse_systems.py`, `algebra_count.py`, `equation_selection.py` | Budgeted finite algebra on a declared \(n\) |
| Appendix D leaf and recovery | `modular_network.py`, `modular_recovery.py` | Direct edge scans and tapered Hermite recovery |
| Collision layout | `modular_packing.py` | Labels, packets, and a child port topology. No tensor evaluation and no recursive call |

**Witnesses.** `reference_verify/` runs the finite checks and writes the JSON reports. One module per `run_*_suite` family, plus `protocol.py` for the COUNTED / CHECKED / REFUSED row and the false completion flags, and `report.py` for the writer. `__init__.py` still exports `REFERENCE_SUITES` and `write_reference_report`. A budget refusal does not return a partial count.

`verify.py` is the CLI orchestrator for Steps 0–21. It is intentionally still one module. Steps 17–21 are not wired into the DIMACS `count` command.

## What is live, and what is a finite witness

The open path toward the reference conclusion is Appendix D recursive contraction:

- D.4.5 weighted modes are not implemented.
- D.4.6 packet-tensor coefficients and immutable creation-time descriptions are not implemented. Step 21 builds the child port topology only.
- D.4.7 recursive child contraction, with edge, degree, and precision invariants, is not implemented. `recursive_contraction` stays false. Step 15's Taylor elimination is a different, finite algorithm.
- D.3 incidence-monomial encoding, exact bit complexity, and the Lean proof (`sorry`) remain open.

Steps 0–19 and the Step 20–21 objects that do run are finite witnesses or scaffolding for that path. A passing placement, a counted small formula, or a checked integer schedule is not an asymptotic bound and is not Theorem 1.1. The Sections 3–6 route has the weaker exponent quoted in `REFERENCE_ALIGNMENT.md`. It is not the Appendix D bound.

## Module size

CI rejects a Python file under `src/counting_sat_algebra/` that exceeds these limits. Sizes are file bytes.

| Limit | Bytes | Rule |
| --- | --- | --- |
| Hard | 65536 | Every module. A file over this fails the build. |
| Soft | 49152 | Every module except the grandfathered files below. A new file over this fails the build. |

Grandfathered until a later split that existing tests cover: `verify.py`, `parameter_scan.py`, `equation_selection.py`. They may exceed 49152 bytes and must stay under 65536. `algebra_count.py` is under the soft limit and is not grandfathered; do not push it over 49152 without splitting it. The old `reference_verify.py` monolith was removed so this gate can pass. Do not put it back together.

## Gates

`.github/workflows/ci.yml` runs a fast `gates` job on pull requests and pushes to `main`, in addition to pytest and `verify --quick`. That job checks UTF-8 text, the size limits above, version mentions in `docs/STATUS.md` and `references/lean/README.md`, the false honesty flags in `benches/reference_alignment.json` and `artifacts/step15-verification.json` through `artifacts/step21-verification.json`, and `ruff check` (`E9`, `F`). `ruff format --check` covers `reference_verify/` and `.github/workflows`. The existing `check_verify_report.py` validators are unchanged in what they already required, and the reference-report check also rejects a row whose honesty flag is true.
