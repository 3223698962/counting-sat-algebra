# Architecture

当前包版本 **0.28.0**。主线按新加入的 Lean `SqrtLogSpace/Final.lean` 移植精确计数程序。Step 23 提供保留声明 n 的流式互不重叠前沿，Step 24 提供原生三元代数查询和计数叶，Step 25 提供日程稀疏码与违反质量贪心选行，Step 26 提供高变量前端与流式秩循环，Step 27 提供计数扩张环与 Hermite 放置条件，Step 28 提供流式修正对与 Taylor 读回。原 Sections 3–6 与 Appendix D 的有限后端继续保留。

新 Lean 库有最终定理的证明源码，静态导入闭包审计通过；本机尚未执行内核编译。旧单文件中的 `sorry` 只描述历史脚手架。本包尚未实现一般输入的完整渐近算法与资源精化，`theorem_1_1_complete`、`lean_proof_complete` 与 `recursive_contraction` 保持 false。具体区别见 [Lean 参考说明](../references/lean/README.md)。

## Layers and step numbers

Step numbers are historical milestones. The code is easier to read as three layers.

**Core.** Formulas, exact rationals, designs, rings and Lean-compatible syntax. These modules do not claim a general runtime bound: `cnf.py`, `dimacs.py`, `replacements.py`, `linear_algebra.py`, `design.py`, `modular_ring.py`, `reference_parameters.py`, `lean_encoding.py`, `lean_parameters.py`, `lean_ring.py` and `lean_squarefree.py`.

**Backends.** Finite procedures that return a count, a checked object, or a refusal before enumeration:

| Backend | Module | What it actually does |
| --- | --- | --- |
| Shared measurement | `common_count.py` | Section 3 local sums and Taylor recovery |
| Disjoint supports | `support_count.py` | Counts one placed formula by replacement support |
| Finite Taylor elimination | `tensor_count.py` | Step 15 tables inside an explicit budget. Not Appendix D recursion |
| DIMACS components | `component_count.py`, `cli.py` | Re-places each component; `common` and `tensor` retain their finite behavior |
| Lean disjoint frontier | `lean_sparsify.py`, `cli.py` | Preflights and replays depth-first leaves on the original n; discards completed leaves and shares the arithmetic work budget |
| Native ternary queries | `lean_query.py` | Newton indicator digits, shared U truncation and squarefree sums; counts history acceptance, marked clause violations or unit-row CNF leaves |
| Sparse-code greedy leaves | `lean_code.py`, `lean_greedy.py` | Scheduled ternary code with slot-regenerated rows and the violation-mass greedy selector; the separate final accepted query counts the leaf. The Lean schedule depth exceeds the Python candidate budgets |
| High-variable loops | `lean_high.py` | Early threshold, ceiling cutoff, ambient-order mask and the streamed conditioned rank loop; the hash leaf can route every query through it |
| Counting extension rings | `lean_ext.py` | Scheduled ring parameters, the searched irreducible monic lift, fixed radix nodes with unit differences and the Q*J > degree Hermite condition with exact Taylor jets |
| Streamed correction pairs | `lean_taylor.py` | Boolean finite-difference pair stream in fold order, per-node J-coefficient blocks and the Newton-Hermite readback |
| Section 5–6 | `sparse_systems.py`, `algebra_count.py`, `equation_selection.py` | Budgeted finite algebra on a declared \(n\) |
| Appendix D leaf and recovery | `modular_network.py`, `modular_recovery.py` | Direct edge scans and tapered Hermite recovery |
| Collision layout | `modular_packing.py` | Labels, packets, and a child port topology. No tensor evaluation and no recursive call |
| Weighted modes | `modular_modes.py` | One packet's finite D.4.5 mode table, exponent budget and one crossing Cauchy residue. No packet tensor |

**Witnesses.** `reference_verify/` runs the finite checks and writes the JSON reports. One module per `run_*_suite` family, plus `protocol.py` for the COUNTED / CHECKED / REFUSED row and the false completion flags, and `report.py` for the writer. `__init__.py` still exports `REFERENCE_SUITES` and `write_reference_report`. A budget refusal does not return a partial count.

`verify.py` is the CLI orchestrator for Steps 0–24. It is intentionally still one module. Steps 17–22 retain separate APIs; Steps 23–24 are wired into DIMACS through `--backend lean` with `--leaf-backend common|tensor|native`.

## What is live, and what is a finite witness

The active path follows the supplied SqrtLogSpace counting program:

- `lean_sparsify.py` stores only the current frontier. Its finite core selection and branching match the tested Lean semantics; general time and all-prefix space transfer remain open.
- `lean_query.py` constructs the indicator itself and performs no Boolean assignment scan. Full squarefree support tables can still grow exponentially. Its present CNF leaf uses J=1; the dense t*J ring is not yet a complete streamed Taylor evaluator.
- `lean_code.py` builds the scheduled ternary code and regenerates rows; `lean_greedy.py` scores canonical candidates by accumulated marked-violation mass and counts the separate final accepted query; `lean_high.py` streams the high-incidence assignments and accumulates the conditioned queries without changing any count; `lean_ext.py` builds the counting extension ring with fixed nodes and checks the Hermite margin; `lean_taylor.py` streams the correction pairs and re-reads the per-node coefficient blocks by Newton Hermite. Switching the production leaf to the streamed path, native placement and whole-program resource accounting remain open.
- Scalar-operation receipts and logical encoding sizes do not measure the Lean primitive machine or Python heap. Those require whole-program resource accounting and refinement.

Appendix D remains a later route. `modular_modes.py` checks one finite D.4.5 packet; D.4.6 packet coefficients, D.4.7 recursive children and D.3 incidence encoding remain open. Step 21 supplies child port topology and Step 22 a frozen mode description, without cache, rows or origins. Step 15's finite Taylor elimination does not supply that recursion.

A counted formula or checked parameter tuple is not an asymptotic bound. The new Lean endpoint's exponent still contains `(log log n)^2 / sqrt(log n)`, while Appendix D has a stronger target. See [REFERENCE_ALIGNMENT.md](REFERENCE_ALIGNMENT.md), [LEAN_NATIVE.md](LEAN_NATIVE.md) and [ROADMAP.md](../ROADMAP.md).

## Module size

CI rejects a Python file under `src/counting_sat_algebra/` that exceeds these limits. Sizes are file bytes.

| Limit | Bytes | Rule |
| --- | --- | --- |
| Hard | 65536 | Every module. A file over this fails the build. |
| Soft | 49152 | Every module except the grandfathered files below. A new file over this fails the build. |

Grandfathered until a later split that existing tests cover: `verify.py`, `parameter_scan.py`, `equation_selection.py`. They may exceed 49152 bytes and must stay under 65536. `algebra_count.py` is under the soft limit and is not grandfathered; do not push it over 49152 without splitting it. The old `reference_verify.py` monolith was removed so this gate can pass. Do not put it back together.

## Gates

`.github/workflows/ci.yml` runs a fast `gates` job on pull requests and pushes to `main`, in addition to pytest and `verify --quick`. That job checks UTF-8 text, the size limits above, version mentions in `docs/STATUS.md` and `references/lean/README.md`, the false completion flags and matching versions in `benches/reference_alignment.json` and `artifacts/step15-verification.json` through `artifacts/step28-verification.json`, and `ruff check` (`E9`, `F`). `ruff format --check` covers `reference_verify/` and `.github/workflows`. `check_verify_report.py` requires the complete Steps 9–28 report and its K–AD success summary. These gates do not invoke Lean.
