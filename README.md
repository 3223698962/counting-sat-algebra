# counting-sat-algebra

根据 `references/` 中的论文和 Lean 程序逐步实现确定性的精确 #SAT 计数。当前版本为 **0.30.0**，主线已转向新加入的 `SubexponentialSAT/SqrtLogSpace`：保留声明变量数，流式遍历互不重叠的稀疏化叶，叶内先装配日程稀疏码，再用违反质量贪心选择历史行，最后单独执行接受查询并累加完整的自然数计数；hash 叶可选地把标记与接受查询路由到流式高变量秩循环，计数扩张环、Hermite 放置、修正对分解与 Taylor 读回、无支持表的流式逐点叶、以及记录端口的叶装配与放置已有有限见证。

Step 23 实现 Lean 二进制语法、精确二进参数搜索、规范化和深度优先计数前沿。Step 24 新增模 `3^t` 的 Newton 指示系数、共享标记截断、布尔平方自由多项式求和，以及真实的接受/标记违反计数查询。Step 25 按日程整数构造三元稀疏码并逐行再生，用真实标记查询累加全部子句的违反质量做贪心选择，核对零质量后单独执行接受查询。Step 26 计算早期高出现阈值与高变量前端，并可按秩流式路由全部查询。Step 27 构造模 `3^t` 的计数扩张环、固定互异节点与 `Q*J > degree` 的 Hermite 条件。Step 28 流式枚举修正对 `(E,U)` 并用布尔有限差分恒等式精确求和，再经逐节点 J 系数流与 Newton–Hermite 读回计数。Step 29 增加无平方自由支持表的逐点修正对叶（`--leaf-backend stream`）：每个修正对只做一次点求值并精确累加。Step 30 把环常量、节点表、插值缓存与成帧行流按记录的 `QueryReady` 端口放置，并用一条贯通运行从日程走到计数字。`--backend lean --leaf-backend native|hash|stream` 已接通这些部件。完整的 Lean 叶程序移植仍需物理机器、分组替换与整机资源。

新的 Lean 参考源码含精确计数和全执行前缀的多项式空间定理。静态审计覆盖 **3969 个本地模块**，没有缺失的本地导入、代码级 `sorry`/`admit` 或 `axiom` 声明。审计不等于 Lean 内核验证；本机尚未编译该库。旧的单文件脚手架另行保留，详见 [Lean 参考说明](references/lean/README.md)。

本包尚未建立一般输入的渐近时间与空间界。`theorem_1_1_complete`、`lean_proof_complete`、`recursive_contraction` 仍为 false；这些收据不把新参考库称为未完成的旧脚手架。本项目目前没有一般的 `2^{o(n)}` 实现，也不据此宣称否定 ETH。

## 运行

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra.verify --reference --output benches/reference_alignment.json
python -m counting_sat_algebra.verify --step24 --output artifacts/step24-verification.json
python -m counting_sat_algebra.verify --step25 --output artifacts/step25-verification.json
python -m counting_sat_algebra.verify --step26 --output artifacts/step26-verification.json
python -m counting_sat_algebra.verify --step27 --output artifacts/step27-verification.json
python -m counting_sat_algebra.verify --step28 --output artifacts/step28-verification.json
python -m counting_sat_algebra.verify --step29 --output artifacts/step29-verification.json
python -m counting_sat_algebra.verify --step30 --output artifacts/step30-verification.json
python -m counting_sat_algebra.lean_reference --output artifacts/lean-reference-audit.json
```

`--reference` 运行 Steps 9–24；`--quick` 同时运行原 Steps 0–8 的较小网格。`--step9` 至 `--step24` 可单独运行并保存报告。报告只记录实际执行的检查，失败保留已执行结果。

## DIMACS 计数

```bash
python -m counting_sat_algebra count examples/count_demo.cnf --backend lean --leaf-backend native --json
python -m counting_sat_algebra count examples/connected_path.cnf --backend tensor --json
python -m counting_sat_algebra count INPUT.cnf --backend lean --leaf-backend native --sparsify-M 2 --output benches/input-count.json
```

第一个示例精确计数 **196**，第二个示例精确计数 **610**。解析器保留声明变量、跨行子句和空子句；计数结果包括未使用变量的全部赋值。`--sparsify-M` 是显式有限参数；省略时按 Lean 的整数条件选择第一个通过测试的二的幂。

| 入口 | 当前计算方式 |
| --- | --- |
| `--backend common`，默认 | 按原 CNF 的独立分量重新放置，进行共享测量和 Hermite 恢复 |
| `--backend tensor` | 按独立分量运行预算内的有限 Taylor 消元 |
| `--backend lean`，叶选择 `common` 或 `tensor` | 流式稀疏化前沿，调用上述有限叶后端 |
| `--backend lean --leaf-backend native` | 同一前沿，使用新构造的三元指示系数和平方自由代数计数叶 |
| `--backend lean --leaf-backend hash` | 同一前沿，先装配日程三元稀疏码，再用违反质量贪心选择历史行并单独计数接受查询；`--high-loop` 可把查询路由到高变量秩循环 |
| `--backend lean --leaf-backend stream` | 同一前沿，用逐点修正对流直接求值，不构造平方自由支持表 |

原后端预算继续由 `--max-component-variables`、`--vector-budget`、`--total-vector-budget`、`--max-tensor-width` 和 `--max-tensor-degree` 控制。新后端使用 `--max-ring-precision`、`--max-algebra-terms`、`--max-algebra-products`、`--max-algebra-operations` 和 `--max-algebra-bits`；乘法与标量操作预算由整棵前沿共享。`hash` 叶另受 `--hash-depth`、`--max-hash-candidates` 与逐叶标记查询预算限制，Lean 日程的候选深度 `distanceDenom` 超出 Python 预算时直接拒绝；`--high-loop` 把每次标记/接受查询展开为高变量秩枚举，秩数受 `HighBudgets` 限制。当前平方自由支持表仍可能指数增长。

`COUNTED` 表示已完成精确计数，`CHECKED` 表示核对了有限对象，`REFUSED` 表示预算阻止完成。拒绝的计数与可满足性均为 null，退出码为 2；已经做过的工作仍可见。大计数以精确十六进制输出，并标记 `count_encoding`。

## 新的查询 API

```python
from counting_sat_algebra import CNF, TernaryHashRow, count_ternary_query

F = CNF(5, [[(v, False)] for v in range(3)])
history = (TernaryHashRow(((0, 1), (1, 1), (2, 1))),)
accepted = count_ternary_query(F, history)
marked = count_ternary_query(F, history, marked_clause=0)
assert accepted.count == 8
assert marked.count == 4
```

接受查询计算历史哈希方程的公共零点。标记查询再要求某个原始子句被违反，用于后续贪心选行的精确评分。它们和原 CNF 的 #SAT 有不同语义；`count_cnf_ternary(F)` 以逐子句的单位哈希行计数原公式。`count_cnf_greedy(F)` 装配日程稀疏码并用上述标记查询做贪心选择，零质量后单独执行接受查询得到原 CNF 的精确计数。`count_cnf_high_loop(F, rows)` 用同一批真实查询逐秩枚举高变量赋值并精确累加，秩和恒等于直接查询计数。`extension_schedule(ambient, ...)` 复现日程整数并核对 `Q*J > degree`；`counting_extension(t, d)` 搜索不可约单首提升并给出固定节点与唯一自然数读回。`streamed_boolean_sum(terms, n, J)` 按修正对流式求和；`streamed_taylor_count(...)` 经逐节点 J 系数流与 Hermite 读回同一计数。定义、预算和数学恒等式见 [LEAN_NATIVE.md](docs/LEAN_NATIVE.md)。

## 项目文件

- [ROADMAP.md](ROADMAP.md)：按新 Lean 依赖安排的剩余主线。
- [docs/STATUS.md](docs/STATUS.md)：当前状态和可运行入口。
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)：核心、后端、见证与模块大小门禁。
- [docs/REFERENCE_ALIGNMENT.md](docs/REFERENCE_ALIGNMENT.md)：区分新 Lean 的时间界、旧有理数路线和 Appendix D。
- [docs/LEAN_SPEC_CHECKLIST.md](docs/LEAN_SPEC_CHECKLIST.md)：Lean 声明到 Python 的有限义务。
- [docs/steps/completed.md](docs/steps/completed.md)：Steps 0–24 已完成清单。
- [docs/PROGRESS.md](docs/PROGRESS.md)：历史结果与各步骤的操作统计。
- [docs/APPENDIX_D.md](docs/APPENDIX_D.md)、[PACKET_PACKING.md](docs/PACKET_PACKING.md)、[WEIGHTED_MODES.md](docs/WEIGHTED_MODES.md)：已实现的有限 Appendix D 部件。
- [docs/APPLICATION_FACTORING.md](docs/APPLICATION_FACTORING.md)：因子分解应用说明、引擎归属与结论距离。
- [docs/RESOURCES_AND_CONCLUSION.md](docs/RESOURCES_AND_CONCLUSION.md)：资源收据、三条路线的结论可达性与开放义务。

CI 核对 UTF-8、模块大小、版本、有限结果标志和报告完整性，并运行 Ruff、pytest 与 `verify --quick`。绿色 CI 表示这些检查通过，不替代 Lean 内核验证或完整 Python 复杂度分析。
