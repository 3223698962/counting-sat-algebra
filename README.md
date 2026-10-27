# counting-sat-algebra

以 `references/` 中的计数算法为对照，逐步实现确定性的精确 #SAT 计数，并核验时间与空间结论。当前版本为 **0.17.0**。

当前已实现参考文献 Sections 3、4 的若干部件：精确条件期望日程、线性代数设计检查、共享测量求和、Taylor 系数恢复、完整参数下的辅助算子条目查询和有限参数的稀疏放置。DIMACS 计数入口通过独立连通分量调用代数计数后端，局部后端为 `--backend common` 或 `--backend tensor`。Step 14 按实际替换多项式的支持集分解，已从同一 16384 变量全局放置恢复精确计数 1。Step 15 在连通支持上做有限 Taylor 张量消元，中间表必须落在显式预算内。它不是 Appendix D 的一般模收缩，也不给出该渐近界。Step 16 按 Lemma 4.8 计算参数日程，并在预算内运行一个通过 Lemma 4.7 全部数值前提的有限元组。该闭式在这些 n 上于形成 D 之前失败；有限元组的运行不是渐近时间或空间界的证明。一般 k-CNF 的整条管线与渐近复杂度结论仍未完成。Theorem 1.1 没有实现，本仓库不否定 ETH，也没有一般的 \(2^{o(n)}\) 算法。`theorem_1_1_complete` 与 `lean_proof_complete` 均为 false。

## 运行

```bash
python -m pip install -e ".[dev]"
pytest -q
python -m counting_sat_algebra.verify --reference --output benches/reference_alignment.json
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra count examples/count_demo.cnf --output benches/count_demo.json
```

`--reference` 运行 Steps 9–16，共 31 项检查：1024 变量实例的日程和设计检查、精确尾界检查、四个小实例的共享测量计数、结构化矩阵查询、16384 变量的稀疏放置、实际 CNF 计数、同一全局放置在启用和关闭复用时的支持集计数与预算拒绝、Step 15 的连通路径计数、完整 Taylor 系数对照和超宽因子在建表前拒绝，以及 Step 16 的 Lemma 4.8 闭式、较小 n 的前提失败、被选择的有限元组、稀疏放置、三个后端的计数或拒绝，和稠密满矩阵仍在分配前拒绝。`--step9` 至 `--step16` 可单独运行，也支持 `--output` 保存实际执行的参考检查。失败报告保留已执行的结果并标记失败步骤。

`--quick` 运行 Steps 0–16 的较小网格，并包含 Step 16 的同一组有限检查。Steps 2–8 的单独命令继续可用。`--step8` 保留原布尔枚举预算模式；默认日程 API 使用精确线性代数。

核验结果中的 `COUNTED` 表示确实恢复了计数，并与独立暴力计数或已知闭式核对；`CHECKED` 表示设计、参数或条目核对；`REFUSED` 表示预算或前提导致停止。大实例的设计通过不会被当成计数完成。

## 新的计数接口

直接读取 DIMACS：

```bash
counting-sat-algebra count examples/count_demo.cnf --json
python -m counting_sat_algebra count INPUT.cnf --no-memo
python -m counting_sat_algebra count examples/connected_path.cnf --backend tensor --output benches/connected_path.json
```

示例返回 **196**：两个独立三文字子句各有 7 个满足赋值，另有两个声明的自由变量。解析器保留声明变量数，支持注释、跨行子句和空子句。计数器去除重言式和重复文字/子句，为每个连通分量重新稀疏放置，以共享测量求和恢复局部计数，再乘以自由变量因子。相同编号化分量可复用已经完成的计数；`--no-memo` 可关闭复用。

默认每个分量最多 20 个变量，每次代数计数最多扫描 500000 个局部布尔向量，总扫描预算为 5000000；可用 `--max-component-variables`、`--vector-budget`、`--total-vector-budget` 调整。`--backend common` 是默认的共享测量后端。`--backend tensor` 改用 Step 15 的有限 Taylor 收缩，并可用 `--max-tensor-width`（默认 12）和 `--max-tensor-degree`（默认 64）限制消除桶与分量次数。预算拒绝返回退出码 2，计数和可满足性均为空。大整数用精确十六进制输出，并标记 `count_encoding`。`examples/connected_path.cnf` 是 13 个变量的连通路径，`--backend tensor` 的计数为 610。

已放置公式仍可直接调用：

```python
from counting_sat_algebra import count_placed_formula, step3_examples

placed = step3_examples()[1].placed
result = count_placed_formula(placed)
assert result.counted
assert result.count == 1
print(result.plan.boolean_vector_bound, result.boolean_vectors)
```

该函数先检查恢复和核条件，再执行 Lemma 3.2 的局部扫描与 Lemma 3.3 的固定替换，最后用 Hermite 插值恢复 H(0)。它返回操作预算、实际局部扫描次数和系数；预算不足时不启动枚举。

它不存可达测量表，也不调用全空间计数作为计算后端。测量元组、局部组和替换对仍可能很多，是否得到次指数界取决于完整参数选择、输入约化和剩余管线。

当替换后的因子支持集可以分解时，可对同一个 `placed` 调用 `count_placed_formula_by_support(placed, memo_entries=0)`。该入口检查每个替换多项式的全部非零系数，按实际支持集形成分量，在 T=1 计算足够阶数的局部 Taylor 系数，再恢复 H(0)。默认每个支持分量最多 12 个变量、次数上界 64，总局部布尔向量预算为 500000；预算按关闭复用计算。`memo_entries=0` 关闭复用，默认缓存上限为 4096。

连通支持可调用 `count_placed_formula_by_contraction(placed)`。它先为全部分量做符号计划，再在 \(\mathbb Q[z]\) 上按最小度消去布尔变量；表是 T=1 处覆盖分量次数上界的 Taylor 喷射，最后恢复 H(0)。默认消除桶最多 12 个变量。计划超出预算时不建数值表，也不返回部分计数。该入口不使用复用。它只处理中间表落在预算内的有限实例，不是 Appendix D 的一般模收缩。

## 对齐结果

| 检查 | 结果 |
| --- | --- |
| 1024 变量，r=2、Q=2、R=1 | 条件期望日程、重数、组大小、核条件、T=0 恢复通过 |
| 同一参数的 Lemma 4.7 尾界 | 7/16 > 1/4，前提不成立 |
| 16384 变量，r=2、Q=3、R=2、τ=1/2 | 全部数值前提、稀疏设计、恢复和系数检查通过；支持集计数精确为 1 |
| 同一全局放置，启用 / 关闭复用 | 实际扫描 96 / 49152 个局部布尔向量，两种模式均计数为 1 |
| 四个 Step 3 小实例 | 共享测量恢复 #F 为 1、1、7、7；Taylor 系数与全空间对照一致 |
| 结构化 M0 | 小矩阵所有条目一致；R=32、Q=192、D=2^384 可作有限条目查询 |
| DIMACS 示例，8 个声明变量 | 精确计数 196，两个连通分量，两个自由变量 |
| 10923 个独立三文字子句，32769 变量 | 精确计数 7^10923；1 次局部代数计数、10922 次复用 |
| 连通路径，13 个变量 | Step 14 支持集路线拒绝；Step 15 计数 **610**，桶宽 2，局部布尔向量 98。`COUNTED` |
| 连通路径，1024 个变量加 2 个自由变量 | 计数 **4·Fibonacci(1026)**，桶宽 2，局部布尔向量 8186，峰值表项 4094；`--backend tensor` 得到同一计数。`COUNTED` |
| 耦合替换的完整 Taylor 喷射 | 计数 1；系数 7/4、1/2、−1/4 与独立全空间 Hasse 喷射一致。计数为 `COUNTED`，系数对照为 `CHECKED` |
| 宽度 13 的单个因子 | `REFUSED`：因子支持 13 超过预算 12，未建数值表，扫描为零 |
| Lemma 4.8，16384 个单元子句 | r=0、Q=0、R=16，D 未形成。`CHECKED` |
| 8192 个单元子句 | 没有元组通过 Lemma 4.7。`CHECKED` |
| 选择器在 16384 个单元子句上 | r=2、Q=3、R=2、τ=5/16，尾概率 5/32，测量允许量 10240。不满足 R≥32 或 6R≤Q。`CHECKED` |
| 同一选择的稀疏放置 | 多项式条目 59392，测量条目 14336，稠密位置 536870912。`CHECKED` |
| 同一放置的共享测量 | `REFUSED`：q=5461、J=5462，扫描为零 |
| 同一放置的支持集复用 | 计数 **1**，扫描 **96**。`COUNTED` |
| 同一放置的有限 Taylor 消元 | 计数 **1**，扫描 **176128**，峰值表项 **72**。`COUNTED`。不是 Appendix D |

16384 变量放置仅存 **59392** 个多项式条目和 **14336** 个测量条目，对应稠密多项式位置数为 **536870912**。手选点使用 τ=1/2；Step 16 的选择器在同一 (r,Q,R) 上改用满足精确尾界的最小 τ=5/16，测量允许量从 16384 降到 10240，独立测量数 4096 仍落在允许量内。共享测量路线的计划为 q=5461、J=5462，超过扫描预算，实际扫描为零。Step 14 与 Step 16 的支持集路线都得到 **10240** 个实际支持分量；启用复用完成 13 次局部求值、10227 次复用，实际扫描 96。Step 16 的张量后端在同一放置上扫描 176128 个局部布尔向量并得到同一计数 1。32769 变量结果来自独立 CNF 分量重新放置；该实例通过数值前提，但未执行全局放置。这些独立分量计数和不相交支持集计数不是一般连通实例的时间界。Lemma 4.8 的闭式在这些 n 上没有形成 D，因此这次运行不是渐近界。

214 项 pytest 与 31 项参考检查全部通过。完整参考报告保存在 `benches/reference_alignment.json`，版本为 0.17.0，步骤为 9–16，只记录计数、检查结果和操作统计。`artifacts/step15-verification.json` 单独保存 Step 15 的四行结果。`artifacts/step16-verification.json` 单独保存 Step 16 的十行结果。

[对齐说明](docs/REFERENCE_ALIGNMENT.md) 区分 Theorem 1.1 的较强时间界、Sections 3–6 的较弱时间界和未完成的 Lean 主定理。当前测试不能作为 ETH 结论、一般 \(2^{o(n)}\) 算法或完整形式化证明。

## 文件

- [ROADMAP.md](ROADMAP.md)：已完成步骤与剩余管线。
- [docs/STATUS.md](docs/STATUS.md)：中文现状。
- [docs/PROGRESS.md](docs/PROGRESS.md)：原始检查记录与新步骤结果。
- [docs/REFERENCE_ALIGNMENT.md](docs/REFERENCE_ALIGNMENT.md)：与参考定理逐项对齐。
- `src/counting_sat_algebra/`：精确有理数实现。
- `tests/`：独立对照及预算检查。
- `references/`：本地预印本与 Lean 脚手架。

核心算术使用 Python 标准库的 `fractions.Fraction`。稠密矩阵切片上限仍为 256；结构化接口提供有限坐标查询，不分配整个巨大矩阵。
