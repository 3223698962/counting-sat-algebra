# counting-sat-algebra

以 `references/` 中的计数算法为对照，逐步实现确定性的精确 #SAT 计数，并核验时间与空间结论。当前版本为 **0.14.0**。

当前已实现参考文献 Sections 3、4 的若干部件：精确条件期望日程、线性代数设计检查、共享测量求和、Taylor 系数恢复、完整参数下的辅助算子条目查询和有限参数的稀疏放置。新增 DIMACS 计数入口，通过独立连通分量调用代数计数后端。一般 k-CNF 的整条管线与渐近复杂度结论仍未完成。

## 运行

```bash
python -m pip install -e ".[dev]"
pytest -q
python -m counting_sat_algebra.verify --reference
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra count examples/count_demo.cnf --output benches/count_demo.json
```

`--reference` 运行 Steps 9–13，共 14 项检查：1024 变量实例的日程和设计检查、精确尾界检查、四个小实例的共享测量计数、结构化矩阵查询、16384 变量的稀疏放置和预算拒绝，以及实际 CNF 计数。`--step9` 至 `--step13` 可单独运行。

`--quick` 运行 Steps 0–13 的较小网格。Steps 2–8 的单独命令继续可用。`--step8` 保留原布尔枚举预算模式；默认日程 API 使用精确线性代数。

核验结果中的 `COUNTED` 表示确实恢复了计数，并与独立暴力计数或已知闭式核对；`CHECKED` 表示设计、参数或条目核对；`REFUSED` 表示预算或前提导致停止。大实例的设计通过不会被当成计数完成。

## 新的计数接口

直接读取 DIMACS：

```bash
counting-sat-algebra count examples/count_demo.cnf --json
python -m counting_sat_algebra count INPUT.cnf --no-memo
```

示例返回 **196**：两个独立三文字子句各有 7 个满足赋值，另有两个声明的自由变量。解析器保留声明变量数，支持注释、跨行子句和空子句。计数器去除重言式和重复文字/子句，为每个连通分量重新稀疏放置，以共享测量求和恢复局部计数，再乘以自由变量因子。相同编号化分量可复用已经完成的计数；`--no-memo` 可关闭复用。

默认每个分量最多 20 个变量，每次代数计数最多扫描 500000 个局部布尔向量，总扫描预算为 5000000；可用 `--max-component-variables`、`--vector-budget`、`--total-vector-budget` 调整。预算拒绝返回退出码 2，计数和可满足性均为空。大整数用精确十六进制输出，并标记 `count_encoding`。

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

## 对齐结果

| 检查 | 结果 |
| --- | --- |
| 1024 变量，r=2、Q=2、R=1 | 条件期望日程、重数、组大小、核条件、T=0 恢复通过 |
| 同一参数的 Lemma 4.7 尾界 | 7/16 > 1/4，前提不成立 |
| 16384 变量，r=2、Q=3、R=2、τ=1/2 | 全部数值前提、稀疏设计、恢复和系数检查通过；全局求和在枚举前拒绝 |
| 四个 Step 3 小实例 | 共享测量恢复 #F 为 1、1、7、7；Taylor 系数与全空间对照一致 |
| 结构化 M0 | 小矩阵所有条目一致；R=32、Q=192、D=2^384 可作有限条目查询 |
| DIMACS 示例，8 个声明变量 | 精确计数 196，两个连通分量，两个自由变量 |
| 10923 个独立三文字子句，32769 变量 | 精确计数 7^10923；1 次局部代数计数、10922 次复用 |

16384 变量放置仅存 **59392** 个多项式条目和 **14336** 个测量条目，对应稠密多项式位置数为 **536870912**。其全局计数计划为 q=5461、J=5462，超过扫描预算，实际扫描为零。32769 变量结果来自独立分量分解；该实例通过数值前提，但未执行全局放置。

[对齐说明](docs/REFERENCE_ALIGNMENT.md) 区分 Theorem 1.1 的较强时间界、Sections 3–6 的较弱时间界和未完成的 Lean 主定理。当前测试不能作为 ETH 结论或完整形式化证明。

## 文件

- [ROADMAP.md](ROADMAP.md)：已完成步骤与剩余管线。
- [docs/STATUS.md](docs/STATUS.md)：中文现状。
- [docs/PROGRESS.md](docs/PROGRESS.md)：原始检查记录与新步骤结果。
- [docs/REFERENCE_ALIGNMENT.md](docs/REFERENCE_ALIGNMENT.md)：与参考定理逐项对齐。
- `src/counting_sat_algebra/`：精确有理数实现。
- `tests/`：独立对照及预算检查。
- `references/`：本地预印本与 Lean 脚手架。

核心算术使用 Python 标准库的 `fractions.Fraction`。稠密矩阵切片上限仍为 256；结构化接口提供有限坐标查询，不分配整个巨大矩阵。
