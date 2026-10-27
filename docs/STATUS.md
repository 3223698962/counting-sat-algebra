# 现状

目标是实现 `references/` 中的确定性精确计数算法，并达到可核验的时间、空间结论。当前版本 **0.17.0** 已推进到 Step 16；完整定理仍未完成。Theorem 1.1 没有实现，本仓库不否定 ETH，也没有一般的 \(2^{o(n)}\) 算法。`theorem_1_1_complete` 与 `lean_proof_complete` 均为 false。

## 本次新增

- **Step 9：** 按 Lemma 4.5 精确计算条件期望，覆盖共享变量与共享因子的相关性。每次选择、初始势及最终势均以 Fraction 核对。设计检查使用缓存的精确行空间，不再默认枚举布尔组。1024 变量的实例完成放置和恢复检查。
- **Step 10：** 实现 Lemma 3.2 的流式共享测量局部求和、Lemma 3.3 的固定替换和式 (3.13) 的 Taylor 系数恢复，再做 Hermite 插值。四个原小实例的计数及所有请求系数与独立对照一致。
- **Step 11：** 通过交换轨道的谱系数实现结构化 M0 和替换条目查询。小切片的所有条目一致；完整参数 R=32、Q=192 的有限条目可查询。
- **Step 12：** 将结构化条目查询接入稀疏放置、测量构造和精确核检查。16384 个独立单元子句在 r=2、Q=3、R=2、τ=1/2 上通过全部 Lemma 4.7 数值前提及设计检查；多项式条目 59392、测量条目 14336，稠密多项式位置数 536870912。共享测量求和在枚举前因预算拒绝，扫描为零。
- **Step 13：** 新增 DIMACS 输入和实际 CNF 计数命令。独立连通分量各自重新放置，以代数后端求局部计数，再乘以声明自由变量因子。8 变量示例得到 196；32769 变量的 10923 个独立三文字子句得到 7^10923，执行 1 次局部计数、10922 次复用。
- **Step 14：** 按实际替换多项式的支持集分解同一全局放置，以局部 Taylor 系数恢复计数。16384 变量实例得到 10240 个支持分量，最大 4 个变量、次数上界 4。启用复用与关闭复用均得到精确计数 **1**；实际扫描分别为 **96**、**49152**，后者恰好达到无复用计划上界。复用模式有 13 次局部求值、10227 次命中；无复用模式有 10240 次局部求值。超过支持分量预算时仍在枚举前拒绝。
- **Step 15：** 在同一替换多项式上做有限 Taylor 张量消元。表是 T=1 处覆盖分量次数上界的喷射，按最小度消去布尔变量；全部分量先做符号计划，超预算则不建数值表、不返回部分计数，也不使用复用。13 变量连通路径在 Step 14 被拒绝，Step 15 计数为 **610**（`COUNTED`），桶宽 2，扫描 **98**。1024 变量连通路径加 2 个自由变量计数为 **4·Fibonacci(1026)**，桶宽 2，扫描 **8186**，峰值表项 **4094**；DIMACS `--backend tensor` 得到同一计数。耦合替换的喷射 **7/4、1/2、−1/4** 与独立全空间 Hasse 喷射一致。宽度 13 的因子为 `REFUSED`，扫描为零。这不是 Appendix D 的一般模收缩，也不给出其渐近界。
- **Step 16：** 按 Lemma 4.8 的证明顺序计算整数日程。16384 个交替单元子句得到 r=0、Q=0、R=16，在形成 D 之前停止（`CHECKED`）。λ=2¹⁹ 的整数为 r=2、Q=215、R=2900，因 Q≤R 同样不形成 D。预算内的搜索只接受通过全部 Lemma 4.7 数值前提的元组，τ 取满足精确尾界的最小值。8192 个单元子句没有这样的元组（`CHECKED`）。16384 个单元子句选择 r=2、Q=3、R=2、τ=5/16，尾概率 5/32，测量允许量 10240；这比手选的 τ=1/2 更紧，但不满足 R≥32 或 6R≤Q。稀疏放置复现 59392 个多项式条目、14336 个测量条目和 536870912 个稠密位置。共享测量在枚举前拒绝，扫描为零（`REFUSED`）。支持集复用计数为 **1**，扫描 **96**。有限 Taylor 消元计数为 **1**，扫描 **176128**，峰值表项 **72**。试验预算为 1 时在日程搜索前拒绝。稠密边长上限不变，R=32、Q=192 仍在分配前拒绝。选择并运行这些有限元组不是渐近时间或空间界的证明。
- **回归：** 214 项 pytest 全部通过，`--reference` 的 31 项检查全部通过。保存报告为 `benches/reference_alignment.json`（版本 0.17.0，Steps 9–16）、`benches/count_demo.json`、`benches/connected_path.json`、`artifacts/step15-verification.json` 和 `artifacts/step16-verification.json`。参考报告支持命令直接重建并记录失败步骤，且保持两个完成标志为 false。

## 运行

```bash
python -m pip install -e ".[dev]"
pytest -q
python -m counting_sat_algebra.verify --reference --output benches/reference_alignment.json
python -m counting_sat_algebra.verify --step14
python -m counting_sat_algebra.verify --step15
python -m counting_sat_algebra.verify --step16
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra count examples/count_demo.cnf --output benches/count_demo.json
python -m counting_sat_algebra count examples/connected_path.cnf --backend tensor --output benches/connected_path.json
```

`--step9` 至 `--step16` 分别检查新步骤，也支持 `--output` 保存参考结果。`--quick` 包含 Step 16。`--step8` 保留原布尔预算拒绝路径；日程 API 默认采用线性代数检查。实际 DIMACS 计数命令支持 `--json`、`--no-memo`、`--backend common|tensor` 和显式分量/枚举预算；张量后端另有 `--max-tensor-width` 与 `--max-tensor-degree`。已放置公式通过 `count_placed_formula_by_support` 或 `count_placed_formula_by_contraction` 计数。`memo_entries=0` 关闭支持分量复用；张量收缩本身不使用复用。拒绝时不返回部分计数。

## 尚未完成

Lemma 4.8 的闭式已经可以计算，但在当前稀疏放置预算能覆盖的 n 上，它在形成 D 之前失败。n=16384 上通过 Lemma 4.7 的有限成员可以放置，并由既有后端计数或在枚举前拒绝；这不证明渐近时间或空间界，也不满足 R≥32 和 6R≤Q。一般输入仍不会因此自动得到该界。Sections 5–6 的输入约化与方程构造、Appendix D 的一般模收缩、完整位复杂度和空间分析，以及 Lean 主定理的实际证明都还未完成。分离的替换支持集可完成这一放置的计数；桶宽不超过预算的连通支持也可由 Step 15 计数。更宽的因子或消除桶仍会在建表前拒绝。独立分量计数、不相交支持集计数和这些有限收缩都不是一般连通实例的时间界。

参考文献 Sections 3–6 与 Lean 范围说明的时间界比 Theorem 1.1 的 Appendix D 时间界弱。要达到用户要求的参考结论，必须区分并完成对应路线，不能把有限测试当成渐近结论。

详见 [REFERENCE_ALIGNMENT.md](REFERENCE_ALIGNMENT.md)、[ROADMAP.md](../ROADMAP.md) 和 [PROGRESS.md](PROGRESS.md)。
