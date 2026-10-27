# 现状

本文件是当前 `main` 的短说明。逐项勾选见 [ROADMAP.md](../ROADMAP.md)，实例、计数和拒绝记录见 [docs/PROGRESS.md](PROGRESS.md)。正文用中文；函数名、定理编号和命令保持与代码相同的英文。

This file is a short orientation snapshot. The checklists and the measured tables stay in `ROADMAP.md` and `docs/PROGRESS.md`.

## 这个仓库是什么

`counting-sat-algebra` 是预印本里计数算法若干片段的机器受限核对原型。成功标准有两条：在小到可以暴力枚举的实例上，结果与暴力计数一致；同时把代价如实记下来。

反驳指数时间假设（ETH）不在成功标准里。

## 诚实边界

本仓库不实现预印本 Theorem 1.1，不得出 “ETH 为假”，也不给出 \(2^{o(n)}\) 或次指数计数。

墙钟时间是一台机器上的一次读数。拒绝（refusal）没有计时。可复现的记录是操作次数：\(2^n\)、测量元组数、局部求值次数、赋值预算。

## Step 0–8（已完成）

一步算完成，是指小实例上与暴力计数或全空间和 \(H(t)\) 一致，并且代价与该步自己的分析相符。下面每行只标这一步核对了什么。

0. 子句多项式与 Lagrange 插值：\(T=0\) 恢复原变量，并得到 \(H(0)=\#F\)。每个样本仍枚举 \(2^n\) 个赋值。
1. 出现次数至多 3 的 3-CNF 试验台，以及次数上界 \(D=m\cdot k\cdot R\) 与 \(D=3\cdot n\cdot R\)。
2. 简化的 \((r,\tau,D)\)-设计与分块求和。手写替换族上，分块和等于全空间 \(H(t)\)。popcount 测量元组数是 \((r+1)^b\)，对块数指数。
3. 在边长至多 64 的四组显式放置上，构造 \(M_0 = V_0 f_I(W_0)\) 的一小片，并核对 \(T=0\) 恢复、次数上界和设计条件。
4. 用 Hasse 系数做 Hermite 插值，从比 Lagrange 更少的相异节点恢复 \(H(0)=\#F\)。条件数仍是 \(D+1\)。每个节点仍走一遍立方体。
5. 蕴涵环复制，使每个变量出现至多 3 次，并在暴力窗口内核对 \(\#F'=\#F\)。\(n'\) 可以大于 \(n\)；\(n'>20\) 只记录、不枚举约化后的立方体。
6. 在一台机器能跑完的固定网格上，把暴力计数与 popcount 分块和对照，其余步骤记为一致或拒绝。
7. 切片上限提到 256。边长 81 与 256 的两处显式放置通过与 Step 3 相同的恒等检查。平衡不等式不成立时，`schedule_placement` 拒绝，不选日程。
8. 不等式成立且切片未超上限时，`schedule_placement` 用逐符号选择选日程。被接纳的玩具实例先满足重数与组大小上界，再因 `design_problems_at` 的布尔点预算而拒绝。

## 当前 tip

`main` 在 PR #8 合入之后，包含 Step 8 的 balancer。包版本见 `pyproject.toml`（0.9.0）。

- 切片上限是 256。\(r=2\)、\(Q=4\) 的边长 256 在上限内。接下来的纯幂 729 与 1024 在存矩阵之前被拒绝。
- `schedule_placement` 做逐符号选择（successive choice）。字母表是 \(\{0,\ldots,r-1\}^Q\)。每个变量、每条子句都试遍当前字母，保留使二次负载势最小的符号；平局取字典序第一个。不抽样，也不枚举日程立方体 \((r^Q)^{n+m}\)。
- 被接纳的玩具实例满足重数与组大小上界，然后拒绝 `design_problems_at` 的布尔点预算。64 条互不相交的单元子句，\(r=2\)、\(Q=1\)、\(D=4\)：重数 \(16 \le 2M/D = 32\)，最大组大小 32，且 \(32\cdot r \le 2n\)。两条子句符号都被用到，因此组外也有 32 个变量，对应 \(2^{32}\) 个布尔点，预算是 4096。`place_formula` 没有被调用。1024 条单元子句在 \(r=2\)、\(Q=2\)、\(R=1\) 上是同一类拒绝，并且 \(1 \le R < Q\) 成立。
- 完整 \(M_0\) 在分配之前被拒绝：参数是 \(R \ge 32\)、\(Q \ge 192\)，\(r=2\) 时边长 \(2^{384}\)。`build_auxiliary_slice(2, 192, 32)` 与 `schedule_placement` 在这一参数上都拒绝，并且不计时。把 `max_dimension` 调到 256 以上也不分配这张矩阵。

## 仍开放的距离

下面五条来自 [docs/PROGRESS.md](PROGRESS.md) 的 “Distance to a full counting algorithm”。Step 7 和 Step 8 没有把它们划掉。

- **完整 \(M_0\)。** \(R \ge 32\) 且 \(6R \le Q\)。最小的这样一对整数是 \(R=32\)、\(Q=192\)。\(r=2\) 时边长是 \(2^{384}\)。`build_auxiliary_slice` 在分配矩阵之前拒绝这次调用。
- **能把已被接纳的公式放进设计预算的 balancer。** 现有搜索在被接纳的行上达到 Lemma 4.5 的重数与组大小上界，然后在 `design_problems_at` 之前停下。它不是 Lemma 4.5 的证明。二项式尾界没有检查。
- **保持 \(n' = n + o(n)\) 的稀疏化。** 出现次数上界 \(B=3\) 时，出现 \(d>3\) 次的变量被换成 \(d\) 个变量。扫描记下 `random_3cnf(5, 12, seed 0)` 的 \(n'=36\)，并没有枚举那个约化立方体。
- **对块数不是指数的测量空间。** popcount 设计枚举 \((r+1)^b\) 个元组。块大小为 2 时这是 \(3^{n/2}\)。本仓库跑过的若干行上，局部求值次数大于立方体大小。
- **不走立方体的 Hermite 求值。** 配对格式把 \(D+1\) 个条件装进更少的相异 \(t\)，每个 \(t\) 仍然对 \(\{0,1\}^n\) 求和。阶为 \(D\) 的完整 jet 在每个赋值上做二次量级的算术。

同一节还记下：跨块子句不被分块和分解，复制也不产生 Step 3 的放置；Theorem 1.1 未实现。把这些写下来，并不等于补上它们。

## 如何运行

```bash
python -m pip install -e ".[dev]"
pytest -q
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra.verify --step2
python -m counting_sat_algebra.verify --step3
python -m counting_sat_algebra.verify --step4
python -m counting_sat_algebra.verify --step5
python -m counting_sat_algebra.verify --step6
python -m counting_sat_algebra.verify --step7
python -m counting_sat_algebra.verify --step8
```

`python -m counting_sat_algebra.verify` 不带 `--quick` 时会多抽随机 3-CNF（\(n\) 到 10），把 Step 2 的组扫描再加一块，跑完整的 Step 6 网格，并跑边长 256 的恒等式。那条路径是精确的，也更慢，因为每个 Step 0 样本仍然花费 \(2^n\)。

`--quick` 包含 Step 3、4、5、较小的 Step 6 网格、便宜的 Step 7 检查（边长 81 和各项拒绝）以及 Step 8，并跳过边长 256。`--step7` 额外跑边长 256。`--step8` 记录被接纳实例上的拒绝、完整 \(M_0\) 的拒绝，以及对八条单元子句的一次直接设计检查；那八条不满足平衡不等式，`schedule_placement` 不会去放置它们。

退出码 0 且出现 `ALL TESTS PASSED`，表示该次运行里核对过的项与暴力一致，或者就是扫描所记录的拒绝。\(n'>20\) 的 Step 5 行记为未检查，并仍算在这一行之下。Step 6、Step 7、Step 8 的拒绝同样是记录下来的停止，不计时。Step 8 被接纳的行不是一次放置。

## 参考文献

离线副本在仓库顶层 [`references/`](../references/README.md)，不参与 `pytest`，也不参与 `verify`。

- [`references/paper/`](../references/paper/README.md)：预印本 *The Exponential Time Hypothesis is False*（OpenAI，2026-10-17，修订于 2026-10-25）的 PDF 与纯文本摘录。这是供离线查阅的本地副本。本仓库不实现该文的 Theorem 1.1。
- [`references/lean/`](../references/lean/README.md)：另一项 Lean 形式化的范围说明 `378.md`，以及脚手架 `SubexponentialSAT.lean`。主定理 `exact_count_and_sat` 以 `sorry` 结束。`sorry` 是未完成证明的占位。本 Python 仓库不是那个 Lean 项目。
