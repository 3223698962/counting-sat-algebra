# Lean 原生三元计数查询

包版本 **0.27.0**。Step 24 将新 Lean 库中的三元指示系数、共享标记截断和查询语义实现为可运行的有限 Python 后端；Step 25 增加日程三元稀疏码、逐行再生和违反质量贪心选行，`--backend lean --leaf-backend hash` 用它们计数流式稀疏化叶；Step 26 计算早期阈值与高变量前端，并可按秩流式路由全部查询；Step 27 构造计数扩张环、固定节点与 `Q*J > degree` 的 Hermite 条件。当前查询直接展开平方自由支持表，完整的流式 Taylor 叶、分组替换和实际放置仍列在 [ROADMAP.md](../ROADMAP.md)。

## 源声明与接口

下表中的 Lean 路径相对于 `references/lean/math/lean/OAI/Algorithms/SubexponentialSAT/`。

| Lean 来源 | Python 接口 | 本次接通的义务 |
| --- | --- | --- |
| `SqrtLog/IdempotentCore.lean`：`newton`、`newton_error`、`iterate_idempotent` | `construct_ternary_indicator(t)` | 实际计算 Newton 迭代，不接收外部指示表 |
| `SqrtLog/IdempotentPolynomial.lean`：`normalized`、`digit`、`expansion` | `TernaryIndicator.digits`、`polynomial` | 模单首多项式归一化，再逐位除法得到 t 个三系数 digit |
| `SqrtLog/IdempotentBase.lean`：`baseCoefficients`、`actual_coefficients_indicator` | `TernaryIndicator.row`、`evaluate`、`encode` | 模 `3^t` 的零剩余类指示器；每个系数编码为 `2*t` 个二进制位 |
| `SqrtLogSpace/TruncatedRing.lean`：`Dense`、`readout` | `TruncatedRing`、`DenseElement` | U 与 Z 分别截断；完整乘积之后才读出 U 系数 |
| `SqrtLog/SquarefreeBasic.lean` | `SquarefreeAlgebra`、`SquarefreePolynomial` | 支持并集乘法和布尔求和，保持 `X_v^2=X_v` |
| `SqrtLogSpace/HashQueryCount.lean`：`payload_accepted_count`、`payload_marked_count` | `count_ternary_query` | 保留原始子句编号、文字重复和声明 n；标记表示违反某个子句 |
| 同一查询的逐子句单位行 | `count_cnf_ternary`、`count_cnf_streaming(..., leaf_backend='native')` | 返回原 CNF 或整个互不重叠前沿的完整自然数计数 |
| `SparseCodeSearch.lean`、`SparseCodeRows.lean`：`localSearch`、`graphSearch`、`serializedRows` | `construct_sparse_code`、`SparseTannerCode.row`、`certify_sparse_code` | 日程整数、两个完整搜索、槽位再生与独立证书 |
| `HashGreedy`、`HashGreedyCNF`：`mass`、`best`、`history`、`loopCap` | `select_history_rows`、`count_cnf_greedy` | 标记质量累加、严格最小保留最早、零质量后的单独接受查询 |
| `HighQueryProgram.lean`、`HighQueryCount.lean`：`highSet`、`threshold`、`rankSum` | `high_frontend`、`count_cnf_high_loop` | 早期阈值、截断与掩码；逐秩条件化查询与精确累加 |
| `RingCore.lean`、`IrreducibleSearch.lean`、`RingModelNodes.lean`：`Base`、`O`、`poly`、`passes`、`residueNode`、`difference_isUnit` | `counting_extension`、`node_element`、`ExtensionSchedule` | 日程整数、不可约单首搜索、固定节点与唯一自然数读回 |

这些是源码和有限行为的对应关系。Python 测试与静态源码审计不提供 Lean 内核编译或 Python 程序精化证明。

## Newton 指示系数

在 \(R=\mathbb Z/3^t\mathbb Z\) 中令

\[
H(Y)=Y(Y-1)(Y-2),\qquad e_0(Y)=\frac{(Y-1)(Y-2)}2.
\]

2 在 R 中可逆。初始误差 \(e_0^2-e_0=H(Y)(Y-3)/4\)。迭代

\[
e_{a+1}=3e_a^2-2e_a^3
\]

使误差中的 H 次数至少翻倍。程序执行 \(\lceil\log_2 t\rceil\) 轮，每轮模 \(H^t\) 归一化，得到次数小于 3t 的多项式，再用单首除法写成

\[
e(Y)=\sum_{s=0}^{t-1}(a_s+b_sY+c_sY^2)H(Y)^s.
\]

任意整数 y 的 H(y) 都能被 3 整除。故模 `3^t` 时，`indicator.evaluate(y)` 在 `y % 3 == 0` 时为 1，在另两个剩余类时为 0。例如 t=2 时三系数表为 `((1,3,5),(3,7,3))`。

## 历史行、标记与共享截断

对原始子句 e，违反多项式 \(q_e(X)\) 等于各文字为假的条件之积；正文字使用 \(1-X_v\)，负文字使用 \(X_v\)。它在布尔点上为 0 或 1。重复文字不改变值，重言子句的违反多项式为零，空子句的违反多项式为一。

`TernaryHashRow(((e, weight), ...))` 给出 \(L_i(X)=\sum_e weight_{i,e}q_e(X)\)，每个 weight 必须是 0、1 或 2。同一子句可以在原始输入或一行中重复出现，编号按调用时的原始子句列表解释。

每行构造

\[
B_i(X,U)=\sum_{s<t}(a_s+b_sL_i+c_sL_i^2)H(L_i)^sU^s.
\]

所有行在同一个模 \(U^t\) 的乘积中累积。最后才将保留的 U 系数求和。布尔点上，累计标记次数达到 t 的项含有至少 t 个因子 3，因此丢弃这些项不改变模 `3^t` 的指示值。

`readout` 是加法映射。一般元素的乘积不能改为先分别读出再相乘：在模 \(U^2\) 的环里，\(U\cdot U=0\)，但两次读出分别都是 1。`TruncatedRing(t,J,modulus)` 对 U 次数保留 `0..t-1`，对 Taylor 变量 Z 次数保留 `0..J-1`；两个截断不能合并成一个总次数守卫。

接受查询计数满足全部 \(L_i(X)=0\pmod3\) 的赋值。标记查询在其上再乘以 \(q_e(X)\)，计数接受集里违反子句 e 的赋值。这两个数不自动等于原 CNF 的 #SAT。`count_cnf_ternary` 流式生成每个原始子句的单位行；此时接受条件恰好是全部子句都满足。

## 平方自由求和与精确读回

支持集 S 表示单项式 \(\prod_{v\in S}X_v\)。乘法将两个支持取并集，故 \(X_v^2=X_v\)。固定坐标先代入多项式；设剩余自由坐标数为 f，则

\[
\sum_{x\in\{0,1\}^{f}}P(x)
 =\sum_S coefficient(S)\,2^{f-|S|}.
\]

原生查询使用右式，不枚举布尔赋值。非平凡查询选择最小的正整数 t 使 \(3^t>2^n\)，并返回最小非负剩余；它唯一确定范围 `0..2^n` 中的自然数计数。显式 `precision` 也必须满足此条件。空历史的未标记查询和含空子句的原 CNF 可以直接返回已知精确计数。

`fixed={v: bit}` 将原始坐标 v 条件化，不重编号、不缩减结果中的 n。计数只包含与这些固定比特一致的赋值；声明中未使用且未固定的变量仍贡献全部赋值。

```python
from counting_sat_algebra import (
    CNF, TernaryHashRow, count_cnf_ternary, count_ternary_query,
)

F = CNF(5, [[(v, False)] for v in range(3)])
history = (TernaryHashRow(((0, 1), (1, 1), (2, 1))),)
accepted = count_ternary_query(F, history)
marked = count_ternary_query(F, history, marked_clause=0)
assert (accepted.count, marked.count, count_cnf_ternary(F).count) == (8, 4, 4)
assert accepted.marker_sums == (8, 36, 45, 0)  # modulus = 81
```

## 稀疏码与贪心选行

日程整数直接来自 `SparseCodeParameters`：`logSize = floor(log2 m) + 1`、`gapDenom = 200·B·logSize`、`degreeBound = 8·gapDenom²·logSize`、`degree = min(m, degreeBound)`、`localRows = degree/(8B) + 1`、局部距离 `degree/(100·B·logSize) + 1`、`distanceDenom = 20000·B²·logSize²` 和最小立方维数 `r`（`degree·2^r ≥ m`）。每次搜索都是完整的递增里程枚举：先搜 Walsh 满足 `q·|Σ±1| ≤ D` 的移位元组，再搜局部矩阵的拒绝谓词（存在非零码字、权重低于距离）；保存首个通过的候选。独立核对包括相对距离恒等式 `(1/(100 B ℓ))·((1/(100 B ℓ)) − 1/gapDenom) = 1/distanceDenom`（精确有理数）、核下界 `distanceDenom·min_weight ≥ m`、`depth ≤ distanceDenom`，以及手工构造码的接受/拒绝证书。

行按 `serializedRows` 的逆序再生：第 i 行对应时间序 `row_count − 1 − i`，列 `= (side ? vertex ^ shifts[port] : vertex)·degree + port`，仅保留 `column < m`，系数取局部矩阵对应项。`SparseTannerCode.row` 与存储序列一致，不保存候选库。

贪心从码行历史开始。质量 `mass(H) = Σ_e #(H 接受且违反子句 e)`，正是全部标记查询计数之和；对规范候选族逐里程枚举，每个候选评分 `mass(candidate :: H)`，严格小于保留最早，然后加入胜者。机器的轮数是固定 `loopCap = 2(n + size m + 1)`；质量为零后其余轮只前置零质量行、不改变接受集，因此 Python 在首个零质量历史停止。零质量即接受集恰为原 CNF 模型集，最后单独执行一次只含历史行的接受查询。候选深度是显式有限参数；Lean 日程深度为 `distanceDenom`，超出候选预算时在分配前拒绝。

## 高变量前端与秩循环

阈值来自 `EarlyThreshold`：`L(n) = floor(log2(n + 16)) + 1`（`Asymptotic.L`）、`threshold = isqrt(L) + 1`、截断 `cutoff = ceil(totalScope · threshold / ambient)`；高变量是度数严格大于截断的坐标，掩码按原环境次序给出。独立核对包括两条整数上界 `card·(cutoff + 1) ≤ totalScope` 与 `card·threshold ≤ ambient`，以及度数恰等于截断时被排除的边界。

秩循环把高坐标排序后按二进制里程枚举全部 `2^|H|` 个赋值；每个秩执行一次真实条件化查询（高比特固定、其余坐标保留），把精确计数累加。`rankSum_fullCount` 使累计和等于原声明 n 上的无条件计数：hub 公式（8 变量、7 子句）直接计数 **129**，两秩之和相同；两 hub 公式 4 秩。`rows=None` 走逐子句单位行；空历史保留全部比特。`--high-loop` 让 hash 叶的全部标记/接受查询都经过该循环，直接与路由计数一致。秩数与代数工作均有预算；拒绝不返回部分计数。

## 计数扩张环与 Hermite 放置

日程整数按 `LeafParametersProgram.produce` 精确复现：`L = floor(log2(n+16)) + 1`、`w = floor(log2 L) + 1`、`s = isqrt(L)/w²`、`r = max(2,s)`、`layoutDen = 64·(bit_length r)`、`Q = L/layoutDen`、`R = isqrt(128·Q·E) + 4`（`E = bit_length(64·occ·r·L²)`）、`t = ambient + size(clauses) + 1`、`d = size(Q) + 1`、`N = r^(2Q)`、`J = width·R·(2·factors + 3(t−1) + 1)/Q + 1`。守卫在 `Q ≥ 2`、`4R ≤ Q`、`16·Q·r² ≤ freeCount` 等整数条件处判定；有限见证在 K=2^24 处按 `4R > Q` 拒绝、K=2^25 处通过，`Q*J > width·R·(2·factors + 3(t−1) + 1)` 由整除加一的恒等式保证。

扩张环 `O = (Z/3^t)[X]/(F)`：`F` 由完整里程枚举的首个通过候选给出（独立乘法分解器核对每个被拒前缀），系数按 2t 位小端字编码。节点是基数 3 数字向量 1..Q，两两差为非零剩余类、必为可逆元；单位可逆性用牛顿提升（剩余域 Fermat 反演 3^d−2 起步）到达全精度。计数在 `3^t` 以内时由常数系数唯一读回（`natCast_injective_below`）。

Hermite 放置：先算目标多项式在节点处的 J 个精确 Taylor 喷射 `target(X + θ) mod X^J`（整数二项系数，无除法），再用 Newton–Hermite 重建——除节点差外不做除法，`Q*J > degree` 时唯一且与目标逐系数一致；条件数以内的目标精确恢复、达到条件数的目标被检出不等。

## 预算、收据与下一步

`LeanAlgebraBudgets` 限制变量、子句、原始文字与行项数、编码长度、精度、单个支持表项数、活跃系数、逻辑编码位数、单项式乘法和标量操作。`CodeBudgets` 限制矩阵项、图维数、候选、字检查、操作和编码位数；`GreedyBudgets` 限制候选深度、候选数、轮数、标记查询数与历史行/项数。`HighBudgets` 限制单个查询的秩数与条件化查询数；`--high-loop` 只改变查询的执行路径，不改变计数。`ExtensionBudgets` 限制扩张度数、候选、除法检查与环操作；搜索与环构造的拒绝不返回任何见证。拒绝不返回部分码或部分计数，已完成工作保留在收据中。默认 `max_precision=64`、`max_polynomial_terms=16384`、`max_monomial_products=2000000`、`max_scalar_operations=8000000`。超预算返回 `status='refused', count=None`；调用参数无效时抛出 `ValueError`。

原生查询的计量器累计整个查询的操作。流式前沿让全部叶共享乘法和标量操作预算；峰值存储取各叶的最大值，前沿存储单独报告。若第一叶已完成而第二叶拒绝，总计数与 SAT 结果仍为空，已完成工作保留在收据中。独立的 `TruncatedRing` 检查每次稠密运算的工作上界；它不记录任意外部调用序列的累计工作。

活跃的不可变支持表通过弱引用跟踪，另计临时表、输出系数与准备数据。`peak_algebra_bits` 使用每项 n 位支持掩码、固定宽度系数和 framing 的逻辑编码；这些值不是 Python 对象的堆占用、位运算次数或实际磁带空间。

当前 native 计数相当于没有 Z 展开的 `J=1` 查询。虽已实现一般 `t*J` 稠密环，它还没有接到节点求值、修正对和 Hermite 恢复。完整平方自由支持表仍可能指数增长；零布尔扫描不能建立一般多项式空间或次指数时间。下一阶段依次接通流式 Taylor 叶、分组放置，以及整条机器的资源精化。

独立核对入口为 `python -m counting_sat_algebra.verify --step24`（17 项：Newton 系数、共享截断、查询语义、前沿计数和预算拒绝）、`--step25`（15 项：日程整数、距离与 Walsh 首选见证、行再生、候选流、贪心轨迹与流式 hash 叶）、`--step26`（12 项：阈值、截断与掩码、秩和恒等式、路由一致性与预算拒绝）与 `--step27`（12 项：日程整数、守卫边界、不可约搜索、环运算、单位与读回、节点、Hermite 恢复与拒绝）。单元测试还独立核对小精度的 Hermite/Hasse 条件、稠密卷积和小布尔立方语义。`theorem_1_1_complete`、`lean_proof_complete` 与 `recursive_contraction` 均保持 false。
