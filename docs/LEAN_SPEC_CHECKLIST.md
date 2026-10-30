# Lean 声明到 Python 义务

包版本 **0.28.0**。新加入的 Lean 库以 `SqrtLogSpace/Final.lean` 为精确计数终点。下表按源码依赖列出 Python 已实现的有限义务和下一阶段连接；路线顺序见 [ROADMAP.md](../ROADMAP.md)，代数定义见 [LEAN_NATIVE.md](LEAN_NATIVE.md)。Lean 路径均相对于 `references/lean/math/lean/OAI/Algorithms/SubexponentialSAT/`。

| 主线义务 | Lean 来源 | 当前 Python 状态与验收条件 |
| --- | --- | --- |
| 原始编码与计数字 | `Encoding.lean`、`SqrtLogSpace/DisjointReuse/CountModel.lean` | Step 23 已实现有限往返、原声明 n 和小端计数字；计数字长度不超过 n+1。继续核对全前端解析与机器端口。 |
| 二进参数与互不重叠前沿 | `SqrtLogSpace/DisjointReuse/` 下的 `SparseBudget.lean`、`SparsePriority.lean`、`SparseCountProgram.lean` | Step 23 已实现首个通过精确整数条件的 M、优先核心和深度优先计数；叶满足集互不重叠、并集等于输入，已完成叶不保留。结构预算先预检，叶拒绝不返回部分总和。 |
| 实际三元指示系数 | `SqrtLog/` 下的 `IdempotentCore.lean`、`IdempotentPolynomial.lean`、`IdempotentBase.lean` | Step 24 已计算 Newton 迭代、单首 digit 除法和 2t 位系数编码；用独立小精度 Hermite/Hasse 条件核对，禁止外部指示表替代构造。 |
| 共享标记与布尔求和 | `SqrtLogSpace/TruncatedRing.lean`、`SqrtLog/SquarefreeBasic.lean` | Step 24 已实现独立 U/Z 截断、完整乘积后读出、支持并集乘法与单项式权重求和。当前查询为 J=1，完整支持表仍可能指数增长。 |
| 接受/标记查询 | `SqrtLogSpace/HashQueryCount.lean` | Step 24 已保留原子句编号、重复和 ambient n，`fixed` 只条件化原坐标；核对接受数、标记违反数及 `3^t>2^n` 唯一读回。单位行计数才直接给出原 CNF #SAT。 |
| 可再生稀疏码 | `SqrtLogSpace/` 下的 `SparseCodeSearch.lean`、`SparseCodeRows.lean`、`NativeCodeAssembly.lean` | Step 25 已按日程整数生成并独立核对宽度、距离和出现次数；逐行从槽位再生，预算失败不保留部分码。Lean 日程深度超出 Python 预算时显式拒绝。 |
| 查询评分与贪心选择 | `SqrtLogSpace/` 下的 `HashGreedyProgram`、`NativeLeafRun.lean`、`HashQueryCount.lean` | Step 25 已接通：对候选行逐子句累加 `count_ternary_query` 标记违反计数，严格最小只保留当前最佳；零违反质量时单独执行只含历史行的接受查询，独立暴力核对接受集恰为原 CNF 模型集。 |
| 高变量与查询准备 | `SqrtLogSpace/HighQueryProgram.lean`、`SqrtLogSpace/NativeQueryPreparation.lean`、`SqrtLog/LeafQueryFrontend` | Step 26 已接通有限版本：早期阈值、截断、高变量掩码与两条上界；同一原 n 上逐秩流式遍历高赋值并完整累加，`--high-loop` 将 hash 叶的全部查询路由到该循环。完整参数前端与出现次数守卫仍随整机精化待办。 |
| 扩张环、放置与 Taylor 叶 | `SqrtLog/RingCore.lean` 与 `SqrtLogSpace/` 下的 `BooleanDifferences.lean`、`QueryStreamProgram.lean`、`QueryJetSemantics.lean`、`QueryCountReconstruction.lean` | Step 27 已接通环与节点部分（日程整数、不可约单首搜索、模 `3^t` 扩张、固定互异节点、`Q*J > degree` 与唯一自然数读回，并有精确 Taylor/Newton–Hermite 有限见证）；Step 28 已接通修正对分解（`fold` 顺序流）、逐节点 J 系数块与 Newton–Hermite 读回。把生产叶切换到该流式路径（含立即 `U` 特化与完整修正对再生）仍待接通。 |
| 原始机器与全前缀资源 | `SqrtLogSpace/Final.lean`：`compiled_exact_count_and_sat`、`exact_count_and_sat`、`exact_count_eth` | **未完成 Python 精化。** 标量操作与逻辑编码位数不等于机器步数或 visited-tape 空间。需核对整条成功/失败路径并在锁定环境内编译入口、保存实际 `#print axioms` 输出。 |

静态审计已经覆盖 3969 个本地导入模块，没有代码级占位符或本地公理声明；`kernel_verified=false`。`lean_proof_complete=false` 记录待执行的本地核验与 Python 精化，不表示新库终点仍是旧脚手架中的 `sorry`。`theorem_1_1_complete` 与 `recursive_contraction` 也保持 false。

## Appendix D 后续义务

SqrtLog/SqrtLogSpace 与论文 Appendix D 的对象和时间指数分别记录。下表保留原有有限模式的形状提示；它们不构成新主线移植的先决条件。

| 缺口 | Lean 文件提示的形状 | Python 将断言的有限义务 |
| --- | --- | --- |
| D.4.5 加权模式 | `WeightedSupport.lean` 的 `Admissible`：关联多项式上 \(\deg X \le \deg U + \deg Z\)。这不是包模式本身。`admissible_zero`、`admissible_scalar` 和 `retained_degree_bound`（\(j<\mu\)、\(r<D\) 时次数 \(\le \mu+D-2\)）只当作拒绝清单。 | `modular_modes.py`。权重 \(w_Z=b+1\)，源组/目标组，同一权重的模式数受 \(6M/D\) 约束。只保留 \(\sum w_Z\beta_Z < p_t\) 的指数元组。交叉项的截断几何级数与 Step 20 的 Cauchy 条目同余则 `COUNTED`。超预算或次数不合格在分配前 `REFUSED`，`count` 为 null。 |
| D.4.6 包张量 | `SubexpFactorPacketData.lean` 的 `Input` / `Packet`：`labels`、`arity`、`degree`、`cache`、`rows`、`origins`、`sourceRows` 在创建后不变。 | **未实现。** Step 22 只冻结模式描述（符号、创建精度、指数预算、模式表），`tensor_evaluated=false`。不构造 cache、rows、origins，不求值包张量。Step 21 的 `child_topology` 仍只是端口几何。 |
| D.4.7 递归收缩 | 同一 Lean 文件里的 `factor_value_eq` / `FactorReference` 是 SqrtLog 因子语义，不是 (D.16) 的子网络收缩。 | **未实现。** 没有子调用。`recursive_contraction` 保持 false。 |
| D.3 关联编码 | `incidenceArgument` 与 `linearForm` 仍是 SqrtLog 的关联单项式，不是 Appendix D.3 的稀疏系统网络。 | **未实现。** Section 6 的有限方程选择是另一条路线。 |

SqrtLog 的 `Admissible` 通过，只说明声明的次数满足那条不等式。它不证明 Appendix D 的模式恒等式，也不证明 ETH 为假。
