# Roadmap

当前版本 **0.25.0**。主目标是按已加入的 `references/lean/math/lean/OAI/Algorithms/SubexponentialSAT/SqrtLogSpace/Final.lean` 移植完整精确计数程序，并逐项接通语义、时间和全执行前缀空间义务。Step 23 的流式互不重叠前沿、Step 24 的三元计数查询和 Step 25 的日程稀疏码与贪心选行已经落地；完成清单见 [docs/steps/completed.md](docs/steps/completed.md)。

新 Lean 的 `exact_count_and_sat` 提供一个固定原始机器、完整计数字、SAT 位和多项式 visited-tape 空间。其时间尺度的指数仍含 `(log log n)^2 / sqrt(log n)`，详见 `SqrtLog/StatementScale.lean` 与 `SqrtLogSpace/Asymptotics.lean`。Appendix D 的更强界保留为后续路线，不作为这条移植主线的先决条件。旧单文件中的 `sorry` 不再代表新库状态。

## Remaining work — Lean 主线，按依赖顺序

1. [x] **稀疏码和可再生候选行。** 对齐 `SparseCodeSearch.lean`、`SparseCodeRows.lean`、`NativeCodeAssembly.lean`：`lean_code.py` 按日程整数构造三元码并独立核对宽度、距离和出现次数，逐行从槽位再生，不保存候选库；预算失败不保留部分码。Step 25 落地。
2. [x] **真实标记查询驱动的贪心选择。** 对齐 `HashGreedyProgram`、`HashQueryCount.lean`、`NativeLeafRun.lean`：`lean_greedy.py` 用 `count_ternary_query` 的自然数评分累加全部子句的违反计数，逐候选比较并只保留当前最佳，零质量后单独执行只含历史行的接受查询；独立暴力核对零质量时接受集恰好为原 CNF 满足集。Step 25 落地。
3. [ ] **参数前端和高出现变量循环。** 接通 `LeafQueryFrontend`、`HighQueryProgram.lean` 与原声明 n 上的高变量条件化。Step 24 只实现了单个 `fixed` 条件的查询；仍需逐高赋值流式执行、完整累加和次数/出现次数界。
4. [ ] **计数扩张环、分组和替换。** 接通 `RingCore`、`NativeQueryPreparation.lean` 与实际放置：搜索不可约多项式、构造模 `3^t` 的扩张、固定互异插值节点，并建立 `Q*J > degree` 和唯一自然数读回条件。当前只实现基环和有限 `t*J` 截断数组。
5. [ ] **流式 Taylor 查询叶。** 对齐 `BooleanDifferences.lean`、`QueryStreamProgram.lean`、`QueryJetSemantics.lean`、`QueryCountReconstruction.lean`：直接再生 `(E,U)` 修正对，在标量数组内立即求值，完成各节点的全部 J 系数，再做 Hermite 恢复。替换当前可能指数增长的完整平方自由支持表，核对支持次数与截断不变量。
6. [ ] **整条机器的资源与精化。** 将源多项式、固定行/布局和再生游标接到同一实际运行。核对所有成功与失败分支的计数字宽度、算术临时量、输入/前沿/叶存储和每个前缀的空间界。Python 的标量操作统计与逻辑编码位数不是 Lean 原始机器步数或已访问磁带格数，需要独立的精化证明。
7. [ ] **Lean 内核核验与可重复环境。** 使用库声明的 `leanprover/lean4:v4.34.1` 和锁定依赖编译具体入口，保存 `#print axioms` 输出；再核对 Python 实现和对应声明。静态源码审计已经通过，但当前环境尚未运行 Lean/Lake。

## Remaining work — Appendix D 后续路线

- [ ] 将 D.4.5 从单个有限包模式扩展到一般模式展开。
- [ ] 实现 D.4.6 包张量系数及不可变创建时描述。
- [ ] 实现 D.4.7 子网络递归收缩并核对边、度数和精度不变量。
- [ ] 将 D.3 关联单项式编码与一般密度查询接到该递归。
- [ ] 单独建立 Appendix D 的更强位时间与空间界。

每个里程碑必须给出可运行入口、源声明映射、独立有限对照和无部分计数的预算拒绝。`theorem_1_1_complete`、`lean_proof_complete`、`recursive_contraction` 当前均为 false。
