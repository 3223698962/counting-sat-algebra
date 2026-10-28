# Lean 形状到 Python 义务

包版本 **0.23.0**。这一页只列还能测试的缺口。`theorem_1_1_complete`、`lean_proof_complete` 和 `recursive_contraction` 保持 false。本仓库不导入 Lean，也不把这些检查当作 Theorem 1.1。

Upstream OpenAI Lean（SqrtLog，Sections 3–6 那一层）和论文 Appendix D **不是同一个对象**。下表里的 Lean 文件只提供拒绝条件和字段形状，不翻译成 `Algorithm.lean`。

| 缺口 | Lean 文件提示的形状 | Python 将断言的有限义务 |
| --- | --- | --- |
| D.4.5 加权模式 | `WeightedSupport.lean` 的 `Admissible`：关联多项式上 \(\deg X \le \deg U + \deg Z\)。这不是包模式本身。`admissible_zero`、`admissible_scalar` 和 `retained_degree_bound`（\(j<\mu\)、\(r<D\) 时次数 \(\le \mu+D-2\)）只当作拒绝清单。 | `modular_modes.py`。权重 \(w_Z=b+1\)，源组/目标组，同一权重的模式数受 \(6M/D\) 约束。只保留 \(\sum w_Z\beta_Z < p_t\) 的指数元组。交叉项的截断几何级数与 Step 20 的 Cauchy 条目同余则 `COUNTED`。超预算或次数不合格在分配前 `REFUSED`，`count` 为 null。 |
| D.4.6 包张量 | `SubexpFactorPacketData.lean` 的 `Input` / `Packet`：`labels`、`arity`、`degree`、`cache`、`rows`、`origins`、`sourceRows` 在创建后不变。 | **未实现。** Step 22 只冻结模式描述（符号、创建精度、指数预算、模式表），`tensor_evaluated=false`。不构造 cache、rows、origins，不求值包张量。Step 21 的 `child_topology` 仍只是端口几何。 |
| D.4.7 递归收缩 | 同一 Lean 文件里的 `factor_value_eq` / `FactorReference` 是 SqrtLog 因子语义，不是 (D.16) 的子网络收缩。 | **未实现。** 没有子调用。`recursive_contraction` 保持 false。 |
| D.3 关联编码 | `incidenceArgument` 与 `linearForm` 仍是 SqrtLog 的关联单项式，不是 Appendix D.3 的稀疏系统网络。 | **未实现。** Section 6 的有限方程选择是另一条路线。 |

SqrtLog 的 `Admissible` 通过，只说明声明的次数满足那条不等式。它不证明 Appendix D 的模式恒等式，也不证明 ETH 为假。
