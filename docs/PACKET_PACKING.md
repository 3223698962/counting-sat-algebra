# Step 21：碰撞标签、分包与子网络拓扑

版本 0.22.0 在 `modular_packing.py` 中实现 Appendix D.4.2 的整数日程、精确条件期望标签和统一填充分包，并构造 D.4.6 的子网络端口拓扑。自动生成的 `PacketLayout` 可以直接交给 [Step 20 的局部恢复](APPENDIX_D.md)。分包和拓扑构造均不求值张量条目。

这一阶段尚未实现 D.4.5 的加权模式、D.4.6 的包张量系数及其不可变创建时描述、D.4.7 的递归收缩，以及 D.3 的关联单项式编码。`recursive_contraction`、`theorem_1_1_complete` 和 `lean_proof_complete` 均保持 false。

## 接口

以下接口均从 `counting_sat_algebra` 导出。

| 接口 | 返回内容 |
| --- | --- |
| `appendix_d_integers(lambda_s, depth=0)` | D.4 的整数 `u`、`a`、`previous_a`；不构造输入规模 `s` |
| `appendix_d_proposal(s, edges, depth=0)` | 叶阈值、度数上界、活动调用守卫和预算内的 D.6 检查 |
| `packing_hypotheses(edges, r, a, delta)` | 三个数值前提及逐项失败原因 |
| `choose_collision_labels(network, r, a, trace=False)` | 确定性标签、初始期望、最终势、标签对负载及操作统计 |
| `pack_network(network, r, a, delta=None, require_hypotheses=True, trace=False)` | 分包结果、预检计划和实际证书检查 |
| `check_packing(packed)` | 重新检查负载、填充、槽位方向、组端口和子图规模 |
| `child_topology(packed, coordinate, group)` | 指定坐标与符号组的原顶点、包顶点及有序端口 |

上述构造接口除 `appendix_d_integers` 和 `check_packing` 外都接受 `budgets=PackingBudgets(...)`。`pack_network` 默认要求 D.6 全部成立；`delta=None` 时使用 `max(network.max_degree, 2*L)`。传入低于原网络最大度的 `delta` 是无效参数。

`require_hypotheses=False` 显式请求有限实验，仍须通过实际包数、组端口数和度数检查。`PackingResult.packed` 才表示返回了可用分包；`refused` 时 `packing=None`。独立标签选择不检查这些分包结论，超预算会抛出 `ModularBudgetRefusal`；子拓扑接口也以该异常报告预算拒绝。

## 整数日程与 D.6

代码按整数算术计算

\[
\lambda=\lceil\log_2s\rceil,\quad
H=\left\lceil\frac{s}{\lfloor\sqrt\lambda\rfloor}\right\rceil,\quad
\Delta=\lceil s^{9/10}\rceil,\quad r=64,\quad\epsilon=1/131072,
\]

\[
u_d=\left\lfloor\frac{\lambda}{131072\,2^d}\right\rfloor,\quad
a_d=2\sum_{j=0}^d u_j,\quad L=2^{a_d},\quad D_{\rm arr}=r^{2a_d}.
\]

`Delta` 通过整数十次根求出，不使用浮点幂。求和利用 `sum floor(t/2**j) = 2*t - popcount(t)` 及其尾和，不循环遍历可能很大的深度。`edges <= H` 只返回叶阈值说明；提案接口不执行直接求和。活动调用若不满足 `1 <= u <= M`，会在形成 `L` 和 `D_arr` 前停止。

D.6 的三个条件是

\[
M\ge 2D_{\rm arr}^{\,2}\Delta,\qquad
M\ge LD_{\rm arr},\qquad 2L\le\Delta.
\]

第一个条件包含 **`D_arr` 的平方**。例如 `r=64, a=1, delta=4` 时，`D_arr=4096`，该条件要求 `M >= 134217728`。`M=8192` 不通过此条件，即使该图的实际分包和子图规模都符合结论。

```python
from counting_sat_algebra import appendix_d_integers, appendix_d_proposal

proposal = appendix_d_proposal(16, 20)
assert (proposal.H, proposal.delta, proposal.integers.lambda_s) == (8, 13, 4)
assert proposal.integers.u == proposal.integers.a == 0
assert not proposal.leaf and proposal.hypotheses is None

witness = appendix_d_integers(131072)
assert (witness.u, witness.a, witness.previous_a) == (1, 2, 0)
deeper = appendix_d_integers(262144, depth=1)
assert (deeper.u, deeper.a, deeper.previous_a) == (1, 6, 4)
```

这些对数整数见证没有构造相应的巨大 `s`，也没有执行它的递归。一般有限元组通过 D.6，并不使它成为 D.4 日程的成员；D.4 固定 `r=64`，且正的 `a` 必为偶数。

## 精确条件期望

设 `K=r**a`，`X_gamma` 为标签对的边负载，势为

\[
\Phi=\sum_\gamma(X_\gamma-M/K^2)^2.
\]

平行边共享两个端点，其他相邻边共享一个端点，不能视为独立。令 `m_uv` 是平行边重数、`d_v` 是端口度数、`P=sum m_uv**2`，则初始期望满足

\[
K^2\mathbb E\Phi=(K^2-1)P+(K-1)\left(\sum_vd_v^2-2P\right).
\]

实现维护已完成标签对的负载和已赋标签邻居的重数直方图。对于当前行顶点 `v`，记其度数为 `d`、未赋标签的邻居端口数为 `z`、已赋邻居直方图为 `d_b`、完成负载为 `C_ab`、仅行端已赋标签的边数为 `R_a`。除去与候选无关的项，候选标签 `a` 的分数为

\[
c_a=K\sum_b C_{ab}d_b+z\sum_b C_{ab}+dR_a
 +(K-1)\sum_{w\text{ 未赋标签}}m_{vw}d_w(a).
\]

选择该候选后，`K**2 * E[Phi]` 的变化为 `2*(K*c_a - sum(c))`。列顶点使用转置公式。按源网络顶点顺序选择，平局取最小标签；内部更新使用整数，公开期望值使用 `Fraction`。这避免了为每个候选重新枚举所有边对。

```python
from fractions import Fraction
from counting_sat_algebra import (
    BooleanTensorNetwork, ModularRing, TensorVertex, choose_collision_labels,
)

def no_entry_evaluation(bits):
    raise AssertionError("labeling must not evaluate tensors")

ends = ((0, 2), (0, 2), (0, 3), (1, 2), (1, 3), (1, 3))
ring = ModularRing(5)
network = BooleanTensorNetwork(ring, len(ends), tuple(
    TensorVertex(str(v), int(v >= 2),
                 tuple(e for e, pair in enumerate(ends) if v in pair),
                 ring, no_entry_evaluation)
    for v in range(4)
))
labels = choose_collision_labels(network, 2, 1, trace=True)
assert labels.initial_potential == Fraction(23, 2)
assert labels.final_potential == 1
assert labels.labels == (0, 1, 0, 1)
assert labels.pair_loads == ((2, 1), (1, 2))
```

`trace=True` 保存每一步所有候选的精确期望；默认关闭，并有独立存储预算。测试通过枚举尚未赋标签的顶点逐项核对每个候选，覆盖交错排列的二部顶点、共享端点、平行边和固定种子的随机图。

## 统一填充与子网络

每个标签对填充到统一的 `H_pad`，即最大实际负载向上取整到 `L` 的倍数；空槽用 `None` 表示，不产生端口。标签的第 `i` 个符号是 `(label // r**i) % r`。固定每一坐标上的无序符号对后，槽号的第 `i` 位决定不同符号的方向；当两个符号相等时，该位选择同一标签对中的不同副本。后者防止边被重复使用。

生成的包数恰为 `H_pad * D_arr / L`。证书检查每条原边恰出现一次、两端槽号相同、方向正确、每侧每组端口不超过 `3M/r`、包数不超过 `3M/L`，以及实际最大度不超过 `Delta`。它同时重新计算势、标签对负载和 D.6 标志，拒绝这些字段之间的不一致。

`child_topology` 将原网络中符号等于 `group` 的顶点放在一侧，所有包放在另一侧。每个原端口连接到其所在包；`port_origins` 记录 `(source_vertex, source_edge)`，`original_ports` 保留原端口顺序。原孤立顶点与零度包顶点也会保留。子拓扑边数是该组两侧原端口数之和，受 `6M/r` 约束；包顶点度数最多为 `2L`。`child_shrink` 另行记录是否实际达到 `M/8`，其为 false 并不会使较小 `r` 的有限分包失效。

```python
from counting_sat_algebra import (
    BooleanTensorNetwork, ModularRing, TensorVertex, child_topology, pack_network,
)

def no_entry_evaluation(bits):
    raise AssertionError("packing and topology must not evaluate tensors")

ring = ModularRing(1)
network = BooleanTensorNetwork(ring, 128, tuple(
    TensorVertex(f"{side}-{edge}", side, (edge,), ring, no_entry_evaluation)
    for side in (0, 1) for edge in range(128)
))
result = pack_network(network, 2, 1, delta=4)
assert result.packed and result.packing.hypotheses.admitted
assert result.packing.H_pad == 32
assert len(result.packing.layout.packets) == 64
child = child_topology(result.packing, 0, 0)
assert (child.edges, child.max_degree) == (128, 4)
assert result.check.ok and not result.check.child_shrink
assert result.tensor_evaluations == 0
```

| 匹配网络 | 实际结果 | 前提与范围 |
| --- | --- | --- |
| `M=128, r=2, a=1, delta=4` | 势 `96→0`，64 个包，最大子拓扑 128 条边、度数 4 | D.6 成立；不满足 `M/8` 缩减；不是 D.4 日程成员 |
| `M=32, r=2, a=2` | `H_pad=4`，16 个包；覆盖 0、1、2 个相等坐标 | 显式有限模式；D.6 不成立 |
| `M=8192, r=64, a=1, delta=4` | 势 `8190→0`，4096 个标签对各负载 2，4096 个包；每侧每组 128 个端口；最大子拓扑 256 条边、度数 4 | 默认在标签搜索前拒绝；显式有限模式通过实际结论，且 `256 <= M/8=1024`；D.6 仍不成立，也不是 D.4 日程成员 |

8192 边实例执行 1048576 次候选标签试验、1572864 次分数更新，张量求值为零。子拓扑本身没有包张量系数、模式权重或收缩值。

## 自动布局接入局部恢复

```python
from counting_sat_algebra import (
    embed_network, encode_cnf_network, extension_ring_for_size,
    pack_network, recover_network_modular,
)
from counting_sat_algebra.cnf import CNF

cnf = CNF(2, [[(0, False)], [(0, False)]])
encoded = encode_cnf_network(cnf)
assert encoded.network is not None and encoded.free_variables == 1
field = extension_ring_for_size(2, encoded.network.ring.precision)
network = embed_network(encoded.network, field.polynomial)
packed = pack_network(network, 2, 1, require_hypotheses=False)
assert packed.packed
recovered = recover_network_modular(network, packed.packing.layout, 1)
assert recovered.recovered
assert field.base_coefficients(recovered.value) == (1,)
count = field.base_coefficients(recovered.value)[0] << encoded.free_variables
assert count == 2
assert recovered.port_pairs_scanned == 6656
assert recovered.tensor_evaluations == 33280
```

该 CNF 的网络有 4 条边，局部恢复剩余为 1，声明的自由变量使最终计数为 2。此处执行的是 Step 20 的端口赋值枚举与 Hermite 恢复；自动分包尚未接入递归求值。DIMACS `count` 命令继续提供原有 `common`、`tensor` 后端。

## 预算和复现

`PackingBudgets` 默认限制 65536 条边、131072 个顶点、4096 个标签、65536 个标签对类别、每包 256 个槽位、5000000 次标签试验、20000000 次分数更新、2000000 个填充槽位和 1000000 个包。可选跟踪最多 4096 个候选期望；子图最多 131072 条边；参数中间整数最多 4096 位。

标签与包的资源检查全部在搜索前完成。设初始期望为 `E0`，使用

\[
b=\min\left(M,\left\lceil M/K^2\right\rceil+\left\lceil\sqrt{E_0}\right\rceil\right),\qquad
B_{\rm slots}=L\lceil b/L\rceil K^2
\]

作为填充槽位上界。候选试验上界为 `K * |vertices|`，分数更新上界为 `K * |vertices| + K*M + sum(d_v**2)`。这些是保守界：128 边匹配的填充预检上界为 168 个槽位，实际只需 128 个；预算为 167 时仍在标签搜索前拒绝。

预算或默认 D.6 拒绝返回零标签试验，且不分配标签表或包。显式有限模式若在构造后未通过实际结论，也返回 `packing=None`，但会如实记录已完成的标签试验、分数更新及分配标志。子拓扑在已有分包的基础上单独检查边数、顶点数和包数预算，再构造端口表。

计数器不包含所有 Python 对象开销、元数据或算术临时对象；这些局部预算不是 Theorem D.1 的位复杂度或空间证明。

```bash
python -m pytest -q tests/test_step21_modular_packing.py
python -m counting_sat_algebra.verify --step21 --output artifacts/step21-verification.json
python -m counting_sat_algebra.verify --reference --output benches/reference_alignment.json
```

Step 21 有 33 项 pytest 和 14 项参考检查。版本 0.22.0 的完整回归为 322 项 pytest、84 项参考检查，覆盖 Steps 9–21；`verify --quick` 与两项 CI 报告校验均通过。报告只记录结果和操作统计，不记录运行时间。剩余实现见 [ROADMAP.md](../ROADMAP.md)。
