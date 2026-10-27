# Step 20：有限模网络与局部恢复

版本 0.21.0 实现 Appendix D.1–D.2 的布尔张量网络、D.4.1 的有限系数环，以及 D.4.3–D.4.4 的 Cauchy 变形和分层精度 Hermite 恢复。网络叶节点枚举边赋值；变形求值枚举行端口和列端口的独立赋值。它们都有显式预算。

本文保留 Step 20 的实现范围。版本 0.22.0 的 [Step 21](PACKET_PACKING.md) 已补充 D.4.2 的整数日程与碰撞分包，以及 D.4.6 的子网络端口拓扑，并可将自动布局交给这里的恢复接口。D.4.5 的加权模式、包张量系数和递归收缩，以及 D.3 的稀疏系统网络编码仍未实现。Theorem D.1 的一般收缩算法、Theorem 1.1 和相应的时间空间界尚未完成；`theorem_1_1_complete`、`lean_proof_complete` 和 `recursive_contraction` 保持 false。

## 接口与表示

| 文件 | 内容 |
| --- | --- |
| `modular_ring.py` | 模系数环、不可约多项式搜索、单位求逆、逐坐标提升、有限域插值 |
| `modular_network.py` | 带有序端口的二部多重图、隐式张量、直接叶求和、CNF 编码与计数 |
| `modular_recovery.py` | 给定的包布局、Cauchy 局部系数、D.8 求值、D.5 精度计划与 D.9 恢复 |

`ModularRing(p, ell, polynomial)` 表示

\[
A_p=((\mathbb Z/2^p\mathbb Z)[U]/f)[V]/(V^\ell).
\]

`polynomial` 的二进制位给出首一多项式 \(f\) 的 0/1 系数；模 2 后必须不可约。默认 `0b11` 的扩张次数为 1，给出基环 \(B_{p,\ell}\)。`extension_ring_for_size(s, p, ell)` 按论文取 \(e=4\lceil\log_2s\rceil\)，在预算内确定性搜索 \(f\)。例如 \(s=2\) 选出 \(U^4+U+1\)。

每个元素是长度 \(e\ell\) 的整数元组。下标 `j*e+i` 是 \(U^iV^j\) 的系数，取值为 `0..2**p-1`。`from_base` 嵌入 \(V\) 系数；`base_coefficients` 在所有非恒定 \(U\) 坐标为零时提取它们。扩张中的约化使用整数减法；例如 \(f=U^2+U+1\) 时 \(U^2=-U-1\)，不能在模 \(2^p\) 的运算中使用异或代替减法。

单位求逆先在剩余域求逆，再做 Newton 提升并逆转截断的 \(V\) 级数。插值点是剩余域元素的二进制坐标提升，其两两之差是单位。普通整数点 `0, 2` 的差不能在模 \(2^p\) 下求逆。

## CNF 网络

```python
from counting_sat_algebra import count_cnf_modular
from counting_sat_algebra.dimacs import read_dimacs

result = count_cnf_modular(read_dimacs("examples/count_demo.cnf"))
assert result.counted and result.count == 196
```

每个子句有一个满足张量，每个不同变量占一个端口。重复文字合并，重言式删除，空子句直接返回零。出现一次的变量接一元全 1 张量，出现两次接二元相等张量；出现 \(d\ge3\) 次的变量使用 \(d-2\) 个三价相等顶点组成树。每条原始边再插入一个二元相等顶点，构成二部网络。端口顺序被保存，平行边使用不同编号。

网络在模 \(2^{n+1}\) 下收缩，再乘以声明的自由变量因子。8 变量示例有 6 个子句端口，细分后有 12 条边和 14 个顶点，最大度为 3。它扫描 4096 个边赋值、求值 57344 次张量条目，得到 196。3 个声明变量上的四条重复正单元子句使用两个三价相等顶点；10 条边、11 个顶点的网络给出计数 4。

`count_cnf_modular` 是独立的 Python 接口。DIMACS `count` 命令仍通过原有的 `common` 或 `tensor` 后端计数。

## 创建时精度

`TensorVertex` 保存创建时的环和不接收当前精度的条目函数。每次求值都先在创建精度下归约，随后逐坐标取规范整数代表，并提升到当前调用的精度。条目函数应当确定性求值，其捕获的描述应不可变。

例如创建精度为 2 时，条目 `-1` 与 `6` 的代表分别为 3 和 2。提升到精度 8 后，每个边赋值贡献 \(3\cdot2=6\)，两项之和是 12。不能先在创建精度下乘为 2，再把乘积提升。规范提升不保持加法和乘法。

## Cauchy 变形与恢复

`PacketLayout` 显式指定每个边在一个 \(L=2^a\) 槽位包中的位置；`None` 是空槽。边在行列两端使用同一槽号。单独给定布局没有符号分组，也不声称满足 Lemma D.3 的负载界。Step 21 的 `pack_network` 可以生成带有标签和实际负载证书的布局；它分别检查 D.6 数值前提和有限分包结论。

在 \(S=\{-L,\ldots,L-1\}\)、\(\alpha=2^{a+1}\) 上，代码计算

\[
W_{xy}(T)=\frac{\prod_{q\in S}(T-q)}{C'(0)(T-(y-x))}.
\]

先消去分子中的指定因子，再从 \(t+\alpha z-q\) 提取已知的 2 的幂，只对分母的奇数部分求逆。这样得到的局部多项式系数属于 \(A_p\)。在 \(a=1,t=1,x=0,y=1\) 时，系数是 \((-3,-22,-48,-32)\)。\(W(0)=I\)。

`evaluate_deformation` 按 D.8 扫描 \(2^{2M}\) 个端口赋值对。`recover_network_modular` 对每个活动节点使用 \(M(2L-1)+1\) 个剩余域插值点恢复整个局部多项式，再保留 D.5 所需的前 \(q_t\) 项。精度为

\[
h=\lceil3M/u\rceil,\quad q_t=h\max(0,u-v_2(t)),\quad
p_t=p+8hL/2^{v_2(t)}+8M(a+1).
\]

D.9 的商级数在模环内计算；最后以 `Fraction` 合并有理数权重，再约去 2 的幂并归约到原精度。不会试图在模环内直接求逆偶数分母。

```python
from counting_sat_algebra import (
    BooleanTensorNetwork, ModularRing, PacketLayout, TensorVertex,
    embed_network, extension_ring_for_size, recover_network_modular,
)

base = ModularRing(5, 2)
network = BooleanTensorNetwork(base, 2, (
    TensorVertex("row", 0, (0, 1), base, lambda b: (1 + b[0], 2 - b[1])),
    TensorVertex("column", 1, (0, 1), base, lambda b: (1 + b[1], -b[0])),
))
field = extension_ring_for_size(2, 5, 2)
network = embed_network(network, field.polynomial)
result = recover_network_modular(network, PacketLayout(2, ((0, None, 1, None),)), 2)
assert result.recovered
assert field.base_coefficients(result.value) == (9, 4)
assert result.port_pairs_scanned == 1440
```

这个例子在节点 \((-3,-2,-1,1,2,3)\) 的重数为 \((6,3,6,6,3,6)\)，精度为 \((149,101,149,149,101,149)\)。恢复结果为 \(9+4V\)（模 32，\(V^2=0\)），与独立直接收缩一致。\(a=u=1\) 时，两个节点均使用 133 位精度，扫描 224 对端口赋值。另一个非线性行张量产生确实达到 6 次的变形，恢复 \(8+5V\)，与独立有理数多项式逐系数核对。

## 预算与复现

`ModularBudgets` 分别限制边、顶点、扫描量、张量条目求值、精度、系数坐标、扩张次数、包槽位、插值点、喷射存储和有理数恢复位数。CNF 计数在网络分配前检查扫描预算；恢复流程先计划全部节点，随后才分配核和喷射、调用条目函数。拒绝返回 `count=None` 或 `value=None`，相应扫描计数为零。

喷射坐标预算不包含 Python 对象开销、规划元数据和算术临时对象；有理数位数预算使用所有恢复权重分母之积的保守界。这些局部资源上界不是整个递归算法的位复杂度或空间证明。

```bash
python -m pytest -q tests/test_step20_modular_recovery.py
python -m counting_sat_algebra.verify --step20 --output artifacts/step20-verification.json
python -m counting_sat_algebra.verify --reference --output benches/reference_alignment.json
```

Step 20 有 14 项参考检查，报告区分实际 CNF 计数 `COUNTED`、环与恢复恒等式 `CHECKED`，以及枚举或存储预算拒绝 `REFUSED`。完整回归和当前剩余工作见 [STATUS.md](STATUS.md) 与 [ROADMAP.md](../ROADMAP.md)。
