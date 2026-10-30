# 与 references 的实现对齐

当前版本 **0.29.0**。主目标是移植新加入的完整 Lean `SqrtLogSpace` 精确计数程序，并核验 Python 的语义、时间和空间连接。Step 23 已实现原声明 n 上的流式互不重叠前沿，Step 24 已实现三元指示系数、接受/标记查询及原生计数叶，Step 25 已实现日程稀疏码、逐行再生与违反质量贪心选行，Step 26 已实现早期阈值、高变量前端与流式秩循环，Step 27 已实现计数扩张环、固定节点与 Hermite 放置条件，Step 28 已实现流式修正对与 Taylor 读回，Step 29 已实现无支持表的逐点修正对叶。原论文 Sections 3–6 和 Appendix D 的有限结果保留；更强的 Appendix D 路线列为后续工作。

本包尚未实现一般输入的完整渐近算法，`theorem_1_1_complete`、`lean_proof_complete` 与 `recursive_contraction` 均为 false。新 Lean 库具有终点证明源码；这些标志记录 Python 实现、本地内核核验和精化的未完成状态，不把新库误写成旧的 `sorry` 脚手架。

## 需要达到的结论

参考文献的 Theorem 1.1 要求：对固定宽度 k 的任意 k-CNF，准确返回满足赋值数，并在位运算模型下满足

\[
(n+m)^{O_k(1)}2^{O_k(n\log\log n/\sqrt{\log n})}
\]

时间与 \((n+m)^{O_k(1)}\) 空间。当子句密度固定时，指数改善为 \(O_{k,D}(n/\sqrt{\log n})\)。这些更强的界来自 Appendix D 的张量收缩。

新 Lean 的 `SqrtLogSpace/Final.lean` 给出 `compiled_exact_count_and_sat`、`exact_count_and_sat` 和 `exact_count_eth`：一个固定程序返回完整自然数计数字及 SAT 位，计数字长度不超过 n+1，并给出全部执行前缀的多项式 visited-tape 空间界。`SubexpOrderedRawMainTheorem.lean/mainTheorem_explicit` 保留原始外部 SAT 编码签名。

时间尺度采用较弱的指数 \(O_k(n(\log\log n)^2/\sqrt{\log n})\)。`SqrtLog/StatementScale.lean` 的具体定义为

\[
runtimeScale(n)=n\frac{(\log\log(n+16))^2}{\sqrt{\log(n+16)}}.
\]

`SqrtLogSpace/Asymptotics.lean` 则以底数 2 的对数定义 rho，并把它连接到该尺度。它与原 Sections 3–6 的指数同阶，不能写成 Appendix D 的更强指数。主路线按这套计数程序的依赖实施，无需先补完 Appendix D 递归。

静态审计覆盖 3969 个本地导入模块，没有缺失的本地导入、代码级 `sorry`/`admit` 或 `axiom` 声明；所列入口具有证明体或证明项。`kernel_verified=false`，本机尚未执行 Lean 编译。旧 `references/lean/SubexponentialSAT.lean` 的 `sorry` 只属于历史脚手架。审计边界和可重复编译入口见 [Lean 参考说明](../references/lean/README.md)。

因此，准确计数、小实例检查、一般输入的复杂度分析和完整 Lean 证明分别需要证据。Python 的测试通过不代表时间界或 Lean 定理已经完成。

## 当前完成的连接

| 参考位置 | 代码 | 实际完成内容 |
| --- | --- | --- |
| Lemma 4.5 | `conditional_balance.py` | 精确最小化 \(\mathbb E[\Phi/(M(b+w))+\Xi/(Qn)]\)；包含共享变量和共享因子的碰撞依赖；逐步期望、初值和终值用 Fraction 核对 |
| Definition 3.1 的核条件 | `linear_algebra.py`、`replacement_matrix.py` | 每个组的测量矩阵消元一次；检查目标向量属于行空间，从而覆盖整个有理数核；布尔枚举作为显式可选对照 |
| Lemma 4.7 | `reference_parameters.py` | 检查平衡不等式、R<Q、区间和精确二项式尾界；不分配矩阵 |
| Lemma 4.8 的整数日程 | `reference_parameters.py` | 按证明顺序计算 λ、ν、β、r、Q、R；仅在 r≥2 且 Q>R 之后才形成 τ、D 和中心区间。机器可放置的 n 上该闭式在形成 D 之前失败 |
| 预算内的有限参数成员 | `reference_parameters.py`、`sparse_placement.py` | 只保留 `check_reference_parameters` 接受的元组；τ 取满足精确尾界的最小允许值。通过后再调用既有稀疏放置和计数后端。这不是渐近时间或空间界的证明 |
| Lemma 3.2 | `common_count.py` | 流式枚举整数测量元组；两次扫描每个局部布尔组；代表元不必满足固定比特条件；不存可达测量表 |
| Lemma 3.3、式 (3.13) | `common_count.py` | 枚举 \(U\subseteq E\)、\(|E|\le q\)，恢复截断 Taylor 系数；在 \(\mathbb Q[z]/(z^{q+1})\) 内求值，再做 Hermite 恢复 |
| 式 (4.7)、(4.10) | `structured_auxiliary.py` | 通过交换轨道的 Krawtchouk 系数查询 M0 和 P 的条目；在 R=32、Q=192 上也可查询；不分配完整稠密数组 |
| Lemmas 4.5–4.8 的有限放置 | `sparse_placement.py`、`linear_algebra.py` | 只存占用坐标；有界稀疏消元检查核条件、恢复、分母及系数范数；自动日程入口检查全部 Lemma 4.7 数值前提 |
| 独立公式乘积与实际输入 | `component_count.py`、`dimacs.py`、`cli.py` | 规范化 CNF、分离独立分量、分别放置并代数计数；声明自由变量贡献 2 的幂；支持预算拒绝与精确大整数输出 |
| 不相交支持集乘积与局部 Taylor 恢复 | `support_count.py` | 对同一全局放置的实际替换因子按支持集分解；局部布尔求和与 Taylor 恢复得到全局 H(0)；有限收缩，不是 Appendix D 的一般算法 |
| 连通支持上的有限 Taylor 消元 | `tensor_count.py`；`cli.py` 的 `--backend tensor` | 在 \(\mathbb Q[z]\) 上对实际替换因子做预算内变量消去；表是 T=1 的完整 Taylor 喷射。只处理中间表较小的连通支持。不是 Appendix D 的一般模收缩，也不建立其渐近界。拒绝发生在建表之前，且不返回部分计数 |
| Section 5.1–5.2 的有限系数代数 | `sparse_systems.py`；`verify --step17` | 稀疏 \(\mathbb F_2\) 单项式、整数提升、Lemma 5.1 的精确除法，以及按 \(h\) 次截断的乘积。公共零点在预算内的立方上穷举。不是 Section 5.3 的替换求和，也不是 Theorem 2.2 的运行时间。拒绝发生在项表或扫描之前，且不返回部分计数 |
| Section 5.3–5.5 的有限替换求和 | `algebra_count.py`；`common_count.py` 的 `separated_group_sum`；`verify --step18` | 方程因子 \(B_i\) 在 \(A_\ell\) 中按 (5.4) 构造，并加上值为 1 的中性因子。Section 3 的测量扫描在 \(A_\ell[z]/(z^J)\) 中求和，再对每个 \(V\) 系数做 Hermite 插值。有限设计满足 Definition 3.1，但不是 Lemma 4.8 的行。高出现变量按 Section 5.5 分支。拒绝发生在枚举或分配喷射之前，且不返回部分计数。不是 Section 6，也不是 Theorem 2.2 的运行时间 |
| Section 6 的有限方程选择与互不重叠稀疏化 | `equation_selection.py`；`verify --step19` | 子句违反多项式、Lemma 6.1 的整数和预算内的 \(A(0)\)、贪心选行、同一声明 \(n\) 上的互不重叠叶，以及 Proposition 6.4 的整数。拒绝发生在分配或枚举之前，且不返回部分计数。不是 \(O_k(\rho(n))\)，也不是 Appendix D。Step 5 的蕴涵环复制会增加 \(n\)，不是 Lemma 6.3 |
| Appendix D.1–D.2 的有限网络 | `modular_network.py` | 有序端口、平行边、三价相等树和逐边细分；直接模叶求和保留自由变量因子。示例计数 196。扫描预算在图分配前核对 |
| Appendix D.4.1 的有限系数环 | `modular_ring.py` | 不可约扩张、Newton 单位求逆、剩余域插值、创建时精度的逐坐标提升。环参数和搜索受预算限制 |
| Appendix D.4.3–D.4.4 的局部恢复 | `modular_recovery.py`；`verify --step20` | 整性 Cauchy 坐标、预算内 D.8 扫描和 D.9 分层恢复；值为 9+4V，精度 101/149。给定包布局不证明 Lemma D.3 的分组结论；不是 D.16 的递归收缩 |
| Appendix D.4.2 的有限日程与碰撞分包 | `modular_packing.py`；`verify --step21` | 精确整数日程及 D.6 平方守卫；邻接重数直方图处理共享端点和平行边的条件期望；统一填充和相等坐标副本位；数值前提与实际结论分别检查。标签与包的资源预算在搜索前核对 |
| Appendix D.4.6 的子网络端口拓扑 | `modular_packing.py` 的 `child_topology` | 保留原端口顺序和孤立顶点，检查有限子图的边数和度数。8192 边匹配的最大子拓扑有 256 条边、度数 4；未构造包张量，也未执行递归收缩 |
| Appendix D.4.5 的有限加权模式 | `modular_modes.py`；`verify --step22` | 一个包的模式表、\(6M/D\) 与指数预算、交叉 Cauchy 余数 8。模式通过后 Step 20 仍恢复 9+4V。拒绝不返回计数。不是包张量，也不是 D.16 |
| Lean 编码、参数和流式计数前沿 | `lean_encoding.py`、`lean_parameters.py`、`lean_sparsify.py`；`verify --step23` | 二进制往返与计数字、首个通过精确整数条件的 M、规范化和优先核心；只保留活动前沿，按原声明 n 累加互不重叠的叶。结构先预检，叶拒绝时不发布部分总和 |
| Lean 三元系数与真实查询 | `lean_ring.py`、`lean_squarefree.py`、`lean_query.py`；`verify --step24` | 实际 Newton digits、独立 U/Z 截断和完整乘积后读出；支持并集乘法与布尔单项式权重求和；原始子句编号、重复和 fixed 条件下的接受/标记违反计数。单位行给出 CNF #SAT，原生叶已接到同一前沿 |

Step 24 的完整查询没有布尔赋值扫描，但仍存储平方自由支持表，可能指数增长。当前 CNF 叶对应 J=1；一般 t*J 稠密环还没有接通流式节点求值、修正对和 Hermite 恢复。接受查询是历史行的公共零点数，标记查询是接受集中的单个子句违反数，它们不能直接替换原 CNF #SAT。源声明、恒等式和预算见 [LEAN_NATIVE.md](LEAN_NATIVE.md)。

`count_placed_formula` 使用安全次数上界 \(d_H=M(R-1)\)，取 \(q=\lfloor d_H/Q\rfloor\)、\(J=q+1\)，因此 \(QJ>d_H\)。它只调用局部扫描，不调用 `brute_count`、`eval_H_at` 或 `eval_H_hasse_jet`。这些全空间函数用于独立验证。

测量矩阵的冗余行通过选取原来的独立整数行删除，保留整数性和绝对行和上界。每行的测量范围收紧到负系数之和与正系数之和之间，仍覆盖所有布尔输入。该范围是参考文献 `[-D,D]` 的子集。

完整算子的结构化描述支持有限个条目查询和受预算限制的测量模态枚举。`place_sparse_formula` 已将这些条目和测量描述接入占用坐标上的放置，公共计数入口同时接受稠密和稀疏放置。预算覆盖条目查询、测量候选、消元和枚举。Step 16 会计算 Lemma 4.8 的闭式，并在预算内执行一个通过 Lemma 4.7 全部数值前提的有限元组。该闭式在这些 n 上于形成 D 之前失败；有限元组的执行不是渐近界的证明。

分量入口采用 r=2、Q=2、R=1 的有限局部设计，直接核验其恢复和核条件；这些小分量不满足 Lemma 4.7 的大 n 前提。局部设计的准确性与 Lemma 4.7 的渐近存在性是不同的核验。每个分量都重新放置，不能假设全局放置的 H(t) 自动按原 CNF 分量分解。

`count_placed_formula_by_support` 检查实际替换多项式在所有次数上的非零系数，并合并共享变量的因子。支持集 V_i 互不相交时，H(T)=2^free·∏H_i(T)。每个局部 H_i 在 T=1 计算覆盖其次数上界的完整 Taylor 系数，再恢复 H_i(0)。原 CNF 的独立分量可能被替换耦合，测试包含这种情况以及非恒定 H(T) 的全系数对照。预算按不使用缓存的 ∑2^|V_i| 计算；复用只减少实际扫描。

`count_placed_formula_by_contraction` 使用同一批替换多项式。它先为每个连通支持做最小度消除计划，再构造数值表。消除一个布尔变量保持和积多项式。局部喷射覆盖该分量的次数上界，Hermite 求值恢复 H_i(0)，再乘以 2^free。槽数统计不含输入多项式存储、规划元数据、算术临时量和整数位长。该后端不使用复用。DIMACS 的 `--backend tensor` 仍先按原 CNF 独立分量重新放置，再对每个放置调用这一收缩；`--backend common` 继续使用共享测量。默认消除桶宽度为 12。

## 已记录的核验

- 1024 个独立单元子句，r=2、Q=2、R=1：精确条件期望选择完成，重数 64≤128，最大组 512，设计及 T=0 恢复通过。没有枚举组外的 \(2^{512}\) 个布尔点。这次核验没有计数该大公式。
- 同一参数的尾概率为 7/16，超过 τ/(rb)=1/4，因此 Lemma 4.7 的全部前提不成立。
- 16384 个独立单元子句，r=2、Q=3、R=2、τ=1/2：尾概率 5/32≤1/4，全部数值前提成立。稀疏放置、设计、T=0 恢复、公共分母和范数检查通过；重数 256、最大组 8192，各节点组秩均为 (2048,2048)。存储 59392 个多项式条目、14336 个测量条目，稠密多项式位置数为 536870912。模态坐标查询 163840 次、拉回候选 34816 个、消元条目 10240 个。共享测量计数计划 q=5461、J=5462，其布尔向量上界有 36697 位；该路线在枚举前拒绝，实际扫描为零。
- 同一 16384 变量全局放置按实际替换支持集分为 10240 个分量，最大支持集 4 个变量、局部次数上界 4，总次数上界 14336。启用复用时，13 次局部求值、10227 次复用，实际扫描 96 个局部布尔向量；关闭复用时，10240 次局部求值、零复用，实际扫描 49152 个，恰好达到 ∑2^|V_i| 的计划上界。两者均恢复精确计数 1，并与独立交替极性单元子句的闭式核对。测试禁用原 CNF 分量计数、全空间求和与共享测量后端；将扫描预算收紧到 49151 时，枚举前拒绝且不返回部分计数。
- 四个原 Step 3 放置经共享测量计数得到 1、1、7、7；各节点的全部请求 Taylor 系数与独立全空间求和相同。另有二阶 Möbius 恢复、负文字、空公式、空子句、自由变量和预算拒绝的测试。
- 稠密边长 4、9、16、64、81 上的所有 M0 条目与替换多项式条目均与结构化查询一致。R=32、Q=192 的测试检查恢复、支持、样本消零及 f 在指定区间的值。
- `examples/count_demo.cnf` 声明 8 个变量，两个独立三文字子句、两个自由变量；实际计数入口返回 196，与暴力对照一致。
- 10923 个独立三文字子句、32769 个变量：实际计数为 7^10923，与独立子句闭式一致，1 次局部代数求和和 10922 次复用。增加两个声明自由变量得到 4·7^10923。该输入也通过 Q=3、R=2 的数值前提，但未执行其全局放置或全局求和。该结果是独立分量计数，不是一般连通实例的界。
- 13 个变量的连通路径：同一放置的 Step 14 支持集路线为 `REFUSED`。Step 15 计数为 610（`COUNTED`），与无相邻零的布尔串闭式一致；最大桶宽 2，局部布尔向量 98，因子求值 48，喷射乘法 190，峰值表项 50。`examples/connected_path.cnf` 经 `--backend tensor` 得到同一计数。
- 1024 个变量的连通路径再加 2 个声明自由变量：计数为 4·Fibonacci(1026)（`COUNTED`），最大支持 1024，自由变量 2，桶宽 2，局部布尔向量 8186，峰值表项 4094。DIMACS `--backend tensor` 且允许该分量大小时，计数和扫描次数相同。闭式来自无相邻零布尔串的递推，不是全空间枚举。
- 两个耦合的替换因子：计数为 1（`COUNTED`）。在 T=1 恢复的喷射为 7/4、1/2、−1/4，与独立全空间 Hasse 喷射逐项一致（系数对照为 `CHECKED`）。这次系数对照使用了全空间求和；长路径计数测试会禁用该神谕。
- 一个宽度为 13 的因子：`REFUSED`。因子支持上界 13 超过预算 12，布尔向量扫描为 0，不建数值表，不返回计数。诱导出的过宽消除桶，以及后面才超预算的分量，同样在任何数值表之前拒绝。
- Lemma 4.8 在 16384 个交替单元子句上给出 r=0、Q=0、R=16，λ=14、ν=4、β=1。初步守卫失败，τ 和 D 都没有形成，也没有分配矩阵。标签为 `CHECKED`。λ=2¹⁹ 的整数是 ν=19、r=2、Q=215、R=2900，因为 Q≤R，同样不形成 D。
- 8192 个交替单元子句：预算内没有任何元组通过 Lemma 4.7。r=2、Q=2、R=1 的尾概率 7/16 大于每个 τ≤1/2；r=2、Q=3 需要 2(b+w)D²=16384>M。标签为 `CHECKED`。没有开始放置。
- 16384 个交替单元子句：选择器给出 r=2、Q=3、R=2、τ=5/16，区间 {0,1}，尾概率 5/32，恰好等于 τ/(rb)。这比手选的 τ=1/2 更紧，测量允许量 2τn 从 16384 降到 10240。该元组不满足 R≥32 或 6R≤Q。标签为 `CHECKED`。它不是 Lemma 4.8 的闭式。
- 同一选择上的稀疏放置复现 Step 12 的存储：多项式条目 59392、测量条目 14336、稠密多项式位置数 536870912、重数 256、最大组 8192、三个节点的组秩均为 (2048,2048)。独立测量数 4096，落在 10240 的允许量内。
- 同一放置的三个既有后端：共享测量在枚举前拒绝，q=5461、J=5462，布尔向量上界 36697 位，实际扫描为 0，标签 `REFUSED`。支持集复用计数为 1（`COUNTED`），13 次局部求值、10227 次复用、96 个局部布尔向量。有限 Taylor 消元计数为 1（`COUNTED`），桶宽 4、局部布尔向量 176128、喷射乘法 262144、峰值表项 72。候选试验预算为 1 时，在日程搜索之前拒绝。稠密 `build_auxiliary_slice(2, 192, 32)` 仍在分配前拒绝，边长仍报告为 2³⁸⁴。
- 这些有限元组的选择和运行不是 Lemma 4.8 的渐近时间或空间界，也不是 Sections 5–6 或 Appendix D。Step 15 仍是有限 Taylor 消元。独立分量计数和不相交支持集计数不是一般连通实例的界。
- Lemma 5.1 在 \(\ell=1,2,3,4\) 上对 \(y=-6,\ldots,12\) 成立：偶数同余于 1，奇数同余于 0（`CHECKED`）。\(\ell=2\) 的多项式是 \(1-3Y^2+2Y^3\)，除法结果是 \((1-Y)+(2Y-1)h(Y)\)。两个整数上的截断乘积与文中三项一致，被丢弃的 \(h(y_1)h(y_2)\) 项能被 4 整除（`CHECKED`）。
- \(n=4\)、\(k=2\)、\(c=1\) 的系统 \(x_0x_1+x_2=0\)、\(x_3+1=0\)：截断指示器模 \(2^5\) 的最小非负剩余是 4，16 个布尔点上逐点等于公共零点指示器（`COUNTED`）。空方程列表给出 \(2^5=32\)，常数方程 \(1=0\) 给出 0，两者的布尔扫描都是 0。
- 同一 4 变量系统在项数预算 2、系数位数预算 1 或布尔立方预算 1 下为 `REFUSED`。计数是空的，`boolean_vectors=0`，项表没有分配。系数位数和立方预算也没有构造出被返回的部分系数。
- 这些检查不是 Theorem 2.2 的位运算界，也不是 Section 6 的方程选择。Step 17 没有完成 Theorem 1.1。
- \(n=4\)、\(k=2\) 的系统 \(x_0x_1+x_2=0\)、\(x_3+1=0\)：输入次数上界 \(d=20\)，\(d_H=20\)，出现次数上限 \(b_0=4\)，重变量集为空，四个中性范围（`CHECKED`）。同一系统在 \(Q=21\)、\(R=1\)、斜率 \(1/21\) 的有限设计上，\(S_\ell(H(0))\) 的 \(V\) 系数是 \((-2,6,0,-192,2496)\)，未约化和 2308，模 \(2^5\) 的最小非负剩余是 **4**（`COUNTED`）。它与独立的 \(\mathbb F_2\) 扫描以及不经替换行的布尔立方求和一致。实际扫描 399 个布尔向量，计划上界 714，因子求值 2016 次。\(T=1\) 处的和不是 4。该系统上 Lemma 4.8 为 \(r=1\)、\(Q=0\)、\(R=10\)，\(D\) 没有形成。
- 五个 \(x_0=0\) 与一个 \(x_0+x_1=0\)：\(M=11\)，\(b_0=6\)，重变量是 \(x_0\)。两个分支的计数是 \((4,0)\)，合计 **4**（`COUNTED`），实际扫描 110。限制后的分支有 \(n'=3\)。Lemma 4.8 给出 \(r=1\)、\(Q=0\)、\(R=8\)，初步守卫失败，\(\tau\) 和 \(D\) 都没有形成（`CHECKED`）。求和用的是有限设计，不是这篇证明里的穷举回退，也不是 Section 4 的辅助阵列。
- 同一 4 变量系统在布尔向量预算 1 下为 `REFUSED`：上界 714，扫描 0，喷射没有分配，系数没有构造。六个方程的系统在分支预算 1 下为 `REFUSED`：2 个分支都没有进入循环。系数位数预算 1 在构造 \(p_\ell\) 时拒绝，喷射没有分配。无效核条件在枚举之前返回，计数为空。
- 这些检查不是 Theorem 2.2 的位运算界，也不是 Section 6 的方程选择或 Appendix D。Step 18 没有完成 Theorem 1.1。
- \(n=m=3\)、\(B=1\)：\(\mu=2\)、\(P=80000\)、\(D_0=2560000\)。\(q_G=1\) 且 \(h\le 0\)，\(A(0)\) 是两行偶校验，行重 3，行数 \(2\le 3n\)。每个非零核向量的权重是 2，不低于 \(m/P\)（`CHECKED`）。扫描 8 个向量。
- \(n=m=2048\)、\(B=1\) 的局部码搜索是 524288 比特，在分配矩阵之前拒绝，扫描为零（`REFUSED`）。
- 两个变量、两条单元子句、列表长度 1：6 轮、24 次计数，最终 \(N(A)=1=\#F\)（`COUNTED`）。`count_common_zeros` 在同一系统上返回 1。论文规模 \(P=80000\)、候选预算 1、立方预算 1 都不返回计数，也不枚举（`REFUSED`）。
- 3 变量 2-CNF：叶计数 \((3,1)\)，和为 4，与暴力 \(\#F\) 一致（`COUNTED`）。声明 \(n\) 不变，出现次数不超过 \(B=6\)。节点预算 1 不返回叶，也不返回计数（`REFUSED`）。
- Proposition 6.4 在 \(n=20\)：\(\lambda_n=5\)、\(\nu_n=3\)、\(s=0\)，因此是穷举情形（`CHECKED`）。\(s=2\) 只是 \(\lambda=521284\)、\(\nu=19\) 上的整数见证。\(k=2\)、\(\sigma=1\) 时 Appendix A 的 2 的幂是 \(M=64\)。不声称 \(O_k(\rho(n))\)。
- 这些检查不是 Appendix D，也不是 Theorem 1.1。Step 19 没有完成 Theorem 1.1。

Step 20 的十四项检查覆盖模系数环、CNF 网络计数 196、相等树计数 4、空公式与空子句、创建时精度、Cauchy 坐标、两种行张量的恢复、两层精度，以及网络扫描、恢复扫描、精度、喷射存储和剩余域点数不足的拒绝。模网络的恢复结果为 9+4V 或 8+5V（模 32，V²=0），属于 `CHECKED`；CNF 的计数属于 `COUNTED`。拒绝不返回值，扫描为零。详细接口和对应公式见 [APPENDIX_D.md](APPENDIX_D.md)。

Step 21 的十四项检查覆盖整数日程、对数见证、D.6 平方守卫、共享端点与平行边的精确期望、统一填充、多坐标副本位、实际子拓扑、自动布局恢复，以及标签、跟踪、分数更新、填充和子图预算拒绝。128 边匹配在 `r=2,a=1,Δ=4` 上通过 D.6，得到 64 个包，但不缩减到 `M/8`。8192 边匹配在 `r=64,a=1,Δ=4` 上默认拒绝；显式有限模式生成 4096 个包，最大子拓扑 256 条边、度数 4，张量求值为零。后者实际结论为 `CHECKED`，但 D.6 第一条件要求 `M≥134217728`，依旧不成立，且该元组不是 D.4 日程成员。自动布局接入 Step 20 后，两条重复单元子句加自由变量的精确计数为 2（`COUNTED`），扫描 6656 个端口赋值对，求值 33280 次条目。详见 [PACKET_PACKING.md](PACKET_PACKING.md)。

Step 22 的十三项检查覆盖一个包的模式表、空的局部包、指数元组、加权支持拒绝门、交叉 Cauchy 余数 **8**、模式通过后的 Step 20 恢复 **9+4V**，以及模式个数、\(6M/D\)、元组上界、部分选择上界、不合格次数、支持记录预算和超预算指数的拒绝。拒绝的 `count` 为 null。包张量求值为 0。详见 [WEIGHTED_MODES.md](WEIGHTED_MODES.md)。

Step 25 的十五项检查覆盖日程整数公式、相对距离恒等式、局部距离与 Walsh 的首选见证、槽位行再生、证书算术、接受/拒绝证书、矩阵与候选预算拒绝、候选流与贪心暴力轨迹、零质量模型集、流式 hash 叶、三类拒绝和终端情形。

Step 26 的十二项检查覆盖早期对数与阈值、截断上取整与两条上界、边界排除、单位行/空历史/显式历史/标记/固定坐标的秩和恒等式、双高坐标四秩、rank/query/代数预算拒绝、直接与路由贪心一致（hub 同为 129）以及流式路由叶与终端情形。

Step 27 的十二项检查覆盖日程整数复现、守卫边界（K=2^24 拒绝、K=2^25 通过）与 `Q*J > degree` 恒等、独立乘法分解器核对的首个不可约见证、环运算与可逆元（牛顿提升）、`3^t` 以内唯一自然数读回、节点互异与单位差、精确 Taylor 喷射与 Newton–Hermite 精确恢复/条件数失配、双字长编码与搜索拒绝。

Step 28 的九项检查覆盖修正对计数与 `fold` 顺序枚举、差分系数对支撑表、有符号探针的流式=直接=暴力（116）、修正对/赋值/项数预算拒绝、节点流几何（`Q*J*d` 字）与解码错误、唯一读回范围、星形/双星/hub 模型计数（9/25/129）、一阶终端情形与流字预算拒绝。

Step 29 的九项检查覆盖度数界与阶公式、标记和逐系数一致（含标记与固定坐标）、接受/标记计数 8/4、单位行叶 9/25/129/196、支持表预算拒绝下的替换计数、三类预算拒绝、计数字宽、终端情形，以及 `--leaf-backend stream` 前沿一致与共享预算拒绝。

当前完整报告 `benches/reference_alignment.json` 覆盖 Steps 9–29，新增 Step 29 的九行：度数界、标记和一致、无支持表替换、预算拒绝、终端情形与流式叶后端。`artifacts/step15-verification.json` 至 `artifacts/step29-verification.json` 可分别用对应 `--stepN --output` 重新执行并保存；报告版本与包版本一致。

版本 0.29.0 的完整 pytest 为 **555/555**，完整参考检查为 **185/185**；Steps 15–29 的单步报告也已实际重新执行生成。

报告区分 `COUNTED`、`CHECKED`、`REFUSED`，保留全部完成标志为 false，不记录运行时间，大计数以精确十六进制保存。`--output` 只记录实际执行的检查，失败或空 suite 不会标为全部通过。静态 Lean 审计独立保存在 `artifacts/lean-reference-audit.json`；它不替代计数检查或 Lean 内核核验。

命令：

```bash
python -m pip install -e ".[dev]"
pytest -q
python -m counting_sat_algebra.verify --reference --output benches/reference_alignment.json
python -m counting_sat_algebra.verify --step14
python -m counting_sat_algebra.verify --step15
python -m counting_sat_algebra.verify --step16
python -m counting_sat_algebra.verify --step17
python -m counting_sat_algebra.verify --step18
python -m counting_sat_algebra.verify --step19
python -m counting_sat_algebra.verify --step20
python -m counting_sat_algebra.verify --step21
python -m counting_sat_algebra.verify --step22
python -m counting_sat_algebra.verify --step23
python -m counting_sat_algebra.verify --step24
python -m counting_sat_algebra.verify --step25
python -m counting_sat_algebra.verify --step26
python -m counting_sat_algebra.verify --step27
python -m counting_sat_algebra.verify --step28
python -m counting_sat_algebra.verify --step29
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra.lean_reference --output artifacts/lean-reference-audit.json
python -m counting_sat_algebra count examples/count_demo.cnf --output benches/count_demo.json
python -m counting_sat_algebra count examples/count_demo.cnf --backend lean --leaf-backend native --json
python -m counting_sat_algebra count examples/connected_path.cnf --backend tensor --output benches/connected_path.json
```

`--step9` 至 `--step29` 可分别运行二十一组参考检查，也支持 `--output` 保存所选结果。`--step8` 保留原布尔预算模式。默认 `schedule_placement` 使用线性代数检查。默认 DIMACS 后端是 `common`；`tensor` 提供有限 Taylor 消元，`lean` 提供流式外层及 `common|tensor|native` 叶。`--quick` 运行 Steps 0–8 的小网格和 Steps 9–29 的同一组有限检查。Steps 17–22 的模网络、分包等仍由独立 API 提供。

## 尚未完成的条件

1. （Step 25 已完成）构造并核对 Lean 的实际三元稀疏码，按游标再生候选行，保留距离、宽度和出现次数约束。
2. （Step 25 已完成）用真实标记查询累加候选行的违反质量，接通流式贪心选择，并独立执行最终历史接受查询；`hash` 叶已使用该程序。
3. （Step 26 已完成）接通原声明 n 上的高变量赋值循环：早期阈值、截断、掩码、逐秩条件化累加与 `--high-loop` 路由；完整参数前端与出现次数守卫随整机精化推进。
4. （Step 27 已完成环与节点部分）构造计数扩张环、互异插值节点并建立 QJ>degree；实际分组放置进 ABI 仍待接通。
5. （Step 28–29 已完成）修正对分解、逐节点 J 系数流、Newton–Hermite 读回与无支持表的逐点叶（`--leaf-backend stream`）已落地；完整机器次数约简与放置随整机推进。
6. 对整条 Python 程序核对精确位时间、所有前缀的工作空间、完整计数字和成功/失败路径；将操作统计与 Lean 原始机器资源建立精化关系。
7. 在库锁定的 Lean 环境中实际编译入口、保存公理依赖输出。新增库已有证明源码；尚未完成的是本地内核核验及 Python 对该程序的移植与精化。

Appendix D 的 D.4.6 包张量、D.4.7 子递归和 D.3 编码单列为后续路线，以达到论文 Theorem 1.1 的更强界。当前的有限参数成员、独立分量、支持集收缩、Taylor 消元与模局部恢复都不能建立一般输入的渐近结论。
