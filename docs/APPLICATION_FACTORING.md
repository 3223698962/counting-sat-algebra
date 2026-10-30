# 确定性因子分解：应用说明与结论距离

`references/others/DeterministicFactoring.tex`（Rudich 与 Rousseau，Polytechnique Montréal）给出一个 Lean 证书：对每个 \(k\ge1\) 存在固定确定性图灵机在 \(C_k\,2^{n/k}\) 步内输出 \(\mathrm{minFac}(N)\)，对完整分解至多 \(\log_2 N\) 次调用（定理 `subexponentialRate`；审计语句仅 13 行、只依赖 Mathlib，完整证明约 36,000 行，未随本仓库提供）。`references/others/` 保存本文档所说明的应用论文材料本身：`DeterministicFactoring.tex` 与 arXiv 2610.00600v1 的 PDF/源码压缩包（压缩包内只有同一份 tex）。

## 证明用的是什么引擎

论文自述（“What we use from these papers”）：使用 Claude（Anthropic）2^{o(n)} 3-SAT 论文（[Cla26]；其 Lean 开发在 anthropics/formal-math 仓库的 eth 目录，不在本仓库）的 Sections 3–5——网络分组（Lemma 1）、八平方交叉恒等式（Lemma 2）、网格化重写与递归（Proposition 3、\(e^*=0\)）；**不使用** OpenAI 论文（[Ope26]）的变形与 Hermite 插值路线。它只重新证明所需部分、不导入任何原开发。

因此应用所需的计数引擎与下面三者都不同：

1. 本仓库的 `SqrtLogSpace` 主线（Step 23–26 的有限移植）——另一套计数引擎；
2. Appendix D 的布尔张量网络收缩（主论文 Theorem 1.1 的两条界的证明路线，Step 20–22 的有限部件属于它）；
3. [Cla26] 的网络递归（本仓库尚未实现）。

分解归约本身只用「按最高位前缀计数」：逐比特确定最小因子（`search_eq_minFac`）。对应的条件化计数能力本仓库已具备（`count_ternary_query(fixed=...)` 与 Step 26 的条件化秩循环）。

## 与本仓库部件的对应

| 应用所需对象 | 本仓库现状 |
| --- | --- |
| 按固定前缀计数（条件化） | 已支持：`fixed` 条件化与 `count_cnf_high_loop` 秩循环 |
| 网络值、分组与平衡（Cla26 Lemma 1） | 有限局部：Step 21 碰撞标签、包布局与子拓扑；一般分组搜索未实现 |
| 交叉恒等式与成本分数（Cla26 Lemma 2、八平方） | 有限局部：Step 22 的 D.4.5 加权模式与一条交叉 Cauchy 余数（模 32 见证）；一般模式展开未实现 |
| 递归收缩（Cla26 Proposition 3、\(e^*=0\)） | 未实现：`recursive_contraction=false`；D.4.6 包张量、D.4.7 子网络递归与 D.3 关联编码仍是后续路线 |
| 计数判零（二项式整除、Erdős–Bertrand） | 未实现：仓库只有精确计数与预算拒绝，没有按小模数判零的程序 |
| 三层编译到 Mathlib 图灵机 | 未实现：Python 侧没有编译与精化 |

## 结论能否达到

- 主论文 Theorem 1.1 与因子分解推论都要求完整的计数引擎加资源精化；当前 `theorem_1_1_complete`、`lean_proof_complete`、`recursive_contraction` 全部保持 false，本仓库没有任何一般输入的 \(2^{o(n)}\) 声明。
- 两条证明线各自的完成状态要分开记录：SqrtLogSpace 主线已有 Steps 23–28 的有限见证（含扩张环与流式修正对读回），剩余生产叶切换与整机资源；Appendix D（与因子分解应用同族）已有 Steps 20–22 的有限部件，剩余 D.4.6/D.4.7/D.3 递归。任何一条的「已完成」都不能替代另一条。详见 [RESOURCES_AND_CONCLUSION.md](RESOURCES_AND_CONCLUSION.md)。
- Lean 内核核验在本机不可执行：没有 `elan`/`lean`/`lake` 工具链（网络可达，但 Mathlib 级依赖构建超出单次会话），应用论文的 36,000 行证明也不在本仓库。静态审计（`artifacts/lean-reference-audit.json`）覆盖的是 SqrtLogSpace 参考库，不是本文档描述的应用证明。
- 可推进的下一步（与 [ROADMAP.md](../ROADMAP.md) 的 Appendix D 路线一致）：D.4.6 包张量系数 → D.4.7 子网络递归 → D.3 编码，并用 Step 22 已核对的交叉余数逐层检查乘性缺口；该路线与 SqrtLogSpace 主线的流式 Taylor 叶相互独立。
