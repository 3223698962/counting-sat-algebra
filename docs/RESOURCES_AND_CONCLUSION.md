# 结论可达性与资源审计

本页汇总 Steps 23–28 实际记录的资源收据、尚未完成的精化义务，以及三条证明路线的结论可达性评估。所有完成标志（`theorem_1_1_complete`、`lean_proof_complete`、`recursive_contraction`）保持 false；本页不主张任何一般输入的渐近结论。

## 已记录的资源收据（有限、可重放）

| 部件 | 收据 |
| --- | --- |
| 流式前沿 | 前沿位数、节点/叶数；叶拒绝不返回部分总计 |
| hash 叶 | 候选数、轮数、标记查询数、接受查询数、历史行数、质量轨迹（初始/最终） |
| native 查询 | 精度、模 `3^t`、marker 系数和、单项式乘积、标量操作、活跃系数、逻辑编码位数 |
| 高变量循环 | 阈值、截断、高变量基数、秩数（每查询 `2^|H|`）、条件化查询数、累加位宽 |
| 扩张环 | 日程整数、候选/除法检查/环操作数、相对距离与核下界、`Q*J − m` 间距 |
| 修正对 | 对数、赋值数、流字数、前缀线性长度、`fold` 顺序与剪枝计数 |

这些计数是 Python 侧的操作与逻辑编码统计，不是 Lean 原始机器的步数，也不是已访问磁带格数；两者之间的精化证明尚未建立（ROADMAP 第 6 项）。

## 三条证明路线与结论距离

1. **SqrtLogSpace 主线**（本仓库 Step 23–28）：语义核心已有限移植并逐项见证；缺 (a) 生产叶切换到流式路径（立即 `U` 特化、完整修正对再生）、(b) 放置进 896 栈 ABI、(c) Python 位时间/全前缀空间精化。时间指数的尺度（含 `(log log n)²/√(log n)`）见 `SqrtLog/StatementScale.lean` 与 `SqrtLogSpace/Asymptotics.lean`。
2. **Appendix D 路线**（论文 Theorem 1.1 的更强界）：Step 20–22 落地系数环、网络、Cauchy 变形、分包、加权模式与一条交叉余数；D.4.6 包张量、D.4.7 递归收缩与 D.3 编码仍是明确开放，`recursive_contraction=false`。
3. **应用论文路线**（`DeterministicFactoring.tex`，使用 Cla26 Sections 3–5 引擎）：归约（按最高位前缀计数）已被本仓库的条件化查询覆盖；引擎部分（网络分组、八平方交叉、递归）不在本仓库，映射与结论距离见 [APPLICATION_FACTORING.md](APPLICATION_FACTORING.md)。

三条路线互不替代：任一路线的有限见证都不构成另一条的完成。

## Lean 内核核验（ROADMAP 第 7 项）

本环境没有 `elan`/`lean`/`lake` 工具链（已核查 PATH 与用户目录），网络可达但 Mathlib 级依赖构建超出单次会话范围；应用论文的约 36,000 行证明亦不在本仓库。静态审计（`artifacts/lean-reference-audit.json`：3969 个模块、无缺失本地导入、无代码级占位符或本地公理）继续作为唯一已执行的核验，`kernel_verified=false`。核验步骤（`lake build …`、保存 `#print axioms` 输出）记录在 [references/lean/README.md](../references/lean/README.md)。

## 怎样算“达到结论”

只有同时满足以下条件，才可把任一路线记为完成：

1. 生产入口在一般输入上运行该路线的完整程序（不接受有限参数替代）；
2. 位时间与全前缀空间界经独立分析或精化建立；
3. 对应 Lean 入口在内核中编译并保存公理依赖输出。

当前全部条件均未满足；本仓库提供的是逐步可重放的有限证据链与明确的开放义务清单。
