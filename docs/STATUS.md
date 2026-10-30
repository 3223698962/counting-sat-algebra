# 现状

当前包版本 **0.28.0**。项目主线已经按新加入的完整 Lean 参考程序调整为 `SqrtLogSpace` 的精确计数移植。新库的证明源码与旧单文件脚手架分别记录；本项目还没有完成全部 Python 叶程序、一般输入的位时间和空间分析，也没有执行新库的 Lean 内核编译核验。

## 现在能做什么

- Steps 0–22：原有精确插值、设计、稀疏放置、DIMACS 分量计数、有限 Taylor 消元、Section 5–6 对象，以及有限 Appendix D 模环、网络、分包、模式和恢复。现有 API 保留。
- Step 23：Lean 二进制语法和计数字、首个通过条件的二进参数、规范化、优先核心选择，以及仅保留活动前沿的互不重叠流式计数。
- Step 24：模 `3^t` 的实际 Newton 指示系数、共享标记截断、平方自由乘法与布尔求和、接受/标记违反查询、原 n 上的条件化，以及已接通前沿的 `native` 计数叶。精度、活跃系数、逻辑位数和整棵前沿的算术预算均受限制。
- Step 25：按日程整数构造三元稀疏码并逐行再生，用真实标记查询累加全部子句的违反质量做贪心选择，零质量后单独执行接受查询。`hash` 计数叶已接通前沿；候选深度、候选数、轮数与标记查询预算均受限制。
- Step 26：计算早期高出现阈值与截断，给出环境次序高变量掩码和两条整数上界；逐秩流式执行条件化查询并精确累加，`--high-loop` 可把 hash 叶的全部查询路由到该循环。秩数与代数预算受限制。
- Step 27：按 `LeafParametersProgram` 日程整数（L、w、s、r、Q、R、t、d、J、N）构造模 `3^t` 的计数扩张环：完整枚举搜索不可约单首提升，固定互异节点，核对 `Q*J > degree`、`card·(cut+1) ≤ totalScope` 型上界与唯一自然数读回；精确 Taylor 喷射与 Newton–Hermite 重建（除节点差外不做除法）通过独立见证。
- Step 28：按 `fold` 顺序流式枚举修正对 `(E,U)`（`|E| < J`、`U ⊆ E`，预算耗尽即剪枝，不保存完成对表），用布尔有限差分恒等式精确求和；逐节点输出 `J` 个 2t 位系数字的扁平流（`Q*J*d` 字，含零字），再用 Newton–Hermite 读回唯一计数。修正对、赋值、项数、喷射阶与流字预算受限制。

示例 CNF 的新后端计数为 **196**；4 变量星形公式流式计数为 **9**；三个独立分支的 9 变量公式计数为 **125**。三元历史行接受/标记计数为 **8/4**，共享标记高阶系数参与计算。`hash` 叶直接贪心在星形上用 3 轮、84 条标记查询把质量 12→4→1→0，与逐赋值暴力复现的轨迹逐步一致；流式 hash 叶给出同样的 196、9 和六变量双星 **25**。直接贪心与 `--high-loop` 路由在 hub 公式（8 变量、7 个子句）上同为 **129**；两 hub 公式的秩数为 4。扩张环在 (t,d)=(2,2) 上由 `X²+1` 定义，节点 {1,2} 的差为可逆元，三次目标在 2×2 条件数下精确重建；日程守卫在 K=2^24 边界处因 `4R > Q` 拒绝、在 K=2^25 处通过。修正对流式求和与直接权重和及暴力逐点求和一致（有符号实例 116）；模型指示多项式给出星形 9、双星 25、hub 129。新后端的布尔赋值扫描次数均为 **0**；平方自由支持表仍可能指数增长，不能据此得到一般渐近界。

本轮完整 pytest 的 **529 项测试**全部通过；Steps 9–28 的完整参考验证 **176/176** 通过。`benches/reference_alignment.json` 和 Steps 15–28 的各份 artifact 均通过实际重新执行生成，版本为 0.28.0。

## 入口

```bash
python -m pytest -q
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra.verify --reference --output benches/reference_alignment.json
python -m counting_sat_algebra.verify --step23 --output artifacts/step23-verification.json
python -m counting_sat_algebra.verify --step24 --output artifacts/step24-verification.json
python -m counting_sat_algebra.verify --step25 --output artifacts/step25-verification.json
python -m counting_sat_algebra.verify --step26 --output artifacts/step26-verification.json
python -m counting_sat_algebra.verify --step27 --output artifacts/step27-verification.json
python -m counting_sat_algebra.verify --step28 --output artifacts/step28-verification.json
python -m counting_sat_algebra count examples/count_demo.cnf --backend lean --leaf-backend native --json
python -m counting_sat_algebra count examples/count_demo.cnf --backend lean --leaf-backend hash --json
python -m counting_sat_algebra count examples/count_demo.cnf --backend lean --leaf-backend hash --high-loop --json
python -m counting_sat_algebra.lean_reference --output artifacts/lean-reference-audit.json
```

`--reference` 运行 Steps 9–28。`--quick` 还包含 Steps 0–8 的较小网格。`--step9` 至 `--step28` 支持单独执行与保存报告。`count` 的 `common`/`tensor` 原入口继续可用，`lean` 提供流式外层及 `common|tensor|native|hash` 四个有限叶后端。

标签为 `COUNTED`、`CHECKED`、`REFUSED`。拒绝发生在对应受限操作或分配之前；若先前叶已经完成，其工作会记录，但未完成的总计数和 SAT 答案均为空。

## 证明与完成状态

静态审计核对了 **3969 个本地模块**：本地导入完整，没有代码级 `sorry`/`admit` 和 `axiom` 声明，所列入口包含证明体或证明项。`kernel_verified=false` 表示尚未用 Lean 编译器核验，并不把静态审计视为证明。新库的时间尺度与 Appendix D 的更强目标见 [REFERENCE_ALIGNMENT.md](REFERENCE_ALIGNMENT.md)。

`theorem_1_1_complete=false`：本项目尚未实现论文 Theorem 1.1 所要求的一般输入和完整复杂度结论。`lean_proof_complete=false`：尚未完成本地 Lean 核验与 Python 程序精化。`recursive_contraction=false`：Appendix D 的递归仍未执行。

剩余依赖已列在 [ROADMAP.md](../ROADMAP.md)：把生产叶切换到流式路径、实际放置、完整机器资源和内核核验；资源收据与三条路线的结论可达性评估见 [RESOURCES_AND_CONCLUSION.md](RESOURCES_AND_CONCLUSION.md)。实现范围见 [LEAN_NATIVE.md](LEAN_NATIVE.md)，源声明义务见 [LEAN_SPEC_CHECKLIST.md](LEAN_SPEC_CHECKLIST.md)。
