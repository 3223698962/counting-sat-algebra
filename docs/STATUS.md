# 现状

当前包版本 **0.22.1**。有限检查与 0.22.0 相同，算法行为没有改变。Steps 0–21 的见证仍在；Theorem 1.1 没有实现，本仓库不否定 ETH，也没有一般的 \(2^{o(n)}\) 算法。

`theorem_1_1_complete`、`lean_proof_complete` 与 `recursive_contraction` 均为 false。预算拒绝发生在枚举或分配之前，不返回部分计数。Lean 主定理 `exact_count_and_sat` 仍以 `sorry` 结束。

## 现在能做什么

- Steps 0–16：有限插值、设计、稀疏放置、独立分量 DIMACS 计数、支持集计数、预算内的 Taylor 消元，以及 Lemma 4.8 的整数日程和一个预算内有限成员。
- Steps 17–19：Section 5 的有限稀疏系统与 \(A_\ell\) 替换求和，Section 6 的有限方程选择和同一声明 \(n\) 上的互不重叠稀疏化。不声称 Proposition 6.4 的 \(O_k(\rho(n))\) 界。
- Steps 20–21：Appendix D 的有限系数环、CNF 网络、Cauchy 局部恢复、碰撞分包和子拓扑。直接叶扫描和端口对扫描都有预算。

活跃路线是 Appendix D 的递归收缩。D.4.5 加权模式、包张量系数、D.4.7 递归收缩和 D.3 编码还没有实现。Step 15 的有限 Taylor 消元不是这条路线。层次说明见 [ARCHITECTURE.md](ARCHITECTURE.md)。已完成步骤的清单见 [steps/completed.md](steps/completed.md)，逐步叙述见 [PROGRESS.md](PROGRESS.md) 与 [steps/through-step21.md](steps/through-step21.md)。

## 入口

```bash
python -m pip install -e ".[dev]"
pytest -q
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra.verify --reference --output benches/reference_alignment.json
python -m counting_sat_algebra.verify --step21 --output artifacts/step21-verification.json
python -m counting_sat_algebra count examples/count_demo.cnf --output benches/count_demo.json
python -m counting_sat_algebra count examples/connected_path.cnf --backend tensor --output benches/connected_path.json
```

`--reference` 运行 Steps 9–21。`--quick` 覆盖 Steps 0–21 的较小网格，并包含 Steps 16–21 的同一组有限检查。`--step9` 至 `--step21` 可单独运行，也支持 `--output`。DIMACS `count` 只使用 `--backend common|tensor`；Steps 17–21 仍是独立 Python API。标签只有 `COUNTED`、`CHECKED`、`REFUSED`。检查失败时报告保留已执行的行，并把完成标志留在 false。

## 尚未完成

D.4.5 加权模式、D.4.6 包张量系数与不可变创建时描述、D.4.7 递归收缩及其不变量、D.3 编码、精确位运算与空间界，以及去掉 `sorry` 的 Lean 证明。详见 [ROADMAP.md](../ROADMAP.md)。
