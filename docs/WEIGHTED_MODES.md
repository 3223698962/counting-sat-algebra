# Step 22：有限加权模式

版本 0.23.0 在 `modular_modes.py` 中实现 Appendix D.4.5 的有限模式表。对象是一个包、一个活跃节点 \(t\) 和显式列槽。交叉槽的截断几何级数与 Step 20 的 Cauchy 局部多项式核对。模式表通过之后，才可以把同一布局交给 Step 20 的恢复。

这一阶段不求值 D.4.6 的包张量，不构造 `cache`、`rows` 或 `origins`，也不做 D.4.7 的递归收缩。`recursive_contraction`、`theorem_1_1_complete` 和 `lean_proof_complete` 保持 false。SqrtLog 的加权支持条件和这里的包模式不是同一个对象；对照见 [LEAN_SPEC_CHECKLIST.md](LEAN_SPEC_CHECKLIST.md)。

## 预算

`build_mode_table` 在分配模式表和指数元组之前检查：

- 包长 \(L=2^a\)、精度和实列数不超过预算；
- 模式个数不超过 `max_modes`；
- 同一权重的瓦片数 \(2^{v_2(t)+1}\) 满足 \(2^{v_2(t)+1} D \le 6M\)，其中 \(D=L/2^{v_2(t)}\)；
- 无约束指数积不超过 `max_tuples`；
- 部分选择上界（元组数乘变量数再加一）不超过 `max_partial_choices`。

任一检查失败时返回 `REFUSED`：没有模式表，没有元组，`count` 为 `None`。恢复函数在模式被拒绝时不会调用 Step 20。

次数记录使用 \(\deg X \le \deg U + \deg Z\)。不满足时 `admit_support` 在存储前拒绝。这只是拒绝门，不是 Lean 多项式定理。

## 小实例

\(a=2\)、\(t=1\)、\(M=2\)、精度 4、符号 `(0,1)` 给出 6 个模式，权重 `1,2,3` 各两份。满足权重和的指数元组有 **18** 个，枚举访问 **67** 个部分选择。分配前使用的无约束乘积上界是 256；上界超过预算时直接拒绝，不先数出这 18 个元组。槽 `(0,2)` 在 \(\zeta=1\) 的 Cauchy 余数为 **8**（`COUNTED`）。系数的 2-adic 赋值是 `(3,4,5)`，分别不低于权重。两条边的网络在模式表通过后由 Step 20 恢复为 **9+4V**（`COUNTED`）；模式本身的张量求值为 0。

```bash
python -m pytest -q tests/test_step22_weighted_modes.py
python -m counting_sat_algebra.verify --step22 --output artifacts/step22-verification.json
```
