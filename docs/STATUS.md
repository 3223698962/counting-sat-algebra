# 现状

目标是实现 `references/` 中的确定性精确计数算法，并达到可核验的时间、空间结论。当前版本 **0.14.0** 已推进到 Step 13；完整定理仍未完成。

## 本次新增

- **Step 9：** 按 Lemma 4.5 精确计算条件期望，覆盖共享变量与共享因子的相关性。每次选择、初始势及最终势均以 Fraction 核对。设计检查使用缓存的精确行空间，不再默认枚举布尔组。1024 变量的实例完成放置和恢复检查。
- **Step 10：** 实现 Lemma 3.2 的流式共享测量局部求和、Lemma 3.3 的固定替换和式 (3.13) 的 Taylor 系数恢复，再做 Hermite 插值。四个原小实例的计数及所有请求系数与独立对照一致。
- **Step 11：** 通过交换轨道的谱系数实现结构化 M0 和替换条目查询。小切片的所有条目一致；完整参数 R=32、Q=192 的有限条目可查询。
- **Step 12：** 将结构化条目查询接入稀疏放置、测量构造和精确核检查。16384 个独立单元子句在 r=2、Q=3、R=2、τ=1/2 上通过全部 Lemma 4.7 数值前提及设计检查；多项式条目 59392、测量条目 14336，稠密多项式位置数 536870912。全局求和在枚举前因预算拒绝，扫描为零。
- **Step 13：** 新增 DIMACS 输入和实际 CNF 计数命令。独立连通分量各自重新放置，以代数后端求局部计数，再乘以声明自由变量因子。8 变量示例得到 196；32769 变量的 10923 个独立三文字子句得到 7^10923，执行 1 次局部计数、10922 次复用。
- **回归：** 163 项 pytest 全部通过，`--reference` 的 14 项检查全部通过。保存报告为 `benches/reference_alignment.json` 和 `benches/count_demo.json`。

## 运行

```bash
python -m pip install -e ".[dev]"
pytest -q
python -m counting_sat_algebra.verify --reference
python -m counting_sat_algebra.verify --quick
python -m counting_sat_algebra count examples/count_demo.cnf --output benches/count_demo.json
```

`--step9` 至 `--step13` 分别检查新步骤。`--step8` 保留原布尔预算拒绝路径；日程 API 默认采用线性代数检查。实际计数命令支持 `--json`、`--no-memo` 和显式分量/枚举预算，拒绝时不返回部分计数。

## 尚未完成

满足全部条件的渐近参数族自动选择与执行、Sections 5–6 的输入约化与方程构造、Appendix D 的张量收缩、完整位复杂度和空间分析、以及 Lean 主定理的实际证明。有限参数的稀疏放置已接入公共计数入口，但一般大型连通分量仍会触发预算拒绝；独立分量上的大计数不能代替一般输入的时间界。

参考文献 Sections 3–6 与 Lean 范围说明的时间界比 Theorem 1.1 的 Appendix D 时间界弱。要达到用户要求的参考结论，必须区分并完成对应路线，不能把有限测试当成渐近结论。

详见 [REFERENCE_ALIGNMENT.md](REFERENCE_ALIGNMENT.md)、[ROADMAP.md](../ROADMAP.md) 和 [PROGRESS.md](PROGRESS.md)。
