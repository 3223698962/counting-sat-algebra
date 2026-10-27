"""
k-CNF 表示与穷举计数 / k-CNF representation and brute-force counting.

对应论文 Section 2.2：
  Φ_e = 1 - ∏_v (1 - L_v),  L_v = X_v 或 1-X_v
  #F  = ∑_{x∈{0,1}^n} ∏_e Φ_e(x)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple
import itertools
import random

# 一个文字：(变量下标 0..n-1, 是否取反)
# A literal: (var_index 0..n-1, is_negated)
Literal = Tuple[int, bool]
Clause = List[Literal]


@dataclass
class CNF:
    """k-CNF：clauses 为子句列表；n 为声明变量个数。"""

    n: int
    clauses: List[Clause]

    def __post_init__(self) -> None:
        if self.n < 0:
            raise ValueError("n must be nonnegative")
        for e, cl in enumerate(self.clauses):
            for v, _neg in cl:
                if not (0 <= v < self.n):
                    raise ValueError(f"clause {e}: var {v} out of range [0,{self.n})")

    @property
    def m(self) -> int:
        return len(self.clauses)

    def copy(self) -> "CNF":
        return CNF(self.n, [list(cl) for cl in self.clauses])


def literal_value(x: Sequence[int], lit: Literal) -> int:
    """布尔赋值下文字取值：正文字为 x_v，负文字为 1-x_v。"""
    v, neg = lit
    return (1 - x[v]) if neg else x[v]


def clause_poly_bool(x: Sequence[int], clause: Clause) -> int:
    """
    Φ_e(x) on {0,1}: 1 if clause satisfied else 0.
    Φ_e = 1 - ∏_v (1 - L_v).
    Empty clause → always unsatisfied → 0.
    """
    if not clause:
        return 0
    prod = 1
    for lit in clause:
        prod *= 1 - literal_value(x, lit)
        if prod == 0:
            break
    return 1 - prod


def is_satisfied(x: Sequence[int], cnf: CNF) -> bool:
    return all(clause_poly_bool(x, cl) == 1 for cl in cnf.clauses)


def brute_count(cnf: CNF) -> int:
    """穷举 #F = |{x∈{0,1}^n : F(x)=1}|。n 宜 ≤ 20。"""
    n = cnf.n
    if n > 20:
        raise ValueError(f"brute_count refuses n={n} > 20")
    total = 0
    for bits in itertools.product((0, 1), repeat=n):
        if is_satisfied(bits, cnf):
            total += 1
    return total


def random_3cnf(
    n: int,
    m: int,
    rng: random.Random | None = None,
    allow_tautology: bool = False,
) -> CNF:
    """
    均匀随机 3-CNF：每个子句 3 个不同变量，随机极性。
    若 allow_tautology=False，跳过含 x∨¬x 的退化子句（本生成方式不会产生）。
    """
    if n < 3:
        raise ValueError("need n >= 3 for 3-CNF")
    rng = rng or random.Random()
    clauses: List[Clause] = []
    for _ in range(m):
        vars_ = rng.sample(range(n), 3)
        cl = [(v, bool(rng.getrandbits(1))) for v in vars_]
        if not allow_tautology:
            # 同变量正负同时出现不可能（变量互异）
            pass
        clauses.append(cl)
    return CNF(n, clauses)


def shift_cnf(cnf: CNF, offset: int) -> CNF:
    """将所有变量下标加上 offset（用于拼接独立分量）。"""
    new_clauses = [
        [(v + offset, neg) for v, neg in cl] for cl in cnf.clauses
    ]
    return CNF(cnf.n + offset, new_clauses)  # temporary; caller sets n


def concatenate_disjoint(f1: CNF, f2: CNF) -> CNF:
    """
    拼接两个变量不相交的公式：f2 的变量接到 f1 之后。
    对应论文 §2.2 动机：无跨组子句时 #F = ∏_g N_g。
    """
    shifted = [
        [(v + f1.n, neg) for v, neg in cl] for cl in f2.clauses
    ]
    return CNF(f1.n + f2.n, [list(cl) for cl in f1.clauses] + shifted)



def restrict_clean(cnf: CNF, fixed: dict[int, int]) -> CNF:
    """同上，但更清晰地按 fixed 化简。"""
    new_clauses: List[Clause] = []
    for cl in cnf.clauses:
        simplified: Clause = []
        satisfied = False
        for v, neg in cl:
            if v not in fixed:
                simplified.append((v, neg))
                continue
            bit = fixed[v]
            lit_true = (bit == 0) if neg else (bit == 1)
            if lit_true:
                satisfied = True
                break
            # lit false → drop
        if not satisfied:
            new_clauses.append(simplified)
    return CNF(cnf.n, new_clauses)
