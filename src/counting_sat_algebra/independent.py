"""
方案 A：独立分量乘积恒等式（论文 §2.2 动机情形）。

若变量分成若干组，且每个子句的变量都落在同一组内，则
  #F = ∏_g N_g
其中 N_g 是第 g 组上局部公式的满足赋值数。
"""

from __future__ import annotations

from typing import List, Tuple
import random

from .cnf import CNF, brute_count, concatenate_disjoint, random_3cnf


def product_identity_check(
    n1: int,
    m1: int,
    n2: int,
    m2: int,
    rng: random.Random | None = None,
) -> Tuple[int, int, int, int, bool]:
    """
    构造两个不相交随机 3-CNF，拼接后验证 #F = #F1 * #F2。
    返回 (n_total, m_total, count_combined, product, ok)。
    """
    rng = rng or random.Random()
    f1 = random_3cnf(n1, m1, rng=rng)
    f2 = random_3cnf(n2, m2, rng=rng)
    c1 = brute_count(f1)
    c2 = brute_count(f2)
    comb = concatenate_disjoint(f1, f2)
    c_comb = brute_count(comb)
    product = c1 * c2
    return comb.n, comb.m, c_comb, product, (c_comb == product)


def run_independent_suite(
    trials: int = 20,
    seed: int = 0,
) -> List[dict]:
    rng = random.Random(seed)
    rows: List[dict] = []
    for i in range(trials):
        n1 = rng.randint(2, 5)
        n2 = rng.randint(2, 5)
        # 保证各侧 n>=3 以便 3-CNF；若较小则用 2-CNF 风格：放宽
        n1 = max(n1, 3)
        n2 = max(n2, 3)
        m1 = rng.randint(n1, 2 * n1)
        m2 = rng.randint(n2, 2 * n2)
        n, m, c_comb, product, ok = product_identity_check(n1, m1, n2, m2, rng=rng)
        rows.append(
            {
                "trial": i,
                "n": n,
                "m": m,
                "n1": n1,
                "n2": n2,
                "count": c_comb,
                "product": product,
                "ok": ok,
            }
        )
    return rows
