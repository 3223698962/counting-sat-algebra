"""
自归约找满足赋值：逐位固定，用计数器检查剩余 #F>0（论文引言）。
"""

from __future__ import annotations

from typing import List, Optional

from .cnf import CNF, brute_count, is_satisfied, restrict_clean


def find_assignment_by_self_reduction(cnf: CNF) -> Optional[List[int]]:
    """
    若 #F>0，返回一个满足赋值；否则 None。
    对每位尝试 0，若 restrict 后计数 >0 则取 0，否则取 1。
    """
    fixed: dict[int, int] = {}
    for v in range(cnf.n):
        # try bit 0
        trial0 = dict(fixed)
        trial0[v] = 0
        c0 = brute_count(restrict_clean(cnf, trial0))
        if c0 > 0:
            fixed[v] = 0
            continue
        trial1 = dict(fixed)
        trial1[v] = 1
        c1 = brute_count(restrict_clean(cnf, trial1))
        if c1 > 0:
            fixed[v] = 1
            continue
        return None
    assignment = [fixed[v] for v in range(cnf.n)]
    if not is_satisfied(assignment, cnf):
        # 理论上不应发生
        return None
    return assignment
