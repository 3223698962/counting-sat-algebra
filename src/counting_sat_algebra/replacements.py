"""
替换线性型与 H(T)（论文 §2.3）。

Replacement linear forms and the univariate H(T) from Section 2.3:

  ⟨P⃗_{e,v}(T), X⃗⟩ = ∑_i p_{e,v,i}(T) X_i
  with P⃗_{e,v}(0) = e_v  (standard basis), deg ≤ R

  H(T) = ∑_x ∏_e Φ_e( ⟨P⃗_{e,v}(T), x⟩ )_{v∈I_e}  ∈ ℚ[T]

  H(0) = #F,   deg H ≤ m · k · R   (k=子句宽度上界)

本原型在全部 2^n 赋值上用精确有理数求值 H(1),…,H(D+1)，
再插值恢复 H(0)。``eval_H_at`` 本身仍然枚举全部赋值。

Step 4 lives in ``hermite.py``. It reads several Hasse coefficients at each
sample point and recovers the same H(0) from fewer distinct nodes than the
D+1 value samples used here. Each node is still one pass over the cube.

Step 2 lives in ``blocked_sum.py``. For a hand-built family whose higher
T-coefficients lie in a fixed list of linear measurements, that module
regroups the same sum by blocks. The regrouping is checked against
``eval_H_at``. It does not replace this full-cube evaluator, and it does
not by itself yield a subexponential algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import List, Sequence, Tuple
import itertools
import random

from .cnf import CNF


CoeffPoly = List[Fraction]  # a0 + a1 T + … + aR T^R ; coeffs[j] = coeff of T^j


def eval_poly(coeffs: CoeffPoly, t: Fraction) -> Fraction:
    """Horner evaluation of univariate polynomial at t."""
    acc = Fraction(0)
    for c in reversed(coeffs):
        acc = acc * t + c
    return acc


def pad_poly(coeffs: CoeffPoly, deg: int) -> CoeffPoly:
    out = list(coeffs)
    while len(out) <= deg:
        out.append(Fraction(0))
    return out


@dataclass
class ReplacementFamily:
    """
    每个 (clause_index, local_occurrence_index) 一组系数行 P⃗(T)。
    rows[e][j] = 长度 n 的列表，每项是 deg≤R 的 CoeffPoly，且 P(0)=δ_v。
    """

    n: int
    R: int
    # rows[e][occ][i] = CoeffPoly for coefficient of X_i
    rows: List[List[List[CoeffPoly]]]
    # which variable each occurrence corresponds to (for Φ_e / literal polarity)
    occ_vars: List[List[int]]  # occ_vars[e][j] = v
    occ_neg: List[List[bool]]  # polarity of that occurrence

    def linear_form(self, e: int, occ: int, x: Sequence[int], t: Fraction) -> Fraction:
        """⟨P⃗_{e,occ}(t), x⟩"""
        s = Fraction(0)
        row = self.rows[e][occ]
        for i, coeffs in enumerate(row):
            if x[i]:
                s += eval_poly(coeffs, t)
        return s


def random_replacements(
    cnf: CNF,
    R: int = 2,
    rng: random.Random | None = None,
) -> ReplacementFamily:
    """
    构造随机低次替换：
      ⟨P,X⟩ = X_v + ∑_{j=1..R} T^j · (∑_i a_{j,i} X_i)
    其中 a_{j,i} 为小有理随机数。保证 P(0)=e_v。
    """
    rng = rng or random.Random()
    n = cnf.n
    rows: List[List[List[CoeffPoly]]] = []
    occ_vars: List[List[int]] = []
    occ_neg: List[List[bool]] = []

    def rand_frac() -> Fraction:
        # small rationals in {-2,-1,-1/2,0,1/2,1,2}
        choices = [
            Fraction(-2), Fraction(-1), Fraction(-1, 2),
            Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2),
        ]
        return rng.choice(choices)

    for cl in cnf.clauses:
        cl_rows: List[List[CoeffPoly]] = []
        vs: List[int] = []
        negs: List[bool] = []
        for v, neg in cl:
            # one CoeffPoly per variable i
            row: List[CoeffPoly] = []
            for i in range(n):
                coeffs = [Fraction(1 if i == v else 0)]  # constant term = δ_{v,i}
                for _j in range(1, R + 1):
                    coeffs.append(rand_frac())
                row.append(coeffs)
            cl_rows.append(row)
            vs.append(v)
            negs.append(neg)
        rows.append(cl_rows)
        occ_vars.append(vs)
        occ_neg.append(negs)

    return ReplacementFamily(n=n, R=R, rows=rows, occ_vars=occ_vars, occ_neg=occ_neg)


def phi_clause_numeric(
    values: Sequence[Fraction],
    negs: Sequence[bool],
) -> Fraction:
    """
    Φ_e on (possibly non-Boolean) inputs:
      L_v = val  or  1-val according to polarity
      Φ = 1 - ∏ (1 - L_v)
    Empty clause → 0.
    """
    if not values:
        return Fraction(0)
    prod = Fraction(1)
    for val, neg in zip(values, negs):
        L = (Fraction(1) - val) if neg else val
        prod *= Fraction(1) - L
    return Fraction(1) - prod


def eval_H_at(
    cnf: CNF,
    fam: ReplacementFamily,
    t: Fraction,
) -> Fraction:
    """
    穷举求和 H(t) = ∑_x ∏_e Φ_e( replacements at t )。
    仍是 2^n 规模——仅用于验证插值骨架，非论文的快速求值。
    """
    n = cnf.n
    total = Fraction(0)
    for bits in itertools.product((0, 1), repeat=n):
        prod = Fraction(1)
        for e, cl in enumerate(cnf.clauses):
            if not cl:
                prod = Fraction(0)
                break
            vals = [
                fam.linear_form(e, j, bits, t)
                for j in range(len(cl))
            ]
            prod *= phi_clause_numeric(vals, fam.occ_neg[e])
            if prod == 0:
                break
        total += prod
    return total


def degree_bound(cnf: CNF, R: int, k: int = 3) -> int:
    """deg H ≤ (#clauses) * k * R 的安全上界（含空子句也按 k 计）。"""
    # 每个 Φ_e 对替换输入次数 ≤ |Ie| ≤ k；每个输入次数 ≤ R
    return max(0, cnf.m * k * R)


def lagrange_eval(
    xs: Sequence[Fraction],
    ys: Sequence[Fraction],
    t: Fraction,
) -> Fraction:
    """
    Value at t of the unique polynomial of degree < len(xs) through (xs, ys).
    Nodes must be pairwise distinct. Exact rational arithmetic.
    """
    n = len(xs)
    if n != len(ys):
        raise ValueError("xs/ys length mismatch")
    if n == 0:
        raise ValueError("need at least one interpolation node")
    acc = Fraction(0)
    for i in range(n):
        num = Fraction(1)
        den = Fraction(1)
        xi = xs[i]
        for j in range(n):
            if i == j:
                continue
            num *= t - xs[j]
            den *= xi - xs[j]
        if den == 0:
            raise ValueError("interpolation nodes must be distinct")
        acc += ys[i] * (num / den)
    return acc


def lagrange_eval_at_zero(xs: Sequence[Fraction], ys: Sequence[Fraction]) -> Fraction:
    """
    唯一次数 < len(xs) 的插值多项式在 0 处的值（Lagrange）。
    要求 xs 互异。调用方使用 T=1..D+1，故 0 不是节点。
    """
    return lagrange_eval(xs, ys, Fraction(0))


def interpolate_H0(
    cnf: CNF,
    fam: ReplacementFamily | None = None,
    R: int = 2,
    rng: random.Random | None = None,
    k: int = 3,
) -> Tuple[Fraction, int, ReplacementFamily]:
    """
    在 T=1..D+1 求值 H，插值得 H(0)。
    返回 (H0, D, family)。
    """
    rng = rng or random.Random()
    fam = fam or random_replacements(cnf, R=R, rng=rng)
    D = degree_bound(cnf, R=fam.R, k=k)
    # 需要 D+1 个点确定 deg≤D 的多项式
    xs: List[Fraction] = []
    ys: List[Fraction] = []
    for t_int in range(1, D + 2):
        t = Fraction(t_int)
        xs.append(t)
        ys.append(eval_H_at(cnf, fam, t))
    h0 = lagrange_eval_at_zero(xs, ys)
    return h0, D, fam
