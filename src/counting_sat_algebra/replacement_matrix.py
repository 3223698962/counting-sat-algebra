"""
Small-parameter slice of the low-degree replacement matrix.

The object built here is the auxiliary matrix

    M0 = V0 f_I(W0)

and the coefficient rows that come from it. Coordinates are pairs of length-Q
schedules with entries in {0, ..., r-1}. There are D = r**(2Q) of them. V0
exchanges the two schedules in every coordinate, W0 = (1/2) Σ_s (I - Π_s) where
Π_s exchanges schedule entry s, and f_I is the degree-at-most-(R-1) polynomial
with f_I(j) = (-1)^j on R consecutive integers.

From a row u and a column y with M0[u, y] ≠ 0, write E(u, y) for the positions
s ∈ {1, ..., Q} at which the column's first schedule disagrees with the row's
second schedule. The auxiliary replacement is

    ⟨P_u(T), ζ⟩ = ((I - M0) ζ)_u
                  + Σ_{y : M0[u, y] ≠ 0} M0[u, y] φ_{u,y}(T) ζ_y,

    φ_{u,y}(T) = Π_{s ∈ E(u, y)} (1 - T/s).

At T = 0 every factor φ equals 1, so the row restores the coordinate ζ_u. Each
factor has degree |E(u, y)|, and the support rule on M0 forces that degree to
be at most R - 1.

What this module does not build
-------------------------------
The parameter list associated with a full matrix asks for R ≥ 32 and 6R ≤ Q
before its tail bound is applied. The smallest integers are R = 32 and
Q = 192. At r = 2 the side length is then 2**384. ``SLICE_MAX_DIMENSION``
refuses anything above 256, and Q > 30 is refused before a matrix is stored,
so that matrix is not allocated, not factored, and not run. The Step 3
schedules below are written out by hand. Step 7 checks two further explicit
placements, at side lengths 81 and 256, and a balance-check scaffold that
refuses to choose schedules when 2(b+w)D**2 ≤ M and 2Qr**2 ≤ n fail. Step 8
runs a successive-choice search when those inequalities hold. On the admitted
rows the legacy Boolean mode refuses before design_problems_at because that walk would exceed a
Boolean-point budget. The inequalities need n on the order of D**2.
Sparsification and a subexponential count are not implemented. Hermite
interpolation of H is a separate exact check in ``hermite.py``; it is not
part of this matrix, and each of its sample points still enumerates the cube.
Theorem 1.1 is not implemented. Nothing here runs in 2^{o(n)}.

The same formulas are checked at four tiny parameter points, including the
four-coordinate identity

    ⟨P_{(00,11)}(T), ζ⟩ = ζ_(00,11) - (T/2) ζ_(10,01) - T ζ_(01,10)

at r = 2, Q = 2, R = 2. That point has R = Q, which sits outside the
hypothesis R < Q used for a runtime bound. The identity itself is exact and is
checked directly.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .bounded_occurrence import equispaced_degree
from .cnf import CNF, brute_count
from .replacements import (
    ReplacementFamily,
    eval_H_at,
    eval_poly,
    lagrange_eval,
    lagrange_eval_at_zero,
)

# Side length r**(2Q). Above this the slice refuses to allocate M0.
# 256 admits the Step 7 point r=2, Q=4. The next pure powers, 729 (r=3, Q=3)
# and 1024 (r=2, Q=5), stay above the cap. Q > 30 is refused before allocation
# even if a caller passes a larger max_dimension, so full M0 is not stored.
SLICE_MAX_DIMENSION = 256

# Smallest integers with R ≥ 32 and 6R ≤ Q. At r = 2 the side length is 2**384.
FULL_M0_MIN_R = 32
FULL_M0_MIN_Q = 6 * FULL_M0_MIN_R
FULL_M0_LOG2_DIMENSION_AT_R2 = 2 * FULL_M0_MIN_Q


def _as_index(value: object, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{what} must be an int, got {value!r}")
    return value


def _trim(coeffs: Sequence[Fraction]) -> list[Fraction]:
    out = list(coeffs)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out or [Fraction(0)]


def _add_poly(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    size = max(len(left), len(right))
    out = [Fraction(0) for _ in range(size)]
    for i, coeff in enumerate(left):
        out[i] += coeff
    for i, coeff in enumerate(right):
        out[i] += coeff
    return _trim(out)


def _scale_poly(coeffs: Sequence[Fraction], factor: Fraction) -> list[Fraction]:
    return _trim([coeff * factor for coeff in coeffs])


def _mul_linear(coeffs: Sequence[Fraction], constant: Fraction, slope: Fraction) -> list[Fraction]:
    """Multiply by ``constant + slope T``."""
    out = [Fraction(0) for _ in range(len(coeffs) + 1)]
    for i, coeff in enumerate(coeffs):
        out[i] += coeff * constant
        out[i + 1] += coeff * slope
    return _trim(out)


def _factorial(n: int) -> int:
    if n < 0:
        raise ValueError(f"factorial of {n}")
    return math.factorial(n)


def vector_in_row_span(
    rows: Sequence[Sequence[Fraction]],
    target: Sequence[Fraction],
) -> bool:
    """Exact check: ``target`` is a rational linear combination of ``rows``.

    Empty ``target`` is in every span. A zero target is in the span of no rows.
    """
    width = len(target)
    if any(len(row) != width for row in rows):
        raise ValueError("row widths must match the target")
    reduced = [list(row) for row in rows]
    left = list(target)
    pivot_row = 0
    for col in range(width):
        pivot = None
        for i in range(pivot_row, len(reduced)):
            if reduced[i][col] != 0:
                pivot = i
                break
        if pivot is None:
            if left[col] != 0:
                return False
            continue
        reduced[pivot_row], reduced[pivot] = reduced[pivot], reduced[pivot_row]
        pivot_value = reduced[pivot_row][col]
        factor = left[col] / pivot_value
        if factor != 0:
            for j in range(col, width):
                left[j] -= factor * reduced[pivot_row][j]
        for i in range(len(reduced)):
            if i == pivot_row:
                continue
            row_factor = reduced[i][col] / pivot_value
            if row_factor == 0:
                continue
            for j in range(col, width):
                reduced[i][j] -= row_factor * reduced[pivot_row][j]
        pivot_row += 1
    return all(entry == 0 for entry in left)


def default_interval_start(r: int, Q: int, R: int) -> int:
    """Floor of Q(r-1)/(2r) - (R-1)/2.

    This is the start used for the asymptotic interval. At small parameters it
    can be negative, in which case ``build_auxiliary_slice`` rejects it unless
    the caller passes a start that does lie in range.
    """
    r = _as_index(r, "r")
    Q = _as_index(Q, "Q")
    R = _as_index(R, "R")
    if r < 1 or Q < 1 or R < 1:
        raise ValueError("r, Q, and R must be positive")
    center = Fraction(Q * (r - 1), 2 * r)
    return math.floor(center - Fraction(R - 1, 2))


def scalar_f_coefficients(A0: int, R: int) -> tuple[Fraction, ...]:
    """Coefficients of f_I, constant term first.

    f_I(z) = (-1)^{A0} Σ_{j=0}^{R-1} (-2)^j binom(z - A0, j). On the integers
    A0, ..., A0+R-1 this equals (-1)^j.
    """
    A0 = _as_index(A0, "A0")
    R = _as_index(R, "R")
    if R < 1:
        raise ValueError(f"R must be positive, got {R}")
    sign = Fraction(-1) ** A0
    acc = [Fraction(0)]
    for j in range(R):
        poly = [Fraction(1)]
        for a in range(j):
            poly = _mul_linear(poly, Fraction(-(A0 + a)), Fraction(1))
        scale = sign * (Fraction(-2) ** j) / Fraction(_factorial(j))
        acc = _add_poly(acc, _scale_poly(poly, scale))
    return tuple(_trim(acc))


def _sequences(r: int, Q: int) -> list[tuple[int, ...]]:
    return list(itertools.product(range(r), repeat=Q))


def _labels(r: int, Q: int) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    seqs = _sequences(r, Q)
    return [(a, d) for a in seqs for d in seqs]


def _swap_position(
    label: tuple[tuple[int, ...], tuple[int, ...]],
    position: int,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    first = list(label[0])
    second = list(label[1])
    first[position], second[position] = second[position], first[position]
    return (tuple(first), tuple(second))


def _full_swap(
    label: tuple[tuple[int, ...], tuple[int, ...]],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    return (label[1], label[0])


def _local_eigenvectors(r: int) -> tuple[list[tuple[Fraction, ...]], list[tuple[Fraction, ...]]]:
    """Eigenvectors of one pair-swap. Plus first, then minus.

    Plus: e_{αα}, and e_{αδ}+e_{δα} for α < δ. Minus: e_{αδ}-e_{δα} for α < δ.
    Entries lie in {0, 1, -1}.
    """
    plus: list[tuple[Fraction, ...]] = []
    minus: list[tuple[Fraction, ...]] = []

    def blank() -> list[Fraction]:
        return [Fraction(0) for _ in range(r * r)]

    for alpha in range(r):
        vector = blank()
        vector[alpha * r + alpha] = Fraction(1)
        plus.append(tuple(vector))
    for alpha in range(r):
        for delta in range(alpha + 1, r):
            symmetric = blank()
            symmetric[alpha * r + delta] = Fraction(1)
            symmetric[delta * r + alpha] = Fraction(1)
            plus.append(tuple(symmetric))
            skew = blank()
            skew[alpha * r + delta] = Fraction(1)
            skew[delta * r + alpha] = Fraction(-1)
            minus.append(tuple(skew))
    return plus, minus


def _product_vector(
    choices: Sequence[tuple[Fraction, ...]],
    labels: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
    r: int,
) -> tuple[Fraction, ...]:
    out: list[Fraction] = []
    for first, second in labels:
        value = Fraction(1)
        for position, local in enumerate(choices):
            value *= local[first[position] * r + second[position]]
            if value == 0:
                break
        out.append(value)
    return tuple(out)


def _eigenvectors(
    r: int,
    Q: int,
    labels: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
) -> list[tuple[int, tuple[Fraction, ...]]]:
    plus, minus = _local_eigenvectors(r)
    catalog = [(vector, False) for vector in plus] + [(vector, True) for vector in minus]
    vectors: list[tuple[int, tuple[Fraction, ...]]] = []
    for choice in itertools.product(range(len(catalog)), repeat=Q):
        selected = [catalog[index][0] for index in choice]
        negatives = sum(1 for index in choice if catalog[index][1])
        vectors.append((negatives, _product_vector(selected, labels, r)))
    return vectors


def binomial_exceptional_count(r: int, Q: int, interval: Sequence[int]) -> int:
    """How many product basis vectors have their minus-count outside ``interval``.

    A position has r(r+1)/2 plus eigenvectors and r(r-1)/2 minus eigenvectors.
    The count is exact; it is not a tail estimate.
    """
    inside = set(interval)
    n_plus = r * (r + 1) // 2
    n_minus = r * (r - 1) // 2
    total = 0
    for negatives in range(Q + 1):
        if negatives in inside:
            continue
        total += math.comb(Q, negatives) * (n_minus**negatives) * (n_plus ** (Q - negatives))
    return total


def _mismatch_positions(
    row_label: tuple[tuple[int, ...], tuple[int, ...]],
    column_label: tuple[tuple[int, ...], tuple[int, ...]],
) -> tuple[int, ...]:
    """1-based positions s where the column's first schedule disagrees with the row's second."""
    _row_first, row_second = row_label
    column_first, _column_second = column_label
    return tuple(
        position + 1
        for position, (left, right) in enumerate(zip(column_first, row_second))
        if left != right
    )


def _phi(positions: Sequence[int]) -> list[Fraction]:
    poly = [Fraction(1)]
    for position in positions:
        poly = _mul_linear(poly, Fraction(1), Fraction(-1, position))
    return _trim(poly)


def _zero_matrix(size: int) -> list[list[Fraction]]:
    return [[Fraction(0) for _ in range(size)] for _ in range(size)]


def _horner_matrix(coeffs: Sequence[Fraction], matrix: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    acc = _zero_matrix(size)
    for coeff in reversed(coeffs):
        acc = _matmul(matrix, acc)
        for i in range(size):
            acc[i][i] += coeff
    return acc


def _matmul(
    left: Sequence[Sequence[Fraction]],
    right: Sequence[Sequence[Fraction]],
) -> list[list[Fraction]]:
    size = len(left)
    width = len(right[0]) if right else 0
    mid = len(right)
    out = [[Fraction(0) for _ in range(width)] for _ in range(size)]
    for i in range(size):
        left_row = left[i]
        out_row = out[i]
        for k in range(mid):
            factor = left_row[k]
            if factor == 0:
                continue
            right_row = right[k]
            for j in range(width):
                entry = right_row[j]
                if entry != 0:
                    out_row[j] += factor * entry
    return out


def _matvec(matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(
        sum((entry * value for entry, value in zip(row, vector)), start=Fraction(0))
        for row in matrix
    )


@dataclass(frozen=True)
class AuxiliarySlice:
    """One tiny M0 together with the auxiliary rows P_u(T).

    ``rows[u][y]`` is the coefficient polynomial of coordinate y in row u,
    constant term first. ``C_rows`` are the product basis vectors whose
    minus-count lies outside ``interval``; their entries are in {0, 1, -1}.
    """

    r: int
    Q: int
    R: int
    A0: int
    interval: tuple[int, ...]
    dimension: int
    labels: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]
    M0: tuple[tuple[Fraction, ...], ...]
    C_rows: tuple[tuple[Fraction, ...], ...]
    rows: tuple[tuple[tuple[Fraction, ...], ...], ...]
    f_coefficients: tuple[Fraction, ...]
    used_default_interval: bool

    def __repr__(self) -> str:
        return (
            f"AuxiliarySlice(r={self.r}, Q={self.Q}, R={self.R}, A0={self.A0}, "
            f"D={self.dimension}, measurements={len(self.C_rows)})"
        )

    def index_of(self, first: Sequence[int], second: Sequence[int]) -> int:
        key = (tuple(first), tuple(second))
        for index, label in enumerate(self.labels):
            if label == key:
                return index
        raise KeyError(key)


def _dimension_refusal(r: int, Q: int, dimension: int, cap: int) -> ValueError:
    return ValueError(
        f"auxiliary side length D=r**(2Q)={dimension} for r={r}, Q={Q} exceeds the "
        f"slice cap {cap}. The full matrix is not built: the parameter list with "
        f"R≥{FULL_M0_MIN_R} and 6R≤Q forces Q≥{FULL_M0_MIN_Q}, and at r=2 the side "
        f"length is 2**{FULL_M0_LOG2_DIMENSION_AT_R2}"
    )


def build_auxiliary_slice(
    r: int,
    Q: int,
    R: int,
    A0: int | None = None,
    max_dimension: int = SLICE_MAX_DIMENSION,
) -> AuxiliarySlice:
    """Build M0 and the auxiliary rows, or refuse if the side length is above the cap.

    ``A0`` defaults to ``default_interval_start``. The interval is
    {A0, ..., A0+R-1} and must sit inside {0, ..., Q}.
    """
    r = _as_index(r, "r")
    Q = _as_index(Q, "Q")
    R = _as_index(R, "R")
    max_dimension = _as_index(max_dimension, "max_dimension")
    if r < 2:
        raise ValueError(f"r must be an integer ≥ 2, got {r}")
    if Q < 1:
        raise ValueError(f"Q must be a positive integer, got {Q}")
    if R < 1:
        raise ValueError(f"R must be a positive integer, got {R}")
    if max_dimension < 1:
        raise ValueError(f"max_dimension must be positive, got {max_dimension}")
    if Q > 30:
        raise _dimension_refusal(r, Q, r ** (2 * Q), max_dimension)
    dimension = r ** (2 * Q)
    if dimension > max_dimension:
        raise _dimension_refusal(r, Q, dimension, max_dimension)
    if R > Q + 1:
        raise ValueError(
            f"R={R} does not fit: an interval of {R} integers cannot sit inside {{0, ..., {Q}}}"
        )
    used_default = A0 is None
    if A0 is None:
        A0 = default_interval_start(r, Q, R)
    else:
        A0 = _as_index(A0, "A0")
    if A0 < 0 or A0 + R - 1 > Q:
        raise ValueError(
            f"interval {{{A0}, ..., {A0 + R - 1}}} is not contained in {{0, ..., {Q}}}. "
            "The default start can fall outside that range at small parameters; "
            "pass an explicit A0 that lies inside it, or choose a larger Q"
        )
    interval = tuple(range(A0, A0 + R))
    labels = tuple(_labels(r, Q))
    if len(labels) != dimension:
        raise RuntimeError("label count drifted from r**(2Q)")
    index = {label: i for i, label in enumerate(labels)}

    size = dimension
    W0 = _zero_matrix(size)
    for position in range(Q):
        for label_index, label in enumerate(labels):
            image = index[_swap_position(label, position)]
            if image == label_index:
                continue
            W0[label_index][label_index] += Fraction(1, 2)
            W0[label_index][image] -= Fraction(1, 2)

    f_coefficients = scalar_f_coefficients(A0, R)
    f_of_W = _horner_matrix(f_coefficients, W0)
    M0 = _zero_matrix(size)
    for label_index, label in enumerate(labels):
        M0[label_index] = list(f_of_W[index[_full_swap(label)]])

    inside = set(interval)
    C_rows = tuple(
        vector for negatives, vector in _eigenvectors(r, Q, labels) if negatives not in inside
    )

    polynomial_rows: list[tuple[tuple[Fraction, ...], ...]] = []
    for u, row_label in enumerate(labels):
        coordinate_polys: list[list[Fraction]] = []
        for y in range(size):
            diagonal = Fraction(1) if u == y else Fraction(0)
            coordinate_polys.append([diagonal - M0[u][y]])
        for y, column_label in enumerate(labels):
            entry = M0[u][y]
            if entry == 0:
                continue
            factor = _phi(_mismatch_positions(row_label, column_label))
            coordinate_polys[y] = _add_poly(coordinate_polys[y], _scale_poly(factor, entry))
        polynomial_rows.append(tuple(tuple(poly) for poly in coordinate_polys))

    return AuxiliarySlice(
        r=r,
        Q=Q,
        R=R,
        A0=A0,
        interval=interval,
        dimension=dimension,
        labels=labels,
        M0=tuple(tuple(row) for row in M0),
        C_rows=C_rows,
        rows=tuple(polynomial_rows),
        f_coefficients=f_coefficients,
        used_default_interval=used_default,
    )


def auxiliary_problems(sl: AuxiliarySlice) -> list[str]:
    """Restoration, degree, support, kernel, and denominator checks for one slice.

    An empty list means every check below passed. These are exact rational
    identities on this finite matrix, not a runtime bound.
    """
    problems: list[str] = []
    interval = set(sl.interval)
    for j in sl.interval:
        expected = Fraction(-1) ** j
        if eval_poly(list(sl.f_coefficients), Fraction(j)) != expected:
            problems.append(f"f({j}) is not {expected}")
    if len(sl.M0) != sl.dimension or len(sl.rows) != sl.dimension:
        problems.append("matrix shape drifted")
        return problems

    matrix_denominator = _factorial(sl.R - 1)
    row_denominator = matrix_denominator * _factorial(sl.Q)
    for u in range(sl.dimension):
        if len(sl.M0[u]) != sl.dimension:
            problems.append(f"M0 row {u} has the wrong length")
            return problems
        for entry in sl.M0[u]:
            if (Fraction(matrix_denominator) * entry).denominator != 1:
                problems.append("(R-1)! M0 is not an integer matrix")
                return problems
        for y in range(sl.dimension):
            poly = sl.rows[u][y]
            degree = len(poly) - 1
            if degree > sl.R - 1:
                problems.append(
                    f"row {u} coordinate {y} has degree {degree}, above R-1={sl.R - 1}"
                )
                return problems
            if eval_poly(list(poly), Fraction(0)) != (Fraction(1) if u == y else Fraction(0)):
                problems.append(f"T=0 failed to restore coordinate {u}")
                return problems
            for coeff in poly:
                if (Fraction(row_denominator) * coeff).denominator != 1:
                    problems.append(
                        "(R-1)! Q! does not clear the denominators of the auxiliary rows"
                    )
                    return problems
            if sl.M0[u][y] == 0:
                continue
            mismatches = len(_mismatch_positions(sl.labels[u], sl.labels[y]))
            if mismatches > sl.R - 1:
                problems.append(
                    f"M0[{u}, {y}] is nonzero with {mismatches} mismatches, above R-1"
                )
                return problems

    expected_exceptional = binomial_exceptional_count(sl.r, sl.Q, sl.interval)
    if len(sl.C_rows) != expected_exceptional:
        problems.append(
            f"measurement row count {len(sl.C_rows)} != binomial count {expected_exceptional}"
        )

    fixed = 0
    for negatives, vector in _eigenvectors(sl.r, sl.Q, sl.labels):
        if any(entry not in (Fraction(-1), Fraction(0), Fraction(1)) for entry in vector):
            problems.append("an eigenbasis vector has an entry outside {-1, 0, 1}")
            return problems
        if negatives not in interval:
            continue
        fixed += 1
        if _matvec(sl.M0, vector) != vector:
            problems.append(f"M0 does not fix a basis vector with minus-count {negatives}")
            return problems
        for row_index, row in enumerate(sl.C_rows):
            dot = sum((a * b for a, b in zip(row, vector)), start=Fraction(0))
            if dot != 0:
                problems.append(
                    f"measurement row {row_index} is not orthogonal to a fixed basis vector"
                )
                return problems
    if fixed + len(sl.C_rows) != sl.dimension:
        problems.append("fixed basis vectors plus measurement rows do not cover the dimension")

    for u in range(sl.dimension):
        correction = [Fraction(1 if u == y else 0) - sl.M0[u][y] for y in range(sl.dimension)]
        if not vector_in_row_span(sl.C_rows, correction):
            problems.append(f"row {u} of I-M0 is not in the span of the measurement rows")
            return problems
    return problems


@dataclass
class PlacedFormula:
    """Auxiliary rows pulled back onto one explicit CNF.

    ``variable_schedules[v]`` and ``clause_schedules[e]`` are length-Q tuples
    with entries in {0, ..., r-1}. They are data, not the output of a balancer.
    ``measurement_rows`` are the auxiliary measurement rows after filling
    coordinates from the original variables, one block of rows per copy.
    ``family.R`` is R-1, the degree bound on those rows.
    """

    cnf: CNF
    auxiliary: AuxiliarySlice
    variable_schedules: tuple[tuple[int, ...], ...]
    clause_schedules: tuple[tuple[int, ...], ...]
    family: ReplacementFamily
    measurement_rows: tuple[tuple[Fraction, ...], ...]
    copies: int
    occurrence_mass: int
    max_width: int

    def __repr__(self) -> str:
        return (
            f"PlacedFormula(n={self.cnf.n}, m={self.cnf.m}, copies={self.copies}, "
            f"{self.auxiliary!r})"
        )

    @property
    def cap_occurrences(self) -> int:
        """M · (R - 1), with M the number of literal occurrences."""
        return self.occurrence_mass * (self.auxiliary.R - 1)

    @property
    def cap_step1(self) -> int:
        """m · k · R, the Step 1 shape of cap, using the matrix parameter R and max width k."""
        return self.cnf.m * self.max_width * self.auxiliary.R


def place_formula(
    cnf: CNF,
    auxiliary: AuxiliarySlice,
    variable_schedules: Sequence[Sequence[int]],
    clause_schedules: Sequence[Sequence[int]],
) -> PlacedFormula:
    """Pull auxiliary rows back along an explicit placement of variables and clauses.

    Occurrence (e, v) is written into the coordinate whose schedules are
    (schedule of v, schedule of e). Repeated coordinates use successive copies.
    Unoccupied auxiliary coordinates contribute 0.
    """
    if len(variable_schedules) != cnf.n:
        raise ValueError(
            f"expected {cnf.n} variable schedules, got {len(variable_schedules)}"
        )
    if len(clause_schedules) != cnf.m:
        raise ValueError(
            f"expected {cnf.m} clause schedules, got {len(clause_schedules)}"
        )
    variables = tuple(_normalize_schedule(schedule, auxiliary, "variable") for schedule in variable_schedules)
    clauses = tuple(_normalize_schedule(schedule, auxiliary, "clause") for schedule in clause_schedules)

    occurrences: list[tuple[int, int, int, tuple[tuple[int, ...], tuple[int, ...]]]] = []
    for clause_index, clause in enumerate(cnf.clauses):
        for occ, (var, _neg) in enumerate(clause):
            label = (variables[var], clauses[clause_index])
            occurrences.append((clause_index, occ, var, label))

    counts: dict[tuple[tuple[int, ...], tuple[int, ...]], int] = {}
    assigned_copies: list[int] = []
    for _clause_index, _occ, _var, label in occurrences:
        copy = counts.get(label, 0)
        assigned_copies.append(copy)
        counts[label] = copy + 1
    copies = max(counts.values(), default=0)

    label_index = {label: i for i, label in enumerate(auxiliary.labels)}
    # placements[copy][coordinate] = variable index
    placements: list[dict[int, int]] = [{} for _ in range(copies)]
    for (clause_index, occ, var, label), copy in zip(occurrences, assigned_copies):
        coordinate = label_index[label]
        previous = placements[copy].get(coordinate)
        if previous is not None and previous != var:
            raise RuntimeError(
                f"copy {copy} coordinate {coordinate} was assigned two variables"
            )
        placements[copy][coordinate] = var

    measurement_rows: list[tuple[Fraction, ...]] = []
    for placement in placements:
        for measurement in auxiliary.C_rows:
            acc = [Fraction(0) for _ in range(cnf.n)]
            for coordinate, coeff in enumerate(measurement):
                if coeff == 0:
                    continue
                var = placement.get(coordinate)
                if var is None:
                    continue
                acc[var] += coeff
            measurement_rows.append(tuple(acc))

    rows: list[list[list[list[Fraction]]]] = []
    occ_vars: list[list[int]] = []
    occ_neg: list[list[bool]] = []
    cursor = 0
    for clause_index, clause in enumerate(cnf.clauses):
        clause_rows: list[list[list[Fraction]]] = []
        variables_in_clause: list[int] = []
        negated: list[bool] = []
        for occ, (var, neg) in enumerate(clause):
            _c, _o, _v, label = occurrences[cursor]
            copy = assigned_copies[cursor]
            cursor += 1
            placement = placements[copy]
            coordinate_polys = auxiliary.rows[label_index[label]]
            pulled: list[list[Fraction]] = [[Fraction(0)] for _ in range(cnf.n)]
            for coordinate, poly in enumerate(coordinate_polys):
                owner = placement.get(coordinate)
                if owner is None:
                    continue
                pulled[owner] = _add_poly(pulled[owner], list(poly))
            clause_rows.append(pulled)
            variables_in_clause.append(var)
            negated.append(neg)
        rows.append(clause_rows)
        occ_vars.append(variables_in_clause)
        occ_neg.append(negated)

    family = ReplacementFamily(
        n=cnf.n,
        R=auxiliary.R - 1,
        rows=rows,
        occ_vars=occ_vars,
        occ_neg=occ_neg,
    )
    mass = sum(len(clause) for clause in cnf.clauses)
    width = max((len(clause) for clause in cnf.clauses), default=0)
    return PlacedFormula(
        cnf=cnf,
        auxiliary=auxiliary,
        variable_schedules=variables,
        clause_schedules=clauses,
        family=family,
        measurement_rows=tuple(measurement_rows),
        copies=copies,
        occurrence_mass=mass,
        max_width=width,
    )


def _normalize_schedule(
    schedule: Sequence[int],
    auxiliary: AuxiliarySlice,
    what: str,
) -> tuple[int, ...]:
    if len(schedule) != auxiliary.Q:
        raise ValueError(f"{what} schedule has length {len(schedule)}, expected Q={auxiliary.Q}")
    out: list[int] = []
    for position, entry in enumerate(schedule):
        value = _as_index(entry, f"{what} schedule[{position}]")
        if not (0 <= value < auxiliary.r):
            raise ValueError(
                f"{what} schedule[{position}]={value} is outside {{0, ..., {auxiliary.r - 1}}}"
            )
        out.append(value)
    return tuple(out)


def groups_at(placed: PlacedFormula, t: int) -> tuple[tuple[int, ...], ...]:
    """Variable groups at a 1-based sample point t ∈ {1, ..., Q}."""
    _check_t(placed, t)
    buckets: list[list[int]] = [[] for _ in range(placed.auxiliary.r)]
    position = t - 1
    for var, schedule in enumerate(placed.variable_schedules):
        buckets[schedule[position]].append(var)
    return tuple(tuple(bucket) for bucket in buckets)


def _check_t(placed: PlacedFormula, t: int) -> None:
    t = _as_index(t, "t")
    if not (1 <= t <= placed.auxiliary.Q):
        raise ValueError(f"t must lie in {{1, ..., {placed.auxiliary.Q}}}, got {t}")


def design_problems_at(placed: PlacedFormula, t: int, *, boolean_check: bool = False) -> list[str]:
    """Legal-design failures at one sample point.

    Checked here, for this explicit placement:

    - the schedules partition the variables;
    - each group has size at most 2n/r, read as ``size * r ≤ 2n``;
    - each pulled-back measurement row is integral and has absolute sum at most
      the auxiliary side length D;
    - for every occurrence, at this t, the coefficient vector outside the
      clause's assigned group lies in the span of the measurement rows
      (the rational kernel condition implies the statement on Boolean points).

    The default eliminates the restricted measurement matrix once per group.
    ``boolean_check=True`` additionally runs the original exhaustive Boolean
    comparison, for small-instance validation. Callers must budget that walk.

    For r = 2 the size bound is ``size ≤ n``, so it holds for every partition.
    The kernel condition is the part that can fail. An empty list means this
    point passed. It does not say the measurement set is small enough to beat
    2^n, and it does not choose the schedules.
    """
    _check_t(placed, t)
    problems: list[str] = []
    groups = groups_at(placed, t)
    n = placed.cnf.n
    r = placed.auxiliary.r
    seen = [-1] * n
    for group, variables in enumerate(groups):
        if len(variables) * r > 2 * n:
            problems.append(
                f"t={t}: group {group} has size {len(variables)}, and size * r > 2n"
            )
        for var in variables:
            if seen[var] != -1:
                problems.append(f"t={t}: variable {var} sits in two groups")
            seen[var] = group
    missing = tuple(var for var in range(n) if seen[var] < 0)
    if missing:
        problems.append(f"t={t}: variables {missing} are missing from the partition")

    dimension = placed.auxiliary.dimension
    for row_index, row in enumerate(placed.measurement_rows):
        if len(row) != n:
            problems.append(f"measurement row {row_index} has length {len(row)}, expected {n}")
            continue
        if any(entry.denominator != 1 for entry in row):
            problems.append(f"t={t}: measurement row {row_index} is not an integer vector")
        absolute = sum((abs(entry) for entry in row), start=Fraction(0))
        if absolute > dimension:
            problems.append(
                f"t={t}: measurement row {row_index} has absolute sum {absolute}, above D={dimension}"
            )

    if any(len(row) != n for row in placed.measurement_rows):
        return problems
    from .linear_algebra import RationalRowSpan

    restricted_rows = [
        [[row[index] for index in variables] for row in placed.measurement_rows]
        for variables in groups
    ]
    spans = [RationalRowSpan(rows, len(variables)) for rows, variables in zip(restricted_rows, groups)]

    position = t - 1
    for clause_index, clause in enumerate(placed.cnf.clauses):
        assigned = placed.clause_schedules[clause_index][position]
        for occ, (var, _neg) in enumerate(clause):
            coeff = [
                eval_poly(placed.family.rows[clause_index][occ][i], Fraction(t))
                for i in range(n)
            ]
            for group, variables in enumerate(groups):
                if group == assigned:
                    continue
                target = [coeff[index] for index in variables]
                if not spans[group].contains(target):
                    problems.append(
                        f"t={t}: clause {clause_index} occurrence {occ} (variable {var}) "
                        f"reads group {group} outside the measurement span"
                    )
                    continue
                if not boolean_check:
                    continue
                table: dict[tuple[Fraction, ...], Fraction] = {}
                for bits in itertools.product((0, 1), repeat=len(variables)):
                    measurement = tuple(
                        sum(
                            (row[index] * bit for index, bit in zip(variables, bits)),
                            start=Fraction(0),
                        )
                        for row in placed.measurement_rows
                    )
                    dot = sum(
                        (coeff[index] * bit for index, bit in zip(variables, bits)),
                        start=Fraction(0),
                    )
                    previous = table.get(measurement)
                    if previous is not None and previous != dot:
                        problems.append(
                            f"t={t}: clause {clause_index} occurrence {occ} has two Boolean "
                            f"contributions for one measurement of group {group}"
                        )
                        break
                    table[measurement] = dot
    return problems


def restoration_problems(placed: PlacedFormula) -> list[str]:
    """Coefficient-wise T = 0 check: the row of occurrence (e, v) is the basis vector e_v."""
    problems: list[str] = []
    for clause_index, clause in enumerate(placed.cnf.clauses):
        for occ, (var, _neg) in enumerate(clause):
            for i in range(placed.cnf.n):
                value = eval_poly(placed.family.rows[clause_index][occ][i], Fraction(0))
                expected = Fraction(1) if i == var else Fraction(0)
                if value != expected:
                    problems.append(
                        f"T=0: clause {clause_index} occurrence {occ} "
                        f"has coefficient {value} on X_{i}, expected {expected}"
                    )
                    return problems
            degree = max(len(poly) - 1 for poly in placed.family.rows[clause_index][occ])
            if degree > placed.auxiliary.R - 1:
                problems.append(
                    f"clause {clause_index} occurrence {occ} has degree {degree}, "
                    f"above R-1={placed.auxiliary.R - 1}"
                )
    return problems


@dataclass(frozen=True)
class DegreeCapReport:
    """Degree of H against the occurrence cap and the Step 1 shape of cap.

    ``cap_occurrences`` is M·(R-1). ``cap_step1`` is m·k·R. Both interpolants
    are compared with the brute-force count and with one extra node.
    ``direct_at_zero`` is ``eval_H_at`` at T = 0, the same sum with no interpolation.
    """

    h0_occurrences: Fraction
    h0_step1: Fraction
    direct_at_zero: Fraction
    brute: int
    observed_degree: int
    cap_occurrences: int
    cap_step1: int
    extra_occurrences: bool
    extra_step1: bool

    @property
    def ok(self) -> bool:
        target = Fraction(self.brute)
        return (
            self.h0_occurrences == target
            and self.h0_step1 == target
            and self.direct_at_zero == target
            and self.direct_at_zero.denominator == 1
            and self.observed_degree <= self.cap_occurrences
            and self.observed_degree <= self.cap_step1
            and self.extra_occurrences
            and self.extra_step1
        )


def _recover(samples: Sequence[Fraction], cap: int) -> tuple[Fraction, bool]:
    xs = [Fraction(t) for t in range(1, cap + 2)]
    ys = list(samples[: cap + 1])
    h0 = lagrange_eval_at_zero(xs, ys)
    predicted = lagrange_eval(xs, ys, Fraction(cap + 2))
    return h0, predicted == samples[cap + 1]


def degree_cap_report(placed: PlacedFormula) -> DegreeCapReport:
    """Sample H(1), H(2), ... and test both degree caps.

    The samples are the full-cube sum. This checks the cap; it does not speed
    up the sum. n has to stay inside the brute-force limit of ``eval_H_at``.
    """
    cap_occurrences = placed.cap_occurrences
    cap_step1 = placed.cap_step1
    window = max(cap_occurrences, cap_step1, 0) + 2
    samples = [eval_H_at(placed.cnf, placed.family, Fraction(t)) for t in range(1, window + 1)]
    h0_occurrences, extra_occurrences = _recover(samples, cap_occurrences)
    h0_step1, extra_step1 = _recover(samples, cap_step1)
    return DegreeCapReport(
        h0_occurrences=h0_occurrences,
        h0_step1=h0_step1,
        direct_at_zero=eval_H_at(placed.cnf, placed.family, Fraction(0)),
        brute=brute_count(placed.cnf),
        observed_degree=equispaced_degree(samples),
        cap_occurrences=cap_occurrences,
        cap_step1=cap_step1,
        extra_occurrences=extra_occurrences,
        extra_step1=extra_step1,
    )


@dataclass(frozen=True)
class Step3Example:
    """Named slice used by tests and by the verify CLI."""

    name: str
    notes: str
    placed: PlacedFormula
    sample_points: tuple[Fraction, ...]

    @property
    def auxiliary(self) -> AuxiliarySlice:
        return self.placed.auxiliary


def _two_positive_units() -> CNF:
    return CNF(2, [[(0, False)], [(1, False)]])


def step3_examples() -> list[Step3Example]:
    """Four slices. The first three are small enough to compare with a hand expansion.

    ``swap-4`` is the 4×4 exchange matrix at r=2, Q=1, R=1.
    ``remark-4-4`` is the four-coordinate row at r=2, Q=2, R=2.
    ``cross-r3`` puts one width-3 clause on three groups of size 1.
    ``r2-q3`` is the machine-checked point with R < Q and side length 64.
    """
    swap = build_auxiliary_slice(2, 1, 1)
    # Each clause is assigned to the other variable's group, so the degree-0
    # rows still equal X_v but the outside variable is carried by the measurement.
    swap_placed = place_formula(
        _two_positive_units(),
        swap,
        variable_schedules=((0,), (1,)),
        clause_schedules=((1,), (0,)),
    )

    remark = build_auxiliary_slice(2, 2, 2)
    remark_placed = place_formula(
        _two_positive_units(),
        remark,
        variable_schedules=((0, 0), (1, 0)),
        clause_schedules=((1, 1), (0, 1)),
    )

    cross = build_auxiliary_slice(3, 1, 1)
    cross_placed = place_formula(
        CNF(3, [[(0, False), (1, False), (2, False)]]),
        cross,
        variable_schedules=((0,), (1,), (2,)),
        clause_schedules=((0,),),
    )

    wider = build_auxiliary_slice(2, 3, 2)
    wider_placed = place_formula(
        CNF(3, [[(0, False), (1, False), (2, True)]]),
        wider,
        variable_schedules=((0, 1, 0), (1, 0, 0), (0, 0, 1)),
        clause_schedules=((1, 0, 1),),
    )

    return [
        Step3Example(
            name="swap-4",
            notes=(
                "r=2, Q=1, R=1. M0 is the 4×4 exchange. Each unit clause is assigned "
                "to the other variable's group, so the measurement carries X_v"
            ),
            placed=swap_placed,
            sample_points=(Fraction(0), Fraction(1), Fraction(2)),
        ),
        Step3Example(
            name="remark-4-4",
            notes=(
                "r=2, Q=2, R=2, so R=Q. Hand row ζ0-(T/2)ζ1-T ζ2 on the two unit clauses"
            ),
            placed=remark_placed,
            sample_points=(Fraction(0), Fraction(1), Fraction(2), Fraction(1, 2)),
        ),
        Step3Example(
            name="cross-r3",
            notes="r=3, Q=1, R=1. One 3-clause meets three groups; each group has size 1 ≤ 2n/r",
            placed=cross_placed,
            sample_points=(Fraction(0), Fraction(1), Fraction(2)),
        ),
        Step3Example(
            name="r2-q3",
            notes="r=2, Q=3, R=2, side length 64, R<Q. Checked by the machine, not expanded by hand",
            placed=wider_placed,
            sample_points=(Fraction(0), Fraction(1), Fraction(2)),
        ),
    ]
