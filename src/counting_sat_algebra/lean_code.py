"""Streamed sparse ternary code construction from the supplied Lean program.

The scheduled dimensions are exactly SparseCodeParameters/SparseCodeBounds.
Searches retain one matrix or shift tuple. Local checking visits just the
potential counterexamples of weight below the required distance. Rows are
regenerated from incidence slots, with zero coefficients and parallel ports
retained. Logical storage/operation receipts do not prove Python bit bounds.
"""

from __future__ import annotations

import operator
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, fields
from fractions import Fraction
from itertools import combinations, product

from .lean_query import TernaryHashRow
from .lean_ring import natural


class CodeRefusal(ValueError):
    """A limit prevents completing a constructor or independent checker."""


@dataclass(frozen=True)
class CodeBudgets:
    max_matrix_entries: int | None = 262_144
    max_dimension: int | None = 24
    max_candidates: int | None = 1_000_000
    max_word_checks: int | None = 2_000_000
    max_operations: int | None = 8_000_000
    max_code_bits: int | None = 16_000_000

    def __post_init__(self) -> None:
        for item in fields(self):
            value = getattr(self, item.name)
            if value is not None:
                object.__setattr__(self, item.name, natural(value, item.name))


def _limit(value: int, cap: int | None, reason: str) -> None:
    if cap is not None and value > cap:
        raise CodeRefusal(reason)


@dataclass
class CodeMeter:
    budgets: CodeBudgets
    candidates: int = 0
    word_checks: int = 0
    character_checks: int = 0
    operations: int = 0
    peak_bits: int = 0

    def candidate(self) -> None:
        _limit(self.candidates + 1, self.budgets.max_candidates, 'code candidate budget exceeded')
        self.candidates += 1

    def word(self) -> None:
        _limit(self.word_checks + 1, self.budgets.max_word_checks, 'code word budget exceeded')
        self.word_checks += 1

    def work(self, number: int) -> None:
        _limit(self.operations + number, self.budgets.max_operations, 'code operation budget exceeded')
        self.operations += number

    def storage(self, number: int) -> None:
        _limit(number, self.budgets.max_code_bits, 'code encoding bit budget exceeded')
        self.peak_bits = max(self.peak_bits, number)


@dataclass(frozen=True)
class SparseCodeParameters:
    m: int
    B: int
    log_size: int
    gap_denominator: int
    degree_bound: int
    degree: int
    dimension: int
    local_rows: int
    local_distance: int
    distance_denominator: int

    @property
    def vertices(self) -> int:
        return 1 << self.dimension

    @property
    def total_rows(self) -> int:
        return 2 * self.vertices * self.local_rows


def sparse_code_parameters(m: int, B: int) -> SparseCodeParameters:
    """Exact Lean integers, including floor(D/(8*B))+1 rather than a ceiling."""
    m, B = natural(m, 'm'), natural(B, 'B')
    if not m or not B:
        raise ValueError('the nonempty sparse code requires positive m and B')
    logarithm = m.bit_length()
    gap = 200 * B * logarithm
    bound = 8 * gap * gap * logarithm
    degree = min(m, bound)
    vertices_needed = (m + degree - 1) // degree
    dimension = (vertices_needed - 1).bit_length()
    return SparseCodeParameters(
        m, B, logarithm, gap, bound, degree, dimension, degree // (8 * B) + 1,
        degree // (100 * B * logarithm) + 1, 20_000 * B * B * logarithm * logarithm,
    )


def radix_words(radix: int, length: int) -> Iterator[tuple[int, ...]]:
    """An odometer; the first coordinate is least significant, no word table."""
    radix, length = natural(radix, 'radix'), natural(length, 'length')
    if not radix:
        raise ValueError('radix must be positive')
    digits = [0] * length
    while True:
        yield tuple(digits)
        for position in range(length):
            digits[position] += 1
            if digits[position] < radix:
                break
            digits[position] = 0
        else:
            return


def _matrix(matrix: Sequence[Sequence[int]], columns: int) -> tuple[tuple[int, ...], ...]:
    result = tuple(tuple(operator.index(c) for c in row) for row in matrix)
    if any(len(row) != columns for row in result):
        raise ValueError('every local row must have exactly D columns')
    if any(c not in (0, 1, 2) for row in result for c in row):
        raise ValueError('local coefficients must be canonical ternary residues')
    return result


def _local_distance(
    matrix: tuple[tuple[int, ...], ...], columns: int, distance: int, meter: CodeMeter,
) -> bool:
    # HasDistance is automatic for distance <= 1. Higher-weight words cannot
    # falsify the predicate, so no complete 3**D word table is needed.
    for weight in range(1, min(columns, distance - 1) + 1):
        for support in combinations(range(columns), weight):
            for coefficients in product((1, 2), repeat=weight):
                meter.word()
                killed = False
                for row in matrix:
                    meter.work(2 * weight)
                    if sum(row[i] * c for i, c in zip(support, coefficients)) % 3:
                        killed = True
                        break
                if not killed:
                    return False
    return True


def _walsh_good(dimension: int, shifts: tuple[int, ...], denominator: int, meter: CodeMeter) -> bool:
    for character in range(1, 1 << dimension):
        meter.work(2 * len(shifts) + 1)
        meter.character_checks += 1
        value = sum(-1 if (character & shift).bit_count() % 2 else 1 for shift in shifts)
        if denominator * abs(value) > len(shifts):
            return False
    return True


def _prepare(columns: int, rows: int, dimension: int, meter: CodeMeter) -> None:
    _limit(columns * rows, meter.budgets.max_matrix_entries, 'local matrix entry budget exceeded')
    _limit(dimension, meter.budgets.max_dimension, 'graph dimension budget exceeded')
    # Three matrix/counter copies, retained shifts, and one low-weight checker.
    meter.storage(6 * rows * columns + 3 * max(1, dimension) * columns
                  + 4 * columns * (columns.bit_length() + 2) + 64)


def _search_local(columns: int, rows: int, distance: int, meter: CodeMeter):
    _prepare(columns, rows, 0, meter)
    for word in radix_words(3, rows * columns):
        meter.candidate()
        candidate = tuple(tuple(word[i * columns:(i + 1) * columns]) for i in range(rows))
        if _local_distance(candidate, columns, distance, meter):
            return candidate
    return None


def _search_shifts(dimension: int, degree: int, denominator: int, meter: CodeMeter):
    _prepare(degree, 0, dimension, meter)
    for candidate in radix_words(1 << dimension, degree):
        meter.candidate()
        if _walsh_good(dimension, candidate, denominator, meter):
            return candidate
    return None


@dataclass(frozen=True)
class CodeSearchResult:
    status: str
    witness: tuple | None
    reason: str
    candidates: int
    word_checks: int
    character_checks: int
    operations: int
    peak_bits: int


def _search_result(status: str, witness, reason: str, meter: CodeMeter) -> CodeSearchResult:
    return CodeSearchResult(
        status, witness, reason, meter.candidates, meter.word_checks,
        meter.character_checks, meter.operations, meter.peak_bits,
    )


def search_local_ternary_matrix(
    columns: int, rows: int, distance: int, *, budgets: CodeBudgets | None = None,
) -> CodeSearchResult:
    columns, rows, distance = (natural(v, name) for v, name in
                              ((columns, 'columns'), (rows, 'rows'), (distance, 'distance')))
    meter = CodeMeter(budgets or CodeBudgets())
    try:
        matrix = _search_local(columns, rows, distance, meter)
        return _search_result('checked' if matrix is not None else 'not_found', matrix,
                              'streamed complete local-code search', meter)
    except CodeRefusal as error:
        return _search_result('refused', None, str(error), meter)


def search_walsh_shifts(
    dimension: int, degree: int, denominator: int, *, budgets: CodeBudgets | None = None,
) -> CodeSearchResult:
    dimension, degree, denominator = (natural(v, name) for v, name in
                                     ((dimension, 'dimension'), (degree, 'degree'),
                                      (denominator, 'denominator')))
    if not denominator:
        raise ValueError('Walsh denominator must be positive')
    meter = CodeMeter(budgets or CodeBudgets())
    try:
        shifts = _search_shifts(dimension, degree, denominator, meter)
        return _search_result('checked' if shifts is not None else 'not_found', shifts,
                              'streamed complete Walsh search', meter)
    except CodeRefusal as error:
        return _search_result('refused', None, str(error), meter)


@dataclass(frozen=True)
class SparseTannerCode:
    m: int
    dimension: int
    shifts: tuple[int, ...]
    local_matrix: tuple[tuple[int, ...], ...]
    gap_denominator: int
    local_distance: int

    def __post_init__(self) -> None:
        for name in ('m', 'dimension', 'gap_denominator', 'local_distance'):
            object.__setattr__(self, name, natural(getattr(self, name), name))
        shifts = tuple(natural(s, 'shift') for s in self.shifts)
        object.__setattr__(self, 'shifts', shifts)
        object.__setattr__(self, 'local_matrix', _matrix(self.local_matrix, len(shifts)))
        if not self.m or not shifts or not self.gap_denominator:
            raise ValueError('Tanner m, degree and Walsh denominator must be positive')
        if any(s.bit_length() > self.dimension for s in shifts):
            raise ValueError('a shift is outside its binary cube')
        if self.m > self.degree * (1 << self.dimension):
            raise ValueError('the edge domain does not cover the original m columns')

    @property
    def degree(self) -> int:
        return len(self.shifts)

    @property
    def row_count(self) -> int:
        return 2 * (1 << self.dimension) * len(self.local_matrix)

    def row(self, index: int) -> TernaryHashRow:
        index = operator.index(index)
        if index < 0:
            index += self.row_count
        if index < 0 or index >= self.row_count:
            raise IndexError('sparse code row index out of range')
        chronological = self.row_count - 1 - index
        local_count = len(self.local_matrix)
        vertex, within = divmod(chronological, 2 * local_count)
        side, local_index = divmod(within, local_count)
        output = []
        for port, coefficient in enumerate(self.local_matrix[local_index]):
            source = vertex ^ self.shifts[port] if side else vertex
            column = source * self.degree + port
            if column < self.m:
                output.append((column, coefficient))
        return TernaryHashRow(tuple(output))

    @property
    def rows(self) -> SparseCodeRows:
        return SparseCodeRows(self)


class SparseCodeRows(Sequence[TernaryHashRow]):
    """Only the retained code and one regenerated row persist."""

    def __init__(self, code: SparseTannerCode):
        self.code = code

    def __len__(self) -> int:
        return self.code.row_count

    def __getitem__(self, index):
        if isinstance(index, slice):
            return tuple(self.code.row(i) for i in range(*index.indices(len(self))))
        return self.code.row(index)


@dataclass(frozen=True)
class SparseCodeCertificate:
    code: SparseTannerCode
    relative_distance: Fraction
    minimum_weight: int
    depth: int
    row_term_bound: int
    column_incidence_bound: int


def certify_sparse_code(code: SparseTannerCode, meter: CodeMeter | None = None) -> SparseCodeCertificate:
    """Recheck supplied candidates; no externally asserted distance flag is trusted."""
    used = meter or CodeMeter(CodeBudgets())
    _prepare(code.degree, len(code.local_matrix), code.dimension, used)
    if not _local_distance(code.local_matrix, code.degree, code.local_distance, used):
        raise ValueError('the local matrix fails its claimed distance')
    if not _walsh_good(code.dimension, code.shifts, code.gap_denominator, used):
        raise ValueError('the supplied shifts fail the exact Walsh bound')
    delta = Fraction(code.local_distance, code.degree)
    mixing = Fraction(1, code.gap_denominator)
    if delta <= mixing:
        raise ValueError('local relative distance must exceed the mixing error')
    relative = delta * (delta - mixing)
    # Matrix-relativeDistance plus zero-shortening: weight >= m * relative.
    weight = max(1, (code.m * relative.numerator + relative.denominator - 1)
                 // relative.denominator)
    depth = (code.m + weight - 1) // weight
    return SparseCodeCertificate(code, relative, weight, depth, code.degree,
                                 2 * len(code.local_matrix))


@dataclass(frozen=True)
class SparseCodeResult:
    status: str
    certificate: SparseCodeCertificate | None
    reason: str
    parameters: SparseCodeParameters
    candidates: int
    word_checks: int
    character_checks: int
    operations: int
    peak_bits: int

    @property
    def checked(self) -> bool:
        return self.status == 'checked'


def construct_sparse_code(m: int, B: int, *, budgets: CodeBudgets | None = None) -> SparseCodeResult:
    """The actual two scheduled searches, followed by checked slot assembly."""
    parameters = sparse_code_parameters(m, B)
    meter = CodeMeter(budgets or CodeBudgets())

    def finish(status, certificate, reason):
        return SparseCodeResult(status, certificate, reason, parameters, meter.candidates,
                                meter.word_checks, meter.character_checks, meter.operations,
                                meter.peak_bits)

    try:
        _prepare(parameters.degree, parameters.local_rows, parameters.dimension, meter)
        shifts = _search_shifts(parameters.dimension, parameters.degree,
                                parameters.gap_denominator, meter)
        matrix = _search_local(parameters.degree, parameters.local_rows,
                               parameters.local_distance, meter)
        if shifts is None or matrix is None:
            raise ArithmeticError('a complete scheduled code search found no promised witness')
        code = SparseTannerCode(m, parameters.dimension, shifts, matrix,
                                parameters.gap_denominator, parameters.local_distance)
        certificate = certify_sparse_code(code, meter)
        if parameters.distance_denominator * certificate.minimum_weight < m:
            raise ArithmeticError('the assembled code fails the scheduled kernel lower bound')
        if certificate.depth > parameters.distance_denominator:
            raise ArithmeticError('the optimized hash depth exceeds the Lean depth')
        return finish('checked', certificate, 'both scheduled witnesses checked; rows regenerate')
    except CodeRefusal as error:
        return finish('refused', None, f'{error}; no partial sparse code returned')


def canonical_hash_candidates(m: int, depth: int) -> Iterator[tuple[int, TernaryHashRow]]:
    """Every ternary tuple once, with its first physical binary-odometer rank.

    Normalization only decreases each raw binary digit. Its first occurrence
    is consequently the canonical word itself. Omitting later duplicates
    preserves the full natural-score minimum and first-minimum tie rule.
    """
    m, depth = natural(m, 'm'), natural(depth, 'depth')
    if not m:
        raise ValueError('hash candidates require positive m')
    width = (3 * m).bit_length()
    for digits in radix_words(3 * m, depth):
        rank = sum(digit << (width * index) for index, digit in enumerate(digits))
        terms = tuple((digit % m, digit // m) for digit in digits)
        yield rank, TernaryHashRow(terms)
