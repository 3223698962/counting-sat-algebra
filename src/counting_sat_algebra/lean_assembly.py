"""Finite assembly and placement of one native leaf from the supplied program.

``NativeCodeAssembly`` wires one connected run: the scheduled searches yield
the shifts and local matrix, ``CodeRowsProgram`` serializes the regenerated
rows as framed payloads into the retained history slot, ``LeafRingSetupProgram``
places the ring constants into the fixed port map recorded by ``QueryReady``,
and the leaf appends greedily chosen rows before the final accepted query.
This finite port places the same objects with the same encodings and lengths:

- modulus bits (length <= 2t+1) at port 194, the d^3 structure words of
  ``X^(i+j) mod F`` at 201, the capacity width word at 202, scalar codes of
  0, 1, -1 and 1/2 at 203/207/211/253, the unary degree, precision, jet and
  ``3^d - 2`` exponent words at 204/212/243/209, the ``3t``-coefficient
  idempotent code (``6t^2`` bits) at 281, flat and framed node tables at
  256/300, and the interpolation caches (inverses, negatives, pairwise unit
  differences) at 713/714/715 with the unary widths at 716/717;
- the serialized row stream is framed exactly like ``Encoding.encodeList``,
  decodes back to the regenerated slot lists, and satisfies the
  ``HashGreedyCode.initial_encoding`` identity;
- one state is threaded through schedule, placement, serialize, decode and
  count, with a cross-path score equality between the decoded history and
  the direct rows.

The 896-stack physical machine, building the ring target from the query
heads, grouping replacement and the whole-program prefix-space refinement
remain open.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields

from .cnf import CNF
from .lean_code import CodeBudgets, construct_sparse_code
from .lean_encoding import count_word, decode_count_word
from .lean_ext import (
    CountingExtension,
    ExtensionBudgets,
    ExtensionRefusal,
    counting_extension,
    node_code_words,
    node_element,
)
from .lean_greedy import GreedyBudgets, count_cnf_greedy, occurrence_bound
from .lean_query import TernaryHashRow, count_ternary_query
from .lean_ring import (
    LeanAlgebraBudgets,
    construct_ternary_indicator,
    natural,
)
from .lean_stream import streamed_query_count
from .lean_taylor import TaylorBudgets

STEP30_SCOPE = (
    'Finite assembly and placement of one native leaf: ring constants, node '
    'tables and interpolation caches at the recorded QueryReady ports, framed '
    'row-stream serialization with round trips, and one threaded run from '
    'schedule to count word. The physical 896-stack machine, the query-derived '
    'ring target, grouping replacement and the prefix-space refinement remain '
    'open. No Theorem 1.1, no 2^{o(n)} claim.'
)

# The recorded QueryReady/LeafRingSetupProgram port indices this module places.
PORT = {
    'ambient': 191,
    'modulus': 194,
    'table': 201,
    'width': 202,
    'scalar_zero': 203,
    'degree': 204,
    'one': 207,
    'exponent': 209,
    'minus': 211,
    'precision': 212,
    'jet_order': 243,
    'half': 253,
    'flat_nodes': 256,
    'zero': 280,
    'idempotent': 281,
    'framed_nodes': 300,
    'inverse': 713,
    'negative': 714,
    'differences': 715,
    'row_width': 716,
    'ring_width': 717,
}


class AssemblyRefusal(ValueError):
    """An explicit limit stops the placement before any state is retained."""


@dataclass(frozen=True)
class AssemblyBudgets:
    max_extension_degree: int = 8
    max_exponent_unary: int = 20_000
    max_nodes: int = 4_096
    max_table_entries: int = 40_000
    max_stream_bits: int = 2_000_000

    def __post_init__(self) -> None:
        for item in fields(self):
            object.__setattr__(self, item.name, natural(getattr(self, item.name), item.name))


def _bits(value: int, width: int) -> tuple[int, ...]:
    value %= 1 << width
    return tuple((value >> i) & 1 for i in range(width))


def _minimal_bits(value: int) -> tuple[int, ...]:
    value = natural(value, 'value')
    return _bits(value, max(1, value.bit_length()))


def _ones(count: int) -> tuple[int, ...]:
    return (1,) * natural(count, 'count')


def _frame(bits: tuple[int, ...]) -> tuple[int, ...]:
    output: list[int] = []
    for bit in bits:
        output.extend((1, bit))
    output.append(0)
    return tuple(output)


def _encode_list(items: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    output: list[int] = []
    for item in items:
        output.append(1)
        output.extend(_frame(item))
    output.append(0)
    return tuple(output)


def _decode_list(bits: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    items: list[tuple[int, ...]] = []
    position = 0
    while position < len(bits) and bits[position] == 1:
        position += 1
        payload: list[int] = []
        while position < len(bits) and bits[position] == 1:
            if position + 1 >= len(bits):
                raise ValueError('incomplete escape frame')
            payload.append(bits[position + 1])
            position += 2
        if position >= len(bits):
            raise ValueError('unterminated escape frame')
        position += 1
        items.append(tuple(payload))
    if position != len(bits) - 1 or bits[position] != 0:
        raise ValueError('unterminated list or trailing bits')
    return tuple(items)


def _flat(extension: CountingExtension, element) -> tuple[int, ...]:
    output: list[int] = []
    for coefficient in element:
        output.extend(_bits(coefficient, 2 * extension.precision))
    return tuple(output)


def decode_flat(extension: CountingExtension, bits: tuple[int, ...], words: int):
    width = 2 * extension.precision
    if len(bits) != words * extension.degree * width:
        raise ValueError('the flat word count does not match words*d*2t')
    elements = []
    for index in range(words):
        start = index * extension.degree * width
        coefficients = [
            sum(bits[start + word * width + bit] << bit for bit in range(width))
            for word in range(extension.degree)
        ]
        elements.append(extension.element(coefficients))
    return tuple(elements)


def structure_table(extension: CountingExtension) -> tuple[tuple[int, ...], ...]:
    """The d^3 words ``coefficient k of X^(i+j) mod F``, row-major."""
    d = extension.degree
    words: list[tuple[int, ...]] = []
    for i in range(d):
        for j in range(d):
            power = [0] * d
            power[0] = 1
            for _ in range(i + j):
                shifted = [0] * d
                for index, coefficient in enumerate(power):
                    if not coefficient:
                        continue
                    if index + 1 < d:
                        shifted[index + 1] = (
                            shifted[index + 1] + coefficient
                        ) % extension.modulus
                    else:
                        for low, table_coefficient in enumerate(extension.coefficients):
                            shifted[low] = (
                                shifted[low] - coefficient * table_coefficient
                            ) % extension.modulus
                power = shifted
            for k in range(d):
                words.append(_bits(power[k], 2 * extension.precision))
    return tuple(words)


def payload_to_slots(payload: tuple[int, ...], width: int, m: int):
    if len(payload) % width:
        raise ValueError('the payload length is not a multiple of the width word')
    slots = []
    for start in range(0, len(payload), width):
        value = sum(payload[start + bit] << bit for bit in range(width))
        slots.append((value % m, value // m))
    return tuple(slots)


@dataclass(frozen=True)
class PlacedLeafState:
    ports: tuple[tuple[str, int, int], ...]  # (role, port index, bit length)
    state: dict = field(compare=False, repr=False)

    def port(self, role: str) -> tuple[int, ...]:
        for name, index, _length in self.ports:
            if name == role:
                return self.state[index]
        raise KeyError(role)

    @property
    def lengths(self) -> dict[str, int]:
        return {name: length for name, _index, length in self.ports}


def place_leaf_state(
    cnf: CNF,
    code,
    *,
    precision: int | None = None,
    degree: int | None = None,
    jet_order: int = 2,
    budgets: AssemblyBudgets | None = None,
    extension_budgets: ExtensionBudgets | None = None,
) -> PlacedLeafState:
    """Place one leaf's ring constants, nodes and caches at the recorded ports."""
    limits = budgets or AssemblyBudgets()
    n = natural(cnf.n, 'n')
    m = natural(cnf.m, 'clause count')
    rows = code.row_count
    if rows > limits.max_nodes:
        raise AssemblyRefusal('node budget exceeded')
    t = natural(precision, 'precision') if precision is not None else n + m.bit_length() + 1
    d = natural(degree, 'degree') if degree is not None else rows.bit_length() + 1
    if d > limits.max_extension_degree:
        raise AssemblyRefusal('extension degree budget exceeded')
    if rows >= 3**d:
        raise AssemblyRefusal('the node count must stay below 3^d')
    if 3**d - 2 > limits.max_exponent_unary:
        raise AssemblyRefusal('exponent unary budget exceeded')
    if d * d * d > limits.max_table_entries:
        raise AssemblyRefusal('structure table budget exceeded')
    jet = natural(jet_order, 'jet order')
    extension = counting_extension(t, d, budgets=extension_budgets)
    nodes = [node_element(extension, index) for index in range(1, rows + 1)]
    half = (extension.modulus + 1) // 2
    if 2 * half % extension.modulus != 1:
        raise ArithmeticError('the half word does not invert two')
    indicator = construct_ternary_indicator(t)
    state: dict = {}
    lengths: dict[str, int] = {}

    def put(role: str, bits: tuple[int, ...]) -> None:
        state[PORT[role]] = bits
        lengths[role] = len(bits)

    put('ambient', _ones(n))
    put('modulus', _minimal_bits(extension.modulus))
    put('table', tuple(bit for word in structure_table(extension) for bit in word))
    width_value = (3 * m).bit_length()
    if m and 3 * m > 1 << width_value:
        raise ArithmeticError('the width word misses the row capacity')
    put('width', _bits(width_value, 2 * t))
    put('scalar_zero', _bits(0, 2 * t))
    put('degree', _ones(d))
    put('one', _flat(extension, extension.one))
    put('exponent', _ones(3**d - 2))
    put('minus', _flat(extension, extension.neg(extension.one)))
    put('precision', _ones(t))
    put('jet_order', _ones(jet))
    put('half', _flat(extension, extension.constant(half)))
    put(
        'flat_nodes',
        tuple(
            bit
            for index in range(1, rows + 1)
            for word in node_code_words(extension, index)
            for bit in _bits(word, 2 * t)
        ),
    )
    put('zero', _flat(extension, extension.zero))
    put('idempotent', tuple(indicator.encode()))
    put('framed_nodes', _encode_list(tuple(_flat(extension, node) for node in nodes)))
    put(
        'inverse',
        tuple(bit for node in nodes for bit in _flat(extension, extension.inverse(node))),
    )
    put(
        'negative',
        tuple(bit for node in nodes for bit in _flat(extension, extension.neg(node))),
    )
    difference_bits: list[int] = []
    for left in nodes:
        for right in nodes:
            if left == right:
                difference_bits.extend(_bits(0, d * 2 * t))
            else:
                difference = extension.sub(left, right)
                difference_bits.extend(_flat(extension, extension.inverse(difference)))
    put('differences', tuple(difference_bits))
    put('row_width', _ones(rows * d * 2 * t))
    put('ring_width', _ones(d * 2 * t))
    ports = tuple((name, PORT[name], lengths[name]) for name in PORT)
    return PlacedLeafState(ports, state)


def serialize_leaf_history(code) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    """Frame the regenerated rows exactly like ``CodeRowsProgram.rowStream``."""
    width = (3 * code.m).bit_length()
    payloads = []
    for index in range(code.row_count):
        bits: list[int] = []
        for column, coefficient in code.row(index).terms:
            bits.extend(_bits(column + code.m * coefficient, width))
        payloads.append(tuple(bits))
    return _encode_list(tuple(payloads)), tuple(payloads)


def history_from_payloads(history_rows: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    """``HashGreedyProgram.encodedHistory``: encode stored history rows.

    A stored history row is a payload in reversed orientation, so the
    encoding reverses each row once more; for the initial rows this equals
    the raw ``rowStream`` (the ``initial_encoding`` identity).
    """
    return _encode_list(tuple(row[::-1] for row in history_rows))


def append_payload(history: tuple[int, ...], payload: tuple[int, ...]) -> tuple[int, ...]:
    """``HashGreedyProgram`` append: ``true :: frame (payload.reverse) ++ history``."""
    return (1,) + _frame(payload[::-1]) + history


def decode_history(history: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return _decode_list(history)


@dataclass(frozen=True)
class AssembledLeafResult:
    status: str
    count: int | None
    reason: str
    n: int = 0
    m: int = 0
    row_count: int = 0
    precision: int = 0
    degree: int = 0
    jet_order: int = 0
    stream_bits: int = 0
    port_bits: int = 0
    rows_round_trip: bool = False
    score_paths_equal: bool = False
    rounds: int = 0
    count_word_bits: int = 0
    same_declared_n: bool = True
    polynomial_terms: int = 0
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False
    recursive_contraction: bool = False

    @property
    def counted(self) -> bool:
        return self.status == 'counted'


def run_assembled_leaf(
    cnf: CNF,
    *,
    depth: int = 1,
    occurrence: int | None = None,
    jet_order: int = 2,
    budgets: AssemblyBudgets | None = None,
    code_budgets: CodeBudgets | None = None,
    greedy_budgets: GreedyBudgets | None = None,
    taylor_budgets: TaylorBudgets | None = None,
    native_budgets: LeanAlgebraBudgets | None = None,
) -> AssembledLeafResult:
    """One threaded run: schedule, placement, framed stream, decode and count."""
    limits = budgets or AssemblyBudgets()
    n = natural(cnf.n, 'n')
    m = natural(cnf.m, 'clause count')
    try:
        code = None
        if m:
            if occurrence is None:
                occurrence = occurrence_bound(cnf)
            assembled = construct_sparse_code(
                m, natural(occurrence, 'occurrence'), budgets=code_budgets
            )
            if not assembled.checked:
                raise AssemblyRefusal(f'the scheduled sparse code refused: {assembled.reason}')
            code = assembled.certificate.code
        state = (
            place_leaf_state(cnf, code, jet_order=jet_order, budgets=limits)
            if code is not None
            else None
        )
        stream: tuple[int, ...] = ()
        payloads: tuple[tuple[int, ...], ...] = ()
        rows_round_trip = True
        if code is not None:
            stream, payloads = serialize_leaf_history(code)
            if len(stream) > limits.max_stream_bits:
                raise AssemblyRefusal('row stream bit budget exceeded')
            stored = decode_history(stream)
            if stored != payloads:
                raise ArithmeticError('the framed stream failed its round trip')
            history_rows = tuple(payload[::-1] for payload in payloads)
            if history_from_payloads(history_rows) != stream:
                raise ArithmeticError('the initial encoding identity failed')
            width = (3 * m).bit_length()
            decoded_slots = tuple(
                payload_to_slots(payload, width, m) for payload in stored
            )
            direct_slots = tuple(code.row(index).terms for index in range(code.row_count))
            rows_round_trip = decoded_slots == direct_slots
            if not rows_round_trip:
                raise ArithmeticError('the decoded slots differ from the regenerated rows')
        score_paths_equal = True
        if m and code is not None:
            decoded_rows = tuple(TernaryHashRow(terms) for terms in decoded_slots)
            direct_rows = tuple(code.rows)
            streamed = streamed_query_count(
                cnf, decoded_rows, marked_clause=0, budgets=taylor_budgets
            )
            native = count_ternary_query(cnf, direct_rows, marked_clause=0)
            score_paths_equal = (
                streamed.counted and native.counted and streamed.count == native.count
            )
            if not score_paths_equal:
                raise ArithmeticError('the decoded history changed the marked score')
        greedy = count_cnf_greedy(
            cnf,
            depth=depth,
            occurrence=occurrence,
            budgets=greedy_budgets,
            native_budgets=native_budgets,
        )
        if not greedy.counted:
            raise AssemblyRefusal(f'the assembled leaf refused: {greedy.reason}')
        word = count_word(greedy.count)
        if decode_count_word(word) != greedy.count or len(word) > n + 1:
            raise ArithmeticError('the count word failed its round trip')
        return AssembledLeafResult(
            'counted',
            greedy.count,
            'assembled placement and framed stream threaded through one run',
            n=n,
            m=m,
            row_count=0 if code is None else code.row_count,
            precision=0 if state is None else state.lengths.get('precision', 0),
            degree=0 if code is None else len(code.local_matrix),
            jet_order=jet_order,
            stream_bits=len(stream),
            port_bits=0 if state is None else sum(length for _n2, _i, length in state.ports),
            rows_round_trip=rows_round_trip,
            score_paths_equal=score_paths_equal,
            rounds=greedy.rounds,
            count_word_bits=len(word),
        )
    except (AssemblyRefusal, ExtensionRefusal) as exc:
        return AssembledLeafResult(
            'refused', None, f'{exc}; no partial count returned', n=n, m=m
        )
