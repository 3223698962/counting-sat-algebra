"""Step 30 finite witnesses for leaf assembly and placement."""

from __future__ import annotations

from .protocol import _row, stamp_completion
from ..cnf import CNF, brute_count
from ..lean_assembly import (
    PORT,
    STEP30_SCOPE,
    AssemblyBudgets,
    append_payload,
    decode_flat,
    decode_history,
    history_from_payloads,
    payload_to_slots,
    place_leaf_state,
    run_assembled_leaf,
    serialize_leaf_history,
    structure_table,
)
from ..lean_code import construct_sparse_code
from ..lean_ext import counting_extension, node_code_words, node_element
from ..lean_greedy import occurrence_bound
from ..lean_ring import construct_ternary_indicator

STAR = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
HUB = CNF(8, [[(0, False), (i, False)] for i in range(1, 8)])
FOREST = CNF(
    6,
    [
        [(0, False), (1, False)],
        [(0, False), (2, False)],
        [(3, False), (4, False)],
        [(3, False), (5, False)],
    ],
)
DEMO = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])
CASES = ((STAR, 9), (HUB, 129), (FOREST, 25), (DEMO, 196))


def _frame_reference(bits):
    output = []
    for bit in bits:
        output.extend((1, bit))
    output.append(0)
    return tuple(output)


def _code_for(cnf):
    return construct_sparse_code(cnf.m, occurrence_bound(cnf)).certificate.code


def run_lean_assembly_suite() -> list[dict]:
    print('=== (AF) Leaf assembly and placement ===', flush=True)
    rows = []

    lengths_ok = True
    details = []
    for cnf, _expected in CASES:
        code = _code_for(cnf)
        state = place_leaf_state(cnf, code)
        t = state.lengths['precision']
        d = state.lengths['degree']
        q = code.row_count
        valid = (
            state.lengths['ambient'] == cnf.n
            and state.lengths['modulus'] <= 2 * t + 1
            and state.lengths['table'] == d**3 * 2 * t
            and state.lengths['width'] == 2 * t
            and state.lengths['scalar_zero'] == 2 * t
            and state.lengths['zero'] == d * 2 * t
            and state.lengths['one'] == d * 2 * t
            and state.lengths['minus'] == d * 2 * t
            and state.lengths['half'] == d * 2 * t
            and state.lengths['degree'] == d
            and state.lengths['precision'] == t
            and state.lengths['jet_order'] == 2
            and state.lengths['exponent'] == 3**d - 2
            and state.lengths['flat_nodes'] == q * d * 2 * t
            and state.lengths['idempotent'] == 6 * t * t
            and state.lengths['framed_nodes'] == 1 + (4 * d * t + 2) * q
            and state.lengths['inverse'] == q * d * 2 * t
            and state.lengths['negative'] == q * d * 2 * t
            and state.lengths['differences'] == q * q * d * 2 * t
            and state.lengths['row_width'] == q * d * 2 * t
            and state.lengths['ring_width'] == d * 2 * t
        )
        lengths_ok = lengths_ok and valid
        details.append(
            dict(
                n=cnf.n,
                m=cnf.m,
                t=t,
                d=d,
                nodes=q,
                port_bits=sum(length for _n2, _i, length in state.ports),
            )
        )
    rows.append(
        _row(
            'lean-assembly-placement-port-lengths',
            0,
            0,
            lengths_ok,
            'CHECKED',
            count=None,
            cases=details,
        )
    )

    demo_code = _code_for(DEMO)
    small_state = place_leaf_state(DEMO, demo_code, precision=2, degree=2)
    extension = counting_extension(2, 2)
    word_width = 2 * small_state.lengths['precision']
    placed_words = [
        sum(small_state.port('table')[index * word_width + b] << b for b in range(word_width))
        for index in range(small_state.lengths['table'] // word_width)
    ]
    table_ok = len(placed_words) == 8 and placed_words[0] == 1 and placed_words[1] == 0
    for i in range(2):
        for j in range(2):
            power = extension.one
            for _ in range(i + j):
                power = extension.mul(power, extension.element([0, 1]))
            for k in range(2):
                table_ok = table_ok and placed_words[(i * 2 + j) * 2 + k] == power[k]
    table_ok = table_ok and structure_table(extension) == tuple(
        tuple((word >> b) & 1 for b in range(word_width)) for word in placed_words
    )
    rows.append(
        _row(
            'lean-assembly-structure-table-words',
            0,
            0,
            table_ok,
            'CHECKED',
            count=None,
            words=len(placed_words),
        )
    )

    scalar_ok = True
    for cnf, code in ((STAR, _code_for(STAR)), (HUB, _code_for(HUB))):
        state = place_leaf_state(cnf, code)
        t = state.lengths['precision']
        d = state.lengths['degree']
        ext = counting_extension(t, d)
        scalar_ok = (
            scalar_ok
            and decode_flat(ext, state.port('zero'), 1) == (ext.zero,)
            and decode_flat(ext, state.port('one'), 1) == (ext.one,)
            and decode_flat(ext, state.port('minus'), 1) == (ext.neg(ext.one),)
        )
        (half,) = decode_flat(ext, state.port('half'), 1)
        scalar_ok = scalar_ok and ext.mul(half, ext.constant(2)) == ext.one
        scalar_ok = scalar_ok and state.port('idempotent') == tuple(
            construct_ternary_indicator(t).encode()
        )
        nodes = [node_element(ext, index) for index in range(1, code.row_count + 1)]
        expected_flat = tuple(
            bit
            for index in range(1, code.row_count + 1)
            for word in node_code_words(ext, index)
            for bit in tuple((word >> b) & 1 for b in range(2 * t))
        )
        scalar_ok = scalar_ok and state.port('flat_nodes') == expected_flat
        expected_frames = tuple(
            tuple(
                bit
                for word in node_code_words(ext, index)
                for bit in tuple((word >> b) & 1 for b in range(2 * t))
            )
            for index in range(1, code.row_count + 1)
        )
        scalar_ok = scalar_ok and decode_history(state.port('framed_nodes')) == expected_frames
        inverses = decode_flat(ext, state.port('inverse'), code.row_count)
        negatives = decode_flat(ext, state.port('negative'), code.row_count)
        scalar_ok = scalar_ok and all(
            ext.mul(node, inverse) == ext.one for node, inverse in zip(nodes, inverses)
        )
        scalar_ok = scalar_ok and negatives == tuple(ext.neg(node) for node in nodes)
        differences = decode_flat(ext, state.port('differences'), code.row_count**2)
        expected_differences = []
        for left in nodes:
            for right in nodes:
                if left == right:
                    expected_differences.append(ext.zero)
                else:
                    expected_differences.append(ext.inverse(ext.sub(left, right)))
        scalar_ok = scalar_ok and differences == tuple(expected_differences)
    rows.append(
        _row(
            'lean-assembly-scalar-codes-and-node-caches',
            0,
            0,
            scalar_ok,
            'CHECKED',
            count=None,
        )
    )

    stream_ok = True
    for cnf, _expected in CASES:
        code = _code_for(cnf)
        stream, payloads = serialize_leaf_history(code)
        stored = decode_history(stream)
        width = (3 * cnf.m).bit_length()
        stream_ok = (
            stream_ok
            and stored == payloads
            and history_from_payloads(tuple(payload[::-1] for payload in payloads)) == stream
            and len(stream) == 1 + sum(2 * len(payload) + 2 for payload in payloads)
            and all(
                payload_to_slots(payload, width, cnf.m) == code.row(index).terms
                for index, payload in enumerate(payloads)
            )
        )
    rows.append(
        _row(
            'lean-assembly-row-stream-round-trip',
            0,
            0,
            stream_ok,
            'CHECKED',
            count=None,
        )
    )

    code = _code_for(STAR)
    stream, payloads = serialize_leaf_history(code)
    reference = (1,) + _frame_reference(payloads[0][::-1]) + stream
    extended = append_payload(stream, payloads[0])
    appended_ok = (
        extended == reference and decode_history(extended) == (payloads[0][::-1],) + payloads
    )
    twice = append_payload(extended, payloads[1])
    appended_ok = (
        appended_ok
        and decode_history(twice)
        == (
            payloads[1][::-1],
            payloads[0][::-1],
        )
        + payloads
    )
    rows.append(
        _row(
            'lean-assembly-append-encoding',
            STAR.n,
            STAR.m,
            appended_ok,
            'CHECKED',
            count=None,
            appended_bits=len(extended) - len(stream),
        )
    )

    connected_ok = True
    connected_details = []
    for cnf, expected in CASES:
        result = run_assembled_leaf(cnf)
        valid = (
            result.counted
            and result.count == brute_count(cnf) == expected
            and result.rows_round_trip
            and result.score_paths_equal
            and result.count_word_bits <= cnf.n + 1
            and result.port_bits > 0
            and result.stream_bits > 0
            and result.polynomial_terms == 0
            and result.same_declared_n
        )
        connected_ok = connected_ok and valid
        connected_details.append(
            dict(
                n=cnf.n,
                m=cnf.m,
                count=result.count,
                rows=result.row_count,
                rounds=result.rounds,
                port_bits=result.port_bits,
            )
        )
    rows.append(
        _row(
            'lean-assembly-connected-run-counts',
            0,
            0,
            connected_ok,
            'COUNTED',
            count=None,
            cases=connected_details,
        )
    )

    refusals = 0
    for budgets in (
        AssemblyBudgets(max_nodes=1),
        AssemblyBudgets(max_extension_degree=1),
        AssemblyBudgets(max_exponent_unary=6),
        AssemblyBudgets(max_table_entries=8),
        AssemblyBudgets(max_stream_bits=2),
    ):
        result = run_assembled_leaf(STAR, budgets=budgets)
        if result.status == 'refused' and result.count is None:
            refusals += 1
    oversized = CNF(16, [[(index, False)] for index in range(16)])
    big_code = _code_for(oversized)
    try:
        place_leaf_state(oversized, big_code, degree=1)
    except Exception as exc:
        if '3^d' in str(exc):
            refusals += 1
    rows.append(
        _row(
            'lean-assembly-budget-refusals',
            0,
            0,
            refusals == 6,
            'REFUSED',
            count=None,
            refusals=refusals,
        )
    )

    terminal = (
        run_assembled_leaf(CNF(8, [])).count == 256
        and run_assembled_leaf(CNF(8, [[]])).count == 0
        and run_assembled_leaf(CNF(1, [[(0, False)]])).count == 1
    )
    rows.append(
        _row(
            'lean-assembly-terminal-leaves',
            0,
            0,
            terminal,
            'COUNTED',
            count=None,
            free_bits=256,
        )
    )

    port_ok = (
        PORT['modulus'] == 194
        and PORT['table'] == 201
        and PORT['flat_nodes'] == 256
        and PORT['idempotent'] == 281
        and PORT['framed_nodes'] == 300
        and PORT['inverse'] == 713
        and PORT['negative'] == 714
        and PORT['differences'] == 715
        and PORT['row_width'] == 716
        and PORT['ring_width'] == 717
    )
    rows.append(
        _row(
            'lean-assembly-recorded-port-map',
            0,
            0,
            port_ok,
            'CHECKED',
            count=None,
            ports=len(PORT),
        )
    )

    for row in rows:
        stamp_completion(row)
        print(f'  [AF] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP30_SCOPE}', flush=True)
    return rows
