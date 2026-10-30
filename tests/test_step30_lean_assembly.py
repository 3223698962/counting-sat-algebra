"""Step 30 independent checks for leaf assembly and placement."""

import pytest

from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.lean_assembly import (
    PORT,
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
from counting_sat_algebra.lean_code import construct_sparse_code
from counting_sat_algebra.lean_ext import counting_extension, node_element
from counting_sat_algebra.lean_greedy import count_cnf_greedy, occurrence_bound
from counting_sat_algebra.lean_ring import construct_ternary_indicator

STAR = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
HUB = CNF(8, [[(0, False), (i, False)] for i in range(1, 8)])
FOREST = CNF(
    6,
    [[(0, False), (1, False)], [(0, False), (2, False)],
     [(3, False), (4, False)], [(3, False), (5, False)]],
)
DEMO = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])
CASES = ((STAR, 9), (HUB, 129), (FOREST, 25), (DEMO, 196))


def _code_for(cnf):
    return construct_sparse_code(cnf.m, occurrence_bound(cnf)).certificate.code


@pytest.mark.parametrize('cnf,expected', CASES)
def test_placement_lengths_match_the_recorded_formulas(cnf, expected):
    code = _code_for(cnf)
    state = place_leaf_state(cnf, code)
    t, d, q = state.lengths['precision'], state.lengths['degree'], code.row_count
    assert state.lengths['ambient'] == cnf.n
    assert state.lengths['modulus'] <= 2 * t + 1
    assert state.lengths['table'] == d**3 * 2 * t
    for role in ('width', 'scalar_zero'):
        assert state.lengths[role] == 2 * t
    for role in ('zero', 'one', 'minus', 'half'):
        assert state.lengths[role] == d * 2 * t
    assert state.lengths['degree'] == d
    assert state.lengths['precision'] == t
    assert state.lengths['jet_order'] == 2
    assert state.lengths['exponent'] == 3**d - 2
    assert state.lengths['flat_nodes'] == q * d * 2 * t
    assert state.lengths['idempotent'] == 6 * t * t
    assert state.lengths['framed_nodes'] == 1 + (4 * d * t + 2) * q
    assert state.lengths['inverse'] == q * d * 2 * t
    assert state.lengths['negative'] == q * d * 2 * t
    assert state.lengths['differences'] == q * q * d * 2 * t
    assert state.lengths['row_width'] == q * d * 2 * t
    assert state.lengths['ring_width'] == d * 2 * t


def test_structure_table_words_match_the_power_oracle():
    code = _code_for(DEMO)
    state = place_leaf_state(DEMO, code, precision=2, degree=2)
    extension = counting_extension(2, 2)
    assert extension.coefficients == (1, 0)
    word_width = 4
    placed = [
        sum(state.port('table')[index * word_width + b] << b for b in range(word_width))
        for index in range(state.lengths['table'] // word_width)
    ]
    assert len(placed) == 8 and placed[0] == 1 and placed[1] == 0
    for i in range(2):
        for j in range(2):
            power = extension.one
            for _ in range(i + j):
                power = extension.mul(power, extension.element([0, 1]))
            for k in range(2):
                assert placed[(i * 2 + j) * 2 + k] == power[k]
    assert structure_table(extension) == tuple(
        tuple((word >> b) & 1 for b in range(word_width)) for word in placed
    )


@pytest.mark.parametrize('precision', [2, 3])
def test_scalar_codes_and_idempotent_placement(precision):
    code = _code_for(DEMO)
    state = place_leaf_state(DEMO, code, precision=precision, degree=2)
    extension = counting_extension(precision, 2)
    assert decode_flat(extension, state.port('zero'), 1) == (extension.zero,)
    assert decode_flat(extension, state.port('one'), 1) == (extension.one,)
    assert decode_flat(extension, state.port('minus'), 1) == (extension.neg(extension.one),)
    (half,) = decode_flat(extension, state.port('half'), 1)
    assert extension.mul(half, extension.constant(2)) == extension.one
    assert state.port('idempotent') == tuple(construct_ternary_indicator(precision).encode())
    assert len(state.port('idempotent')) == 6 * precision * precision


def test_node_tables_and_caches_match_independent_oracles():
    code = _code_for(HUB)
    state = place_leaf_state(HUB, code)
    t, d = state.lengths['precision'], state.lengths['degree']
    extension = counting_extension(t, d)
    nodes = [node_element(extension, index) for index in range(1, code.row_count + 1)]
    width = 2 * t
    expected_flat = tuple(
        bit
        for index in range(1, code.row_count + 1)
        for word in (
            (index // 3**position) % 3 for position in range(d)
        )
        for bit in tuple((word >> b) & 1 for b in range(width))
    )
    assert state.port('flat_nodes') == expected_flat
    frames = decode_history(state.port('framed_nodes'))
    assert frames == tuple(
        tuple(
            bit
            for word in ((index // 3**position) % 3 for position in range(d))
            for bit in tuple((word >> b) & 1 for b in range(width))
        )
        for index in range(1, code.row_count + 1)
    )
    inverses = decode_flat(extension, state.port('inverse'), code.row_count)
    negatives = decode_flat(extension, state.port('negative'), code.row_count)
    for node, inverse, negative in zip(nodes, inverses, negatives):
        assert extension.mul(node, inverse) == extension.one
        assert negative == extension.neg(node)
    differences = decode_flat(extension, state.port('differences'), code.row_count**2)
    for left_index, left in enumerate(nodes):
        for right_index, right in enumerate(nodes):
            entry = differences[left_index * code.row_count + right_index]
            if left == right:
                assert entry == extension.zero
            else:
                assert extension.mul(extension.sub(left, right), entry) == extension.one


@pytest.mark.parametrize('cnf,expected', CASES)
def test_row_stream_round_trips(cnf, expected):
    code = _code_for(cnf)
    stream, payloads = serialize_leaf_history(code)
    width = (3 * cnf.m).bit_length()
    assert decode_history(stream) == payloads
    assert history_from_payloads(tuple(p[::-1] for p in payloads)) == stream
    assert len(stream) == 1 + sum(2 * len(payload) + 2 for payload in payloads)
    for index, payload in enumerate(payloads):
        assert len(payload) % width == 0
        assert payload_to_slots(payload, width, cnf.m) == code.row(index).terms


def test_append_encoding_order_and_bits():
    code = _code_for(STAR)
    stream, payloads = serialize_leaf_history(code)
    frame = []
    for bit in payloads[0][::-1]:
        frame.extend((1, bit))
    frame.append(0)
    assert append_payload(stream, payloads[0]) == (1,) + tuple(frame) + stream
    twice = append_payload(append_payload(stream, payloads[0]), payloads[1])
    assert decode_history(twice) == (payloads[1][::-1], payloads[0][::-1]) + payloads


@pytest.mark.parametrize('cnf,expected', CASES)
def test_connected_run_matches_brute_and_the_direct_greedy(cnf, expected):
    result = run_assembled_leaf(cnf)
    direct = count_cnf_greedy(cnf)
    assert result.counted and direct.counted
    assert result.count == direct.count == brute_count(cnf) == expected
    assert result.rows_round_trip and result.score_paths_equal
    assert result.count_word_bits <= cnf.n + 1 and result.port_bits > 0
    assert result.polynomial_terms == 0 and result.same_declared_n


@pytest.mark.parametrize(
    'budgets',
    [
        AssemblyBudgets(max_nodes=1),
        AssemblyBudgets(max_extension_degree=1),
        AssemblyBudgets(max_exponent_unary=6),
        AssemblyBudgets(max_table_entries=8),
        AssemblyBudgets(max_stream_bits=2),
    ],
)
def test_budget_refusals_publish_no_counts(budgets):
    result = run_assembled_leaf(STAR, budgets=budgets)
    assert result.status == 'refused' and result.count is None
    assert result.count_word_bits == 0


def test_node_count_below_three_to_the_d_guard():
    oversized = CNF(16, [[(index, False)] for index in range(16)])
    code = _code_for(oversized)
    assert code.row_count >= 3
    with pytest.raises(ValueError, match='3\\^d'):
        place_leaf_state(oversized, code, degree=1)
    default = place_leaf_state(oversized, code)
    assert default.lengths['degree'] == code.row_count.bit_length() + 1


def test_terminal_leaves_through_the_assembled_run():
    free = run_assembled_leaf(CNF(8, []))
    assert free.count == 256 and free.row_count == 0 and free.port_bits == 0
    empty = run_assembled_leaf(CNF(8, [[]]))
    assert empty.count == 0
    single = run_assembled_leaf(CNF(1, [[(0, False)]]))
    assert single.count == 1


def test_recorded_port_map_and_receipt_flags():
    assert PORT == {
        'ambient': 191, 'modulus': 194, 'table': 201, 'width': 202,
        'scalar_zero': 203, 'degree': 204, 'one': 207, 'exponent': 209,
        'minus': 211, 'precision': 212, 'jet_order': 243, 'half': 253,
        'flat_nodes': 256, 'zero': 280, 'idempotent': 281, 'framed_nodes': 300,
        'inverse': 713, 'negative': 714, 'differences': 715, 'row_width': 716,
        'ring_width': 717,
    }
    result = run_assembled_leaf(STAR)
    assert result.theorem_1_1_complete is False
    assert result.lean_proof_complete is False
    assert result.recursive_contraction is False


def test_malformed_sequences_are_rejected():
    with pytest.raises(ValueError, match='width'):
        payload_to_slots((1, 0, 1), 2, 3)
    with pytest.raises(ValueError):
        decode_history((1, 1))
    code = _code_for(STAR)
    stream, _payloads = serialize_leaf_history(code)
    with pytest.raises(ValueError):
        decode_history(stream[:-1])


def test_product_sanity_on_an_independent_small_case():
    cnf = CNF(3, [[(0, False), (1, False)], [(1, True), (2, False)]])
    result = run_assembled_leaf(cnf)
    assert result.count == brute_count(cnf)
    assert result.rows_round_trip
    extra = CNF(2, [[(0, False), (1, False)], [(0, True), (1, True)]])
    assert run_assembled_leaf(extra).count in (0, 1, 2)
