"""Step 28 independent checks for streamed correction pairs and Taylor readback."""

from itertools import product
from math import comb

import pytest

from counting_sat_algebra.cnf import CNF
from counting_sat_algebra.lean_ext import counting_extension, node_element
from counting_sat_algebra.lean_taylor import (
    NONE, OUT, IN,
    TaylorBudgets,
    TaylorRefusal,
    decode_pair,
    decode_stream,
    difference_coefficient,
    direct_sum,
    encode_pair,
    iter_pairs,
    node_coefficient_stream,
    pair_count,
    reconstruct_from_stream,
    streamed_boolean_sum,
    streamed_taylor_count,
)

STAR = CNF(4, [[(0, False), (v, False)] for v in range(1, 4)])
HUB = CNF(8, [[(0, False), (i, False)] for i in range(1, 8)])


def _mul_terms(left, right):
    output = {}
    for a, x in left.items():
        for b, y in right.items():
            support = a | b
            output[support] = output.get(support, 0) + x * y
    return {support: coefficient for support, coefficient in output.items() if coefficient}


def _model_terms(cnf):
    terms = {frozenset(): 1}
    for clause in cnf.clauses:
        violation = {frozenset(): 1}
        for variable, negated in clause:
            if negated:
                factor = {frozenset((variable,)): 1}
            else:
                factor = {frozenset(): 1, frozenset((variable,)): -1}
            violation = _mul_terms(violation, factor)
        product = _mul_terms(terms, violation)
        keys = set(terms) | set(product)
        terms = {
            s: terms.get(s, 0) - product.get(s, 0) for s in keys
            if terms.get(s, 0) - product.get(s, 0)
        }
    return terms


@pytest.mark.parametrize('variables,order', list(product(range(0, 6), range(1, 4))))
def test_pair_count_matches_the_support_square_majorants(variables, order):
    expected = sum(2**j * comb(variables, j) for j in range(order))
    assert pair_count(variables, order) == expected
    words = list(iter_pairs(variables, order))
    assert len(words) == expected
    assert len(set(words)) == len(words)
    independent = {
        (
            frozenset(i for i, bit in enumerate(e_bits) if bit),
            frozenset(i for i, bit in enumerate(u_bits) if bit),
        )
        for e_bits in product((0, 1), repeat=variables)
        for u_bits in product((0, 1), repeat=variables)
        if all(b <= a for a, b in zip(e_bits, u_bits)) and sum(e_bits) < order
    }
    assert {decode_pair(word) for word in words} == independent


def test_pair_enumeration_order_and_codec_round_trips():
    assert list(iter_pairs(2, 2)) == [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)]
    assert list(iter_pairs(3, 2)) == [
        (NONE, NONE, NONE),
        (NONE, NONE, OUT),
        (NONE, NONE, IN),
        (NONE, OUT, NONE),
        (NONE, IN, NONE),
        (OUT, NONE, NONE),
        (IN, NONE, NONE),
    ]
    for word in iter_pairs(4, 3):
        outer, inner = decode_pair(word)
        assert encode_pair(outer, inner, 4) == word
    with pytest.raises(ValueError, match='subset'):
        encode_pair(frozenset({0}), frozenset({1}), 2)
    with pytest.raises(ValueError, match='coordinate'):
        encode_pair(frozenset({4}), frozenset(), 2)


def test_difference_matches_the_support_table_on_all_masks():
    probe = {frozenset((0, 1)): 5, frozenset((0,)): 2, frozenset(): 7}
    for size in range(4):
        for subset in product((0, 1, 2), repeat=size):
            mask = frozenset(subset)
            assert difference_coefficient(probe, mask) == probe.get(mask, 0)


@pytest.mark.parametrize('variables', [3, 4, 5])
def test_streamed_sum_equals_direct_and_brute(variables):
    rng = __import__('random').Random(2800 + variables)
    order = 3
    terms = {}
    for _ in range(6):
        support = frozenset(rng.sample(range(variables), rng.randrange(0, order)))
        coefficient = rng.randrange(-9, 10)
        if coefficient:
            terms[support] = terms.get(support, 0) + coefficient
    terms = {s: c for s, c in terms.items() if c}
    streamed = streamed_boolean_sum(terms, variables, order)
    assert streamed == direct_sum(terms, variables)
    brute = sum(
        sum(c for s, c in terms.items() if all(bits[v] for v in s))
        for bits in product((0, 1), repeat=variables)
    )
    assert streamed == brute


def test_streamed_refusals_keep_no_partial_sum():
    sample = {(0,): 1, (): 2}
    with pytest.raises(TaylorRefusal, match='order'):
        streamed_boolean_sum(sample, 4, 1)
    with pytest.raises(TaylorRefusal, match='pair budget'):
        streamed_boolean_sum(sample, 4, 2, TaylorBudgets(max_pairs=1))
    with pytest.raises(TaylorRefusal, match='term budget'):
        streamed_boolean_sum(sample, 4, 2, TaylorBudgets(max_terms=1))
    with pytest.raises(TaylorRefusal, match='assignment'):
        streamed_boolean_sum(sample, 8, 2, TaylorBudgets(max_assignments=100))
    with pytest.raises(ValueError):
        pair_count(-1, 2)


def test_node_stream_geometry_and_reconstruction():
    extension = counting_extension(3, 2)
    nodes = [node_element(extension, index) for index in (1, 2)]
    target = [extension.constant(4), extension.one]
    stream, total = node_coefficient_stream(extension, nodes, 3, target)
    assert total == len(stream) == 2 * 3 * 2
    blocks = decode_stream(extension, stream, 2, 3)
    assert len(blocks) == 2 and all(len(block) == 3 for block in blocks)
    recovered = reconstruct_from_stream(extension, nodes, 3, stream)
    assert recovered[0] == extension.constant(4)
    assert extension.decode_count(recovered[0]) == 4
    with pytest.raises(ValueError, match='shorter'):
        decode_stream(extension, stream[:-1], 2, 3)
    with pytest.raises(ValueError, match='trailing'):
        decode_stream(extension, stream + (0,), 2, 3)
    with pytest.raises(ValueError, match='ascending'):
        reconstruct_from_stream(extension, [nodes[1], nodes[0]], 3, stream)
    with pytest.raises(TaylorRefusal, match='stream word'):
        node_coefficient_stream(
            extension, nodes, 3, target, TaylorBudgets(max_stream_words=4)
        )


def test_streamed_taylor_counts_and_readback_range():
    violation = {(0, 1): 1, (0,): -1, (1,): -1, (): 1}
    assert streamed_taylor_count({}, 4, precision=3, degree=2, nodes=2, order=3) == 0
    assert (
        streamed_taylor_count(violation, 4, precision=3, degree=2, nodes=2, order=3) == 4
    )
    assert (
        streamed_taylor_count({frozenset(): 1}, 4, precision=3, degree=2, nodes=2, order=3)
        == 16
    )
    with pytest.raises(TaylorRefusal, match='readback range'):
        streamed_taylor_count(violation, 4, precision=1, degree=2, nodes=2, order=3)
    with pytest.raises(TaylorRefusal, match='at least two coefficients'):
        streamed_taylor_count({frozenset(): 1}, 4, precision=3, degree=2, nodes=1, order=1)


@pytest.mark.parametrize('cnf,expected,precision,order', [(STAR, 9, 3, 7), (HUB, 129, 5, 15)])
def test_cnf_model_counts_through_the_stream(cnf, expected, precision, order):
    terms = _model_terms(cnf)
    streamed = streamed_boolean_sum(terms, cnf.n, order)
    assert streamed == expected
    readback = streamed_taylor_count(
        terms, cnf.n, precision=precision, degree=2, nodes=2, order=order,
        budgets=TaylorBudgets(max_degree=16),
    )
    assert readback == expected


def test_unsatisfiable_indicator_streams_to_zero():
    unsat = CNF(2, [[(0, False)], [(0, True)]])
    terms = _model_terms(unsat)
    assert terms == {}
    assert streamed_boolean_sum(terms, 2, 1) == 0
