"""Step 28 finite witnesses for streamed correction pairs and Taylor readback."""

from __future__ import annotations

from itertools import product

from .protocol import _row, stamp_completion
from ..cnf import CNF, clause_poly_bool
from ..lean_ext import counting_extension, node_element
from ..lean_taylor import (
    NONE,
    OUT,
    IN,
    STEP28_SCOPE,
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
    streamed_boolean_sum,
    streamed_taylor_count,
)

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
SIGNED = {(0, 1): 3, (1, 2): -2, (0,): 5, (3,): 1, (): 4}


def _mul_terms(left: dict, right: dict) -> dict:
    output: dict = {}
    for a, x in left.items():
        for b, y in right.items():
            support = a | b
            output[support] = output.get(support, 0) + x * y
    return {support: coefficient for support, coefficient in output.items() if coefficient}


def _model_terms(cnf: CNF) -> dict:
    """Expand the model indicator ∏_e (1 - q_e) as an integer squarefree polynomial."""
    terms: dict = {frozenset(): 1}
    for clause in cnf.clauses:
        required: dict = {}
        for variable, negated in clause:
            required[variable] = int(negated)
        violation: dict = {frozenset(): 1}
        for variable, bit in sorted(required.items()):
            factor = {frozenset(): 1, frozenset((variable,)): -1}
            if bit == 0:
                factor = {frozenset((variable,)): 1}
                factor[frozenset()] = -1
                factor = {support: -coefficient for support, coefficient in factor.items()}
            violation = _mul_terms(violation, factor)
        product = _mul_terms(terms, violation)
        keys = set(terms) | set(product)
        terms = {
            support: terms.get(support, 0) - product.get(support, 0)
            for support in keys
            if terms.get(support, 0) - product.get(support, 0)
        }
    return terms


def _brute_indicator(cnf: CNF) -> int:
    return sum(
        1
        for bits in product((0, 1), repeat=cnf.n)
        if all(clause_poly_bool(bits, clause) for clause in cnf.clauses)
    )


def run_lean_taylor_suite() -> list[dict]:
    print('=== (AD) Streamed correction pairs and Taylor readback ===', flush=True)
    rows = []

    counts_ok = (
        pair_count(0, 3) == 1
        and pair_count(1, 1) == 1
        and pair_count(2, 2) == 5
        and pair_count(3, 2) == 7
        and pair_count(4, 3) == 33
    )
    enumeration_ok = True
    words_cache = {}
    for variables, order in ((2, 2), (3, 2), (3, 3), (4, 2)):
        words = list(iter_pairs(variables, order))
        words_cache[(variables, order)] = words
        enumeration_ok = enumeration_ok and len(words) == pair_count(variables, order)
        for word in words:
            outer, inner = decode_pair(word)
            enumeration_ok = (
                enumeration_ok
                and len(word) == variables
                and inner <= outer
                and encode_pair(outer, inner, variables) == word
                and all(digit in (NONE, OUT, IN) for digit in word)
            )
    enumeration_ok = (
        enumeration_ok
        and words_cache[(2, 2)] == [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)]
        and words_cache[(3, 2)][1] == (NONE, NONE, OUT)
    )
    rows.append(
        _row(
            'lean-taylor-pair-counts-and-enumeration',
            0,
            0,
            counts_ok and enumeration_ok,
            'CHECKED',
            count=None,
            samples=[pair_count(n, J) for n, J in ((0, 3), (2, 2), (4, 3))],
        )
    )

    probe = {frozenset((0, 1)): 5, frozenset((0,)): 2, frozenset(): 7}
    difference_ok = (
        difference_coefficient(probe, frozenset()) == 7
        and difference_coefficient(probe, frozenset({0})) == 2
        and difference_coefficient(probe, frozenset({1})) == 0
        and difference_coefficient(probe, frozenset({0, 1})) == 5
        and difference_coefficient(probe, frozenset({0, 2})) == 0
        and all(
            difference_coefficient(probe, frozenset(subset)) == probe.get(frozenset(subset), 0)
            for size in range(4)
            for subset in product((0, 1, 2), repeat=size)
        )
    )
    rows.append(
        _row(
            'lean-taylor-difference-extracts-monomial-coefficients',
            0,
            0,
            difference_ok,
            'CHECKED',
            count=None,
        )
    )

    signed_streamed = streamed_boolean_sum(SIGNED, 4, 3)
    brute_signed = sum(
        sum(c for s, c in SIGNED.items() if all(bits[v] for v in s))
        for bits in product((0, 1), repeat=4)
    )
    rows.append(
        _row(
            'lean-taylor-streamed-sum-equals-direct-and-brute',
            4,
            0,
            signed_streamed == direct_sum(SIGNED, 4) == brute_signed == 116,
            'COUNTED',
            count=116,
            pairs=pair_count(4, 3),
        )
    )

    refused = 0
    for call in (
        dict(terms=SIGNED, variables=4, order=2),
        dict(terms=SIGNED, variables=4, order=3, budgets=TaylorBudgets(max_pairs=2)),
        dict(terms=SIGNED, variables=4, order=3, budgets=TaylorBudgets(max_terms=1)),
        dict(terms=SIGNED, variables=4, order=3, budgets=TaylorBudgets(max_assignments=10)),
    ):
        try:
            streamed_boolean_sum(
                call.pop('terms'), call.pop('variables'), call.pop('order'), **call
            )
        except TaylorRefusal:
            refused += 1
    try:
        difference_coefficient(SIGNED, frozenset({0, 1}), TaylorBudgets(max_degree=1))
        refused += 1
    except TaylorRefusal:
        pass
    rows.append(
        _row(
            'lean-taylor-budget-refusals',
            0,
            0,
            refused == 4,
            'REFUSED',
            count=None,
            refusals=refused,
        )
    )

    extension = counting_extension(3, 2)
    nodes = [node_element(extension, index) for index in range(1, 3)]
    stream, total = node_coefficient_stream(
        extension, nodes, 3, [extension.constant(4), extension.one]
    )
    blocks = decode_stream(extension, stream, 2, 3)
    geometry_ok = (
        total == 2 * 3 * 2 == len(stream)
        and len(blocks) == 2
        and all(len(block) == 3 for block in blocks)
        and blocks == reconstruct_stream_blocks(extension, nodes, stream)
        and all(0 <= word < 2 ** (2 * extension.precision) for word in stream)
    )
    errors = 0
    for bad in (stream[:-1], stream + (0,)):
        try:
            decode_stream(extension, bad, 2, 3)
        except ValueError:
            errors += 1
    try:
        encode_pair(frozenset({0}), frozenset({1}), 2)
    except ValueError:
        errors += 1
    rows.append(
        _row(
            'lean-taylor-node-stream-geometry-and-errors',
            0,
            0,
            geometry_ok and errors == 3,
            'CHECKED',
            count=None,
            words=total,
        )
    )

    violation = {(0, 1): 1, (0,): -1, (1,): -1, (): 1}
    empty = streamed_taylor_count({}, 4, precision=3, degree=2, nodes=2, order=3)
    violation_count = streamed_taylor_count(violation, 4, precision=3, degree=2, nodes=2, order=3)
    free = streamed_taylor_count({frozenset(): 1}, 4, precision=3, degree=2, nodes=2, order=3)
    readback_refused = 0
    try:
        streamed_taylor_count(violation, 4, precision=1, degree=2, nodes=2, order=3)
    except TaylorRefusal:
        readback_refused = 1
    rows.append(
        _row(
            'lean-taylor-streamed-counts-and-unique-readback',
            4,
            0,
            empty == 0 and violation_count == 4 and free == 16 and readback_refused == 1,
            'COUNTED',
            count=None,
            violation=violation_count,
            free_bits=free,
        )
    )

    cnf_ok = True
    cnf_details = []
    for cnf, expected, precision, order in (
        (STAR, 9, 3, 7),
        (FOREST, 25, 4, 9),
        (HUB, 129, 5, 15),
    ):
        terms = _model_terms(cnf)
        streamed = streamed_boolean_sum(terms, cnf.n, order)
        readback = streamed_taylor_count(
            terms,
            cnf.n,
            precision=precision,
            degree=2,
            nodes=2,
            order=order,
            budgets=TaylorBudgets(max_degree=16),
        )
        valid = streamed == readback == _brute_indicator(cnf) == expected
        cnf_ok = cnf_ok and valid
        cnf_details.append(dict(n=cnf.n, m=cnf.m, count=readback, pairs=pair_count(cnf.n, order)))
    rows.append(
        _row(
            'lean-taylor-cnf-model-counts-through-the-stream',
            0,
            0,
            cnf_ok,
            'COUNTED',
            count=None,
            cases=cnf_details,
        )
    )

    single = streamed_boolean_sum({frozenset(): 3}, 4, 1)
    single_refused = 0
    try:
        streamed_boolean_sum(violation, 4, 1)
    except TaylorRefusal:
        single_refused = 1
    rows.append(
        _row(
            'lean-taylor-first-order-and-terminal-cases',
            4,
            0,
            single == 3 * 16 and single_refused == 1,
            'COUNTED',
            count=single,
        )
    )

    stream_refused = 0
    try:
        node_coefficient_stream(
            extension,
            nodes,
            3,
            [extension.constant(4), extension.one],
            TaylorBudgets(max_stream_words=4),
        )
    except TaylorRefusal:
        stream_refused = 1
    rows.append(
        _row(
            'lean-taylor-stream-word-budget-refusal',
            0,
            0,
            stream_refused == 1,
            'REFUSED',
            count=None,
        )
    )

    for row in rows:
        stamp_completion(row)
        print(f'  [AD] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP28_SCOPE}', flush=True)
    return rows


def reconstruct_stream_blocks(extension, nodes, stream):
    """A second decoder written from the block description, for the geometry row."""
    blocks = []
    position = 0
    for _ in nodes:
        jet = []
        for _ in range(3):
            jet.append(extension.element(stream[position : position + extension.degree]))
            position += extension.degree
        blocks.append(tuple(jet))
    return blocks
