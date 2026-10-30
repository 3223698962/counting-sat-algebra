"""Executable syntax from SubexponentialSAT/Encoding.lean and CountModel.lean.

The binary alphabet uses unary number fields, escaped list items, and a
little-endian count word. Canonicalization removes multiplicities and sorts
both finite-set levels; it does not perform semantic CNF simplification.
"""

from __future__ import annotations

import operator
from collections.abc import Sequence

from .cnf import CNF, Literal

DEFAULT_ENCODING_BITS = 16_000_000


def _natural(value: int, name: str) -> int:
    value = operator.index(value)
    if value < 0:
        raise ValueError(f'{name} must be nonnegative')
    return value


def literal_code(literal: Literal) -> int:
    variable, negated = literal
    variable = _natural(variable, 'variable')
    if not isinstance(negated, bool):
        raise ValueError('literal polarity must be a bool')
    return 2 * variable + int(not negated)


def cnf_numbers(cnf: CNF, *, canonical: bool = True) -> tuple[tuple[int, ...], ...]:
    """Use Lean's polarity convention and finite-set lexicographic order."""
    n = _natural(cnf.n, 'n')
    rows = []
    for clause in cnf.clauses:
        codes = []
        for literal in clause:
            code = literal_code(literal)
            if code >= 2 * n:
                raise ValueError('literal variable exceeds the declared n')
            codes.append(code)
        rows.append(tuple(sorted(set(codes))) if canonical else tuple(codes))
    return tuple(sorted(set(rows))) if canonical else tuple(rows)


def clause_bit_length(clause: Sequence[Literal]) -> int:
    """Contribution of one framed clause to the outer formula list."""
    return 4 + sum(4 * literal_code(literal) + 8 for literal in clause)


def formula_bit_length(cnf: CNF, *, canonical: bool = True) -> int:
    n = _natural(cnf.n, 'n')
    if canonical:
        rows = cnf_numbers(cnf)
        return n + 2 + sum(4 + sum(4 * code + 8 for code in row) for row in rows)
    total = n + 2
    for clause in cnf.clauses:
        total += 4
        for literal in clause:
            code = literal_code(literal)
            if code >= 2 * n:
                raise ValueError('literal variable exceeds the declared n')
            total += 4 * code + 8
    return total


def _frame(bits: bytes) -> bytes:
    framed = bytearray()
    for bit in bits:
        framed.extend((1, bit))
    framed.append(0)
    return bytes(framed)


def _encode_list(items: Sequence[bytes]) -> bytes:
    result = bytearray()
    for item in items:
        result.append(1)
        result.extend(_frame(item))
    result.append(0)
    return bytes(result)


def encode_cnf(
    cnf: CNF, *, canonical: bool = True, max_bits: int = DEFAULT_ENCODING_BITS,
) -> bytes:
    """Return one byte per logical bit, checking the size before bit allocation."""
    limit = _natural(max_bits, 'max_bits')
    size = formula_bit_length(cnf, canonical=canonical)
    if size > limit:
        raise ValueError(f'formula encoding needs {size} bits, above budget {limit}')
    rows = cnf_numbers(cnf, canonical=canonical)
    payload = _encode_list(tuple(
        _encode_list(tuple(b'\x01' * code + b'\x00' for code in row)) for row in rows
    ))
    result = b'\x01' * cnf.n + b'\x00' + payload
    if len(result) != size:
        raise RuntimeError('formula encoding length disagrees with the exact size')
    return result


def _read_nat(bits: bytes, start: int) -> tuple[int, int]:
    end = start
    while end < len(bits) and bits[end] == 1:
        end += 1
    if end == len(bits):
        raise ValueError('unterminated unary field')
    return end - start, end + 1


def _read_frame(bits: bytes, start: int) -> tuple[bytes, int]:
    payload = bytearray()
    while start < len(bits):
        if bits[start] == 0:
            return bytes(payload), start + 1
        if start + 1 == len(bits):
            raise ValueError('incomplete escape frame')
        payload.append(bits[start + 1])
        start += 2
    raise ValueError('unterminated escape frame')


def _read_list(bits: bytes) -> tuple[bytes, ...]:
    items = []
    position = 0
    while position < len(bits) and bits[position] == 1:
        item, position = _read_frame(bits, position + 1)
        items.append(item)
    if position != len(bits) - 1 or bits[position] != 0:
        raise ValueError('unterminated list or trailing bits')
    return tuple(items)


def decode_cnf(bits: bytes, *, max_bits: int = DEFAULT_ENCODING_BITS) -> CNF:
    """Decode the explicit syntax; malformed frames and suffixes are rejected."""
    limit = _natural(max_bits, 'max_bits')
    if len(bits) > limit:
        raise ValueError('formula encoding exceeds the bit budget')
    if any(bit not in (0, 1) for bit in bits):
        raise ValueError('formula encoding must use the binary alphabet')
    n, position = _read_nat(bits, 0)
    clauses = []
    for row in _read_list(bits[position:]):
        literals = []
        for field in _read_list(row):
            code, consumed = _read_nat(field, 0)
            if consumed != len(field) or code >= 2 * n:
                raise ValueError('invalid literal field or variable index')
            literals.append((code // 2, code % 2 == 0))
        clauses.append(literals)
    return CNF(n, clauses)


def count_word(count: int) -> bytes:
    count = _natural(count, 'count')
    return bytes((count >> bit) & 1 for bit in range(max(1, count.bit_length())))


def decode_count_word(word: bytes) -> int:
    if not word or any(bit not in (0, 1) for bit in word):
        raise ValueError('a count word must be a nonempty binary word')
    if len(word) > 1 and word[-1] != 1:
        raise ValueError('a count word must have no redundant high zero')
    return sum(bit << index for index, bit in enumerate(word))
