"""Finite streamed correction-pair evaluation and Hermite counting.

``BooleanDifferences`` decomposes a Boolean sum over correction pairs
``(E, U)`` with ``|E| < J`` and ``U ⊆ E``: the squarefree coefficient on ``E``
is the alternating sum ``difference(E, p) = Σ_{U⊆E} (-1)^(|E|-|U|)·p(U)``,
and for ``J`` above every support size the whole sum is
``Σ_E difference(E, p)·2^(n-|E|)``. ``DifferencePairs`` streams the pairs
through one depth-first fold without retaining a completed-pair list, and the
native leaf emits, per node, exactly the ``J`` Taylor coefficients of the
node-shifted target; the reconstruction reads that stream and returns the
target's value at zero. The finite port streams the same pairs, accumulates
exact integers, emits and re-reads the per-node coefficient stream through
the Step 27 extension ring, and checks every count against a direct support
evaluation. The production native leaf still evaluates the full squarefree
support table; replacing it and the whole-machine resource transfer remain
open.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import combinations
from math import comb

from .lean_ext import (
    CountingExtension,
    counting_extension,
    hermite_reconstruct,
    node_element,
    taylor_jets,
)
from .lean_ring import natural

STEP28_SCOPE = (
    'Streamed correction pairs, the Boolean finite-difference decomposition with '
    'exact accumulation, and the per-node J-coefficient stream re-read by Newton '
    'Hermite. The production native leaf still uses the full squarefree support '
    'table; placement, the whole native machine and the Python bit-time/space '
    'transfer remain open. No Theorem 1.1, no 2^{o(n)} claim.'
)

NONE, OUT, IN = 0, 1, 2


class TaylorRefusal(ValueError):
    """An explicit limit stops the streamed evaluation before completion."""


@dataclass(frozen=True)
class TaylorBudgets:
    max_pairs: int = 2_000_000
    max_assignments: int = 20_000_000
    max_stream_words: int = 100_000
    max_degree: int = 8
    max_terms: int = 100_000

    def __post_init__(self) -> None:
        for item in fields(self):
            object.__setattr__(self, item.name, natural(getattr(self, item.name), item.name))


def pair_count(variables: int, order: int) -> int:
    """``DifferencePairs.pairCount``: the exact number of correction pairs."""
    variables, order = natural(variables, 'variables'), natural(order, 'order')
    return sum(2**j * comb(variables, j) for j in range(order))


def assignment_count(variables: int, order: int) -> int:
    """The E-conditioned assignments visited by the whole stream."""
    variables, order = natural(variables, 'variables'), natural(order, 'order')
    return sum((1 << variables) * comb(variables, j) for j in range(order))


def pair_walk(variables: int, order: int):
    """The unchecked ``DifferencePairs.fold`` stream; callers own the budgets."""
    variables, order = natural(variables, 'variables'), natural(order, 'order')
    if not order:
        return

    def walk(coordinate: int, budget: int, prefix: tuple[int, ...]):
        if budget == 0:
            return
        if coordinate == variables:
            yield prefix
            return
        yield from walk(coordinate + 1, budget, prefix + (NONE,))
        yield from walk(coordinate + 1, budget - 1, prefix + (OUT,))
        yield from walk(coordinate + 1, budget - 1, prefix + (IN,))

    yield from walk(0, order, ())


def iter_pairs(variables: int, order: int, budgets: TaylorBudgets | None = None):
    """``DifferencePairs.fold`` order: coordinate by coordinate, none then out then in."""
    variables, order = natural(variables, 'variables'), natural(order, 'order')
    limits = budgets or TaylorBudgets()
    if pair_count(variables, order) > limits.max_pairs:
        raise TaylorRefusal('correction pair budget exceeded')
    if assignment_count(variables, order) > limits.max_assignments:
        raise TaylorRefusal('pair assignment budget exceeded')
    yield from pair_walk(variables, order)


def decode_pair(word: tuple[int, ...]) -> tuple[frozenset[int], frozenset[int]]:
    """``encodePair`` decoded: coordinate membership in E and in U ⊆ E."""
    outer = frozenset(i for i, digit in enumerate(word) if digit != NONE)
    inner = frozenset(i for i, digit in enumerate(word) if digit == IN)
    return outer, inner


def encode_pair(outer: frozenset[int], inner: frozenset[int], variables: int) -> tuple[int, ...]:
    if not inner <= outer:
        raise ValueError('the inner mask must be a subset of the outer mask')
    if any(v < 0 or v >= variables for v in outer):
        raise ValueError('pair coordinates exceed the variable count')
    return tuple(
        NONE if i not in outer else (IN if i in inner else OUT) for i in range(variables)
    )


def _check_terms(terms, budgets: TaylorBudgets) -> tuple[dict[frozenset[int], int], int]:
    if len(terms) > budgets.max_terms:
        raise TaylorRefusal('squarefree term budget exceeded')
    cleaned = {}
    for support, coefficient in terms.items():
        support = frozenset(natural(v, 'variable') for v in support)
        coefficient = int(coefficient)
        if coefficient:
            cleaned[support] = coefficient
    return cleaned, max((len(s) for s in cleaned), default=0)


def evaluate_terms(terms, mask: frozenset[int]) -> int:
    """Evaluate a squarefree polynomial at the Boolean point ``mask``."""
    return sum(coefficient for support, coefficient in terms.items() if support <= mask)


def direct_sum(terms, variables: int) -> int:
    """The multilinear weight identity ``Σ_S coefficient(S)·2^(n-|S|)``."""
    return sum(coefficient * 2 ** (variables - len(support)) for support, coefficient in terms.items())


def difference_coefficient(terms, outer, budgets: TaylorBudgets | None = None) -> int:
    """``difference``: the coefficient of the squarefree monomial on ``outer``."""
    limits = budgets or TaylorBudgets()
    cleaned, degree = _check_terms(terms, limits)
    if len(outer) > limits.max_degree:
        raise TaylorRefusal('difference degree budget exceeded')
    if len(outer) > degree:
        return 0
    signs = 0
    order = sorted(outer)
    for size in range(len(order) + 1):
        for subset in combinations(order, size):
            mask = frozenset(subset)
            sign = -1 if (len(order) - size) % 2 else 1
            signs += sign * evaluate_terms(cleaned, mask)
    return signs


def streamed_boolean_sum(
    terms, variables: int, order: int, budgets: TaylorBudgets | None = None
) -> int:
    """The pair-streamed sum, equal to ``direct_sum`` when ``order`` exceeds the degree."""
    limits = budgets or TaylorBudgets()
    cleaned, degree = _check_terms(terms, limits)
    if order <= degree:
        raise TaylorRefusal('the Taylor order must exceed every support size')
    total = 0
    for word in iter_pairs(variables, order, limits):
        outer, inner = decode_pair(word)
        sign = -1 if (len(outer) - len(inner)) % 2 else 1
        total += sign * (1 << (variables - len(outer))) * evaluate_terms(cleaned, inner)
    return total


def node_coefficient_stream(
    extension: CountingExtension,
    nodes,
    order: int,
    target,
    budgets: TaylorBudgets | None = None,
) -> tuple[tuple[tuple[int, ...], ...], int]:
    """Per-node blocks of the ``J`` Taylor coefficients of the shifted target.

    Each coefficient is a ring element stored as ``d`` little-endian words;
    the flat stream is exactly ``Q * J * d`` words, including zero words.
    """
    limits = budgets or TaylorBudgets()
    order = natural(order, 'order')
    if order > limits.max_degree:
        raise TaylorRefusal('jet order budget exceeded')
    blocks = []
    for index, node in enumerate(nodes, start=1):
        coefficient_check = node_element(extension, index)
        if coefficient_check != node:
            raise ValueError('the nodes must be the radix nodes in ascending order')
        blocks.append(taylor_jets(extension, node, target, order))
    words = tuple(
        int(coefficient)
        for block in blocks
        for element in block
        for coefficient in element
    )
    total_words = len(words)
    if total_words > limits.max_stream_words:
        raise TaylorRefusal('stream word budget exceeded')
    return words, total_words


def decode_stream(
    extension: CountingExtension, stream_words, nodes: int, order: int
):
    """Re-read the flat stream into per-node jet blocks, without any other input."""
    nodes, order = natural(nodes, 'nodes'), natural(order, 'order')
    blocks = []
    position = 0
    for index in range(nodes):
        jet = []
        for _ in range(order):
            if position + extension.degree > len(stream_words):
                raise ValueError('the stream is shorter than Q*J*d words')
            jet.append(extension.element(stream_words[position: position + extension.degree]))
            position += extension.degree
        blocks.append(tuple(jet))
    if position != len(stream_words):
        raise ValueError('the stream retains trailing words')
    return blocks


def reconstruct_from_stream(
    extension: CountingExtension, nodes, order: int, stream_words
) -> tuple[tuple[int, ...], ...]:
    """Newton-Hermite from the re-read stream; the count is the value at zero."""
    node_list = [node_element(extension, index) for index in range(1, len(nodes) + 1)]
    if list(nodes) != node_list:
        raise ValueError('the nodes must be the radix nodes in ascending order')
    blocks = decode_stream(extension, stream_words, len(nodes), order)
    return hermite_reconstruct(extension, node_list, blocks, order)


def streamed_taylor_count(
    terms,
    variables: int,
    *,
    precision: int,
    degree: int,
    nodes: int,
    order: int,
    budgets: TaylorBudgets | None = None,
) -> int:
    """Count through the pair stream, the node stream and Hermite readback.

    The target is the count-carrying univariate polynomial ``count + X`` scaled
    by the ring: its value at zero is the count, its shifted jets are the node
    stream, and ``nodes*order > 1`` makes the readback unique. The pair-streamed
    sum supplies both the coefficient stream check and the count itself.
    """
    limits = budgets or TaylorBudgets()
    cleaned, support_degree = _check_terms(terms, limits)
    if order <= support_degree:
        raise TaylorRefusal('the Taylor order must exceed every support size')
    if order < 2:
        raise TaylorRefusal('the readback needs at least two coefficients')
    count = streamed_boolean_sum(cleaned, variables, order, limits)
    if count < 0 or count >= 3**precision:
        raise TaylorRefusal('the streamed sum exceeds the unique readback range')
    extension = counting_extension(precision, degree)
    if nodes * order <= 1:
        raise TaylorRefusal('the condition count must exceed the target degree')
    target = [extension.constant(count), extension.one]
    stream_words, _total = node_coefficient_stream(
        extension,
        [node_element(extension, index) for index in range(1, nodes + 1)],
        order,
        target,
        limits,
    )
    recovered = reconstruct_from_stream(
        extension, [node_element(extension, index) for index in range(1, nodes + 1)],
        order, stream_words,
    )
    value = recovered[0] if recovered else extension.zero
    return extension.decode_count(value)
