"""Step 27 finite witnesses for counting extension rings and Hermite placement."""

from __future__ import annotations

from itertools import product
from math import isqrt

from .protocol import _row, stamp_completion
from ..lean_code import radix_words
from ..lean_ext import (
    STEP27_SCOPE,
    CountingExtension,
    ExtensionBudgets,
    ExtensionRefusal,
    counting_extension,
    extension_schedule,
    hermite_reconstruct,
    is_irreducible_f3,
    layout_cutoff,
    layout_den,
    layout_energy,
    layout_q,
    layout_r,
    node_code_words,
    node_element,
    node_stream,
    search_irreducible_polynomial,
    taylor_jets,
)


def _reducible_by_multiplication(coefficients) -> bool:
    """Independent factor oracle: multiply out monic divisors instead of dividing."""
    degree = len(coefficients)
    poly = list(coefficients) + [1]
    for a_degree in range(1, degree):
        for a in product((0, 1, 2), repeat=a_degree):
            divisor = list(a) + [1]
            for b_degree in range(1, degree - a_degree + 1):
                if a_degree + b_degree != degree:
                    continue
                for b in product((0, 1, 2), repeat=b_degree):
                    other = list(b) + [1]
                    product_poly = [0] * (degree + 1)
                    for i, x in enumerate(divisor):
                        for j, y in enumerate(other):
                            product_poly[i + j] = (product_poly[i + j] + x * y) % 3
                    if product_poly == poly:
                        return True
    return False


def _independent_schedule(ambient, clauses, width, factors, occurrence):
    level = (ambient + 16).bit_length()
    w = level.bit_length()
    s = isqrt(level) // w**2
    r = max(2, s)
    den = 64 * r.bit_length()
    q = level // den
    energy = (64 * occurrence * r * level**2).bit_length()
    radius = isqrt(128 * q * energy) + 4
    precision = ambient + clauses.bit_length() + 1
    degree = q.bit_length() + 1 if q else 1
    order = width * radius * (2 * factors + 3 * (precision - 1) + 1) // q + 1 if q else 0
    return level, w, s, r, q, radius, precision, degree, order


def run_lean_ext_suite() -> list[dict]:
    print('=== (AC) Counting extension rings and Hermite placement ===', flush=True)
    rows = []

    samples = (0, 1, 8, 255, 1024, 4096)
    levels_ok = all(
        extension_schedule(n, clauses=7, width=3, factors=2, occurrence=2).level
        == (n + 16).bit_length()
        and extension_schedule(n, clauses=7, width=3, factors=2, occurrence=2).w
        == (n + 16).bit_length().bit_length()
        and layout_r(n)
        == max(2, isqrt((n + 16).bit_length()) // ((n + 16).bit_length().bit_length()) ** 2)
        and layout_den(n) == 64 * layout_r(n).bit_length()
        and layout_q(n) == (n + 16).bit_length() // layout_den(n)
        and layout_energy(n, 2) == (64 * 2 * layout_r(n) * (n + 16).bit_length() ** 2).bit_length()
        and layout_cutoff(n, 2) == isqrt(128 * layout_q(n) * layout_energy(n, 2)) + 4
        for n in samples
    )
    rows.append(
        _row(
            'lean-ext-schedule-early-levels',
            0,
            0,
            levels_ok,
            'CHECKED',
            count=None,
            levels=[(n + 16).bit_length() for n in samples],
            node_counts=[layout_q(n) for n in samples],
        )
    )

    small = extension_schedule(1024, clauses=7, width=3, factors=2, occurrence=2)
    boundary = extension_schedule(1 << (2**24), clauses=7, width=3, factors=2, occurrence=2)
    guarded = extension_schedule(1 << (2**25), clauses=7, width=3, factors=2, occurrence=2)
    exact = guarded.width * guarded.R * (2 * guarded.factors + 3 * (guarded.precision - 1) + 1)
    rows.append(
        _row(
            'lean-ext-schedule-guard-and-hermite-margin',
            0,
            0,
            not small.guarded
            and small.Q == 0
            and not boundary.guarded
            and 4 * boundary.R > boundary.Q
            and guarded.guarded
            and guarded.hermite_conditions_exceed_degree
            and guarded.field_cardinality_cap
            and guarded.Q * guarded.jet_order == exact + (guarded.Q - exact % guarded.Q)
            and guarded.Q * guarded.jet_order > exact,
            'CHECKED',
            count=None,
            small_q=small.Q,
            boundary_q=boundary.Q,
            boundary_radius=boundary.R,
            guarded_q=guarded.Q,
            guarded_radius=guarded.R,
            guarded_precision_bits=guarded.precision.bit_length(),
            guarded_jet_order_bits=guarded.jet_order.bit_length(),
        )
    )

    witnesses = []
    search_ok = True
    for degree, expected in ((1, (0,)), (2, (1, 0)), (3, (1, 2, 0)), (4, None)):
        found, tested, _checks = search_irreducible_polynomial(degree)
        rejected = list(radix_words(3, degree))[: tested - 1]
        search_ok = (
            search_ok
            and (expected is None or found == expected)
            and tested == len(rejected) + 1
            and found not in rejected
            and is_irreducible_f3(found)
            and not _reducible_by_multiplication(found)
            and (degree == 1 or _reducible_by_multiplication((0,) * degree))
            and all(_reducible_by_multiplication(candidate) for candidate in rejected)
        )
        witnesses.append(dict(degree=degree, witness=list(found), tested=tested))
    rows.append(
        _row(
            'lean-ext-irreducible-search-first-witness',
            0,
            0,
            search_ok,
            'CHECKED',
            count=None,
            witnesses=witnesses,
        )
    )

    extension = counting_extension(2, 2)
    x_element = extension.element([0, 1])
    samples_two = [extension.element([2, 1]), extension.element([0, 3]), extension.element([7, 0])]
    axioms_ok = (
        extension.coefficients == (1, 0)
        and extension.mul(x_element, x_element) == extension.constant(8)
        and all(
            extension.mul(a, b) == extension.mul(b, a)
            and extension.add(a, b) == extension.add(b, a)
            for a in samples_two
            for b in samples_two
        )
        and all(
            extension.mul(a, extension.add(b, c))
            == extension.add(extension.mul(a, b), extension.mul(a, c))
            for a in samples_two
            for b in samples_two
            for c in samples_two
        )
        and all(
            all(0 <= coefficient < extension.modulus for coefficient in extension.mul(a, b))
            for a in samples_two
            for b in samples_two
        )
    )
    rows.append(
        _row(
            'lean-ext-extension-ring-arithmetic',
            0,
            0,
            axioms_ok,
            'CHECKED',
            count=None,
            defining_polynomial=list(extension.coefficients),
            x_squared=list(extension.mul(x_element, x_element)),
        )
    )

    inverse_ok = True
    for value in samples_two:
        unit = extension.is_unit(value)
        residue_zero = not any(extension.residue(value))
        inverse_ok = inverse_ok and unit == (not residue_zero)
        if unit:
            inverse_ok = (
                inverse_ok and extension.mul(value, extension.inverse(value)) == extension.one
            )
    refusals = 0
    for bad in (extension.zero, extension.constant(3)):
        try:
            extension.inverse(bad)
        except ValueError:
            refusals += 1
    rows.append(
        _row(
            'lean-ext-residue-units-and-newton-inverse',
            0,
            0,
            inverse_ok and refusals == 2,
            'CHECKED',
            count=None,
            modulus=extension.modulus,
        )
    )

    readback_ok = (
        extension.decode_count(extension.nat_cast(5)) == 5
        and extension.nat_cast(extension.modulus + 5) == extension.nat_cast(5)
        and extension.nat_cast(extension.modulus) == extension.zero
        and extension.nat_cast(4) != extension.nat_cast(5)
        and all(
            extension.nat_cast(a) != extension.nat_cast(b)
            for a, b in ((0, 1), (4, 5), (8, 9 - 1))
            if a != b
        )
        and extension.decode_count(extension.nat_cast(extension.modulus - 1))
        == extension.modulus - 1
    )
    rows.append(
        _row(
            'lean-ext-unique-natural-readback',
            0,
            0,
            readback_ok,
            'CHECKED',
            count=None,
            modulus=extension.modulus,
        )
    )

    nodes = [node_element(extension, index) for index in range(1, 4)]
    node_ok = (
        len(set(nodes)) == 3
        and all(extension.is_unit(extension.sub(a, b)) for a in nodes for b in nodes if a != b)
        and all(
            extension.mul(extension.sub(a, b), extension.inverse(extension.sub(a, b)))
            == extension.one
            for a in nodes
            for b in nodes
            if a != b
        )
        and node_code_words(extension, 3) == (0, 1)
        and node_element(extension, 3) == extension.element([0, 1])
    )
    rejected = 0
    for bad in (0, 9):
        try:
            node_element(extension, bad)
        except ValueError:
            rejected += 1
    rows.append(
        _row(
            'lean-ext-nodes-distinct-unit-differences',
            0,
            0,
            node_ok and rejected == 2,
            'CHECKED',
            count=None,
            nodes=[list(node) for node in nodes],
        )
    )

    c1, c2 = nodes[0], nodes[1]
    square = [extension.zero, extension.zero, extension.one]
    jets_ok = taylor_jets(extension, c1, square, 2) == ((1, 0), (2, 0)) and taylor_jets(
        extension, c2, square, 2
    ) == ((4, 0), (4, 0))
    cubic = [
        extension.element([3, 0]),
        extension.element([1, 2]),
        extension.element([8, 1]),
        extension.element([2, 2]),
    ]
    recovered = hermite_reconstruct(
        extension,
        [c1, c2],
        [taylor_jets(extension, c1, cubic, 2), taylor_jets(extension, c2, cubic, 2)],
        2,
    )
    quartic = cubic + [extension.element([5, 3])]
    mismatch = hermite_reconstruct(
        extension,
        [c1, c2],
        [taylor_jets(extension, c1, quartic, 2), taylor_jets(extension, c2, quartic, 2)],
        2,
    )
    rows.append(
        _row(
            'lean-ext-jets-and-hermite-recovery',
            0,
            0,
            jets_ok
            and recovered == tuple(cubic)
            and len(recovered) <= 2 * 2
            and mismatch != tuple(quartic),
            'CHECKED',
            count=None,
            recovered=[list(coefficient) for coefficient in recovered],
        )
    )

    single = hermite_reconstruct(extension, [c2], [taylor_jets(extension, c2, cubic, 4)], 4)
    rows.append(
        _row(
            'lean-ext-single-node-full-jet-recovery',
            0,
            0,
            single == tuple(cubic),
            'CHECKED',
            count=None,
            order=4,
        )
    )

    stream, bits = node_stream(2, 2, 3)
    stream_ok = (
        bits == 3 * 2 * 2 * 2
        and stream == tuple(node_code_words(extension, index) for index in range(1, 4))
        and all(len(word) == 2 for word in stream)
        and all(all(0 <= entry < 2**4 for entry in word) for word in stream)
    )
    rows.append(
        _row(
            'lean-ext-node-stream-geometry',
            0,
            0,
            stream_ok,
            'CHECKED',
            count=None,
            words=len(stream),
            bits=bits,
        )
    )

    refused = 0
    for call in (
        dict(degree=7),
        dict(degree=2, max_candidates=1),
        dict(degree=6, max_division_checks=3),
        dict(degree=4, max_ring_operations=5),
    ):
        try:
            degree = call.pop('degree')
            search_irreducible_polynomial(degree, budgets=ExtensionBudgets(**call))
        except ExtensionRefusal:
            refused += 1
    rows.append(
        _row(
            'lean-ext-search-refusals-without-witness',
            0,
            0,
            refused == 4,
            'REFUSED',
            count=None,
            refusals=refused,
        )
    )

    terminal = CountingExtension(1, 1, (0,))
    terminal_ok = (
        terminal.mul(terminal.constant(2), terminal.constant(2)) == terminal.constant(1)
        and terminal.decode_count(terminal.nat_cast(2)) == 2
        and node_element(terminal, 2) == terminal.constant(2)
        and search_irreducible_polynomial(1)[0] == (0,)
    )
    rows.append(
        _row(
            'lean-ext-terminal-precision-one',
            0,
            0,
            terminal_ok,
            'CHECKED',
            count=None,
        )
    )

    for row in rows:
        stamp_completion(row)
        print(f'  [AC] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP27_SCOPE}', flush=True)
    return rows
