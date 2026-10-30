"""Step 27 independent checks for counting extension rings and Hermite placement."""

from itertools import product
from math import isqrt

import pytest

from counting_sat_algebra.lean_code import radix_words
from counting_sat_algebra.lean_ext import (
    CountingExtension,
    ExtensionBudgets,
    ExtensionRefusal,
    counting_extension,
    extension_schedule,
    hermite_reconstruct,
    is_irreducible_f3,
    layout_cutoff,
    layout_energy,
    layout_q,
    layout_r,
    node_code_words,
    node_element,
    node_stream,
    search_irreducible_polynomial,
    taylor_jets,
)


@pytest.mark.parametrize('ambient', [0, 1, 16, 255, 1024, 1 << 20])
def test_schedule_levels_match_the_asymptotic_integers(ambient):
    schedule = extension_schedule(
        ambient, clauses=7, width=3, factors=2, occurrence=2
    )
    level = (ambient + 16).bit_length()
    w = level.bit_length()
    s = isqrt(level) // w**2
    r = max(2, s)
    assert schedule.level == level and schedule.w == w and schedule.s == s
    assert schedule.r == r == layout_r(ambient)
    assert layout_q(ambient) == level // (64 * r.bit_length())
    assert layout_energy(ambient, 2) == (64 * 2 * r * level**2).bit_length()
    assert layout_cutoff(ambient, 2) == isqrt(128 * layout_q(ambient) * layout_energy(ambient, 2)) + 4
    assert schedule.precision == ambient + 7 .bit_length() + 1


def test_schedule_guard_fails_below_and_holds_above_the_radius_threshold():
    boundary = extension_schedule(1 << (2**24), clauses=7, width=3, factors=2, occurrence=2)
    guarded = extension_schedule(1 << (2**25), clauses=7, width=3, factors=2, occurrence=2)
    assert not boundary.guarded and 4 * boundary.R > boundary.Q
    assert guarded.guarded
    exact = guarded.width * guarded.R * (2 * guarded.factors + 3 * (guarded.precision - 1) + 1)
    assert guarded.Q * guarded.jet_order == exact + (guarded.Q - exact % guarded.Q)
    assert guarded.Q * guarded.jet_order > exact
    assert guarded.field_cardinality_cap
    tiny = extension_schedule(1024, clauses=3, width=3, factors=1, occurrence=1)
    assert tiny.Q == 0 and not tiny.guarded and tiny.jet_order == 0


def test_irreducible_search_matches_an_independent_multiplication_oracle():
    def reducible(coefficients):
        degree = len(coefficients)
        poly = list(coefficients) + [1]
        for a_degree in range(1, degree):
            for a in product((0, 1, 2), repeat=a_degree):
                for b in product((0, 1, 2), repeat=degree - a_degree):
                    first, second = list(a) + [1], list(b) + [1]
                    product_poly = [0] * (degree + 1)
                    for i, x in enumerate(first):
                        for j, y in enumerate(second):
                            product_poly[i + j] = (product_poly[i + j] + x * y) % 3
                    if product_poly == poly:
                        return True
        return False

    for degree in (1, 2, 3, 4):
        witness, tested, _checks = search_irreducible_polynomial(degree)
        assert witness is not None and not reducible(witness)
        earlier = list(radix_words(3, degree))[: tested - 1]
        assert len(earlier) == tested - 1
        assert all(reducible(candidate) for candidate in earlier)
        assert is_irreducible_f3(witness)


def test_extension_arithmetic_and_cardinality_are_consistent():
    for precision, degree in ((1, 1), (2, 2), (3, 2), (2, 3)):
        extension = counting_extension(precision, degree)
        modulus = 3**precision
        elements = [extension.element(coefficients) for coefficients in product(range(3), repeat=degree)]
        assert len(set(elements)) == 3**degree
        for left in elements[:8]:
            for right in elements[:8]:
                product_element = extension.mul(left, right)
                assert all(0 <= coefficient < modulus for coefficient in product_element)
                assert extension.mul(left, right) == extension.mul(right, left)
                assert extension.add(left, right) == extension.add(right, left)
        assert extension.mul(extension.one, elements[1]) == elements[1]


@pytest.mark.parametrize('precision,degree', [(1, 1), (2, 2), (3, 2), (2, 4)])
def test_residue_units_and_newton_inverse_round_trip(precision, degree):
    extension = counting_extension(precision, degree)
    checked = 0
    for coefficients in product(range(3**precision), repeat=degree):
        value = extension.element(coefficients)
        if extension.is_unit(value):
            assert extension.mul(value, extension.inverse(value)) == extension.one
            checked += 1
        else:
            assert not any(extension.residue(value))
            with pytest.raises(ValueError):
                extension.inverse(value)
    assert checked >= 1


def test_unique_natural_readback_below_the_modulus():
    extension = counting_extension(2, 2)
    modulus = extension.modulus
    for value in range(0, modulus + 5):
        assert extension.decode_count(extension.nat_cast(value)) == value % modulus
        assert extension.nat_cast(value) == extension.nat_cast(value % modulus)
    assert extension.nat_cast(4) != extension.nat_cast(5)
    assert extension.nat_cast(modulus) == extension.zero


def test_nodes_are_distinct_with_unit_differences_and_physical_words():
    extension = counting_extension(2, 2)
    nodes = [node_element(extension, index) for index in range(1, 3**extension.degree)]
    assert len(set(nodes)) == len(nodes)
    for left in nodes:
        for right in nodes:
            if left != right:
                difference = extension.sub(left, right)
                assert extension.is_unit(difference)
                assert extension.mul(difference, extension.inverse(difference)) == extension.one
    assert node_code_words(extension, 3) == (0, 1)
    stream, bits = node_stream(2, 2, 4)
    assert bits == 4 * 2 * 2 * 2
    assert stream[2] == node_code_words(extension, 3)
    with pytest.raises(ValueError):
        node_element(extension, 0)
    with pytest.raises(ValueError):
        node_element(extension, 9)


@pytest.mark.parametrize('order', [1, 2, 3, 4])
def test_taylor_jets_match_closed_forms(order):
    extension = counting_extension(2, 2)
    node = node_element(extension, 5)
    square = [extension.zero, extension.zero, extension.one]
    jets = taylor_jets(extension, node, square, order)
    expected_value = extension.mul(node, node)
    expected_slope = extension.scale(2, node)
    assert jets[0] == expected_value
    if order >= 2:
        assert jets[1] == expected_slope
    assert jets == taylor_jets(extension, node, square + [extension.zero], order)


@pytest.mark.parametrize('precision,degree,nodes,order', [(2, 2, 2, 2), (2, 2, 1, 4), (1, 1, 2, 1), (3, 2, 2, 2)])
def test_hermite_recovery_is_exact_below_the_condition_count(precision, degree, nodes, order):
    extension = counting_extension(precision, degree)
    node_list = [node_element(extension, i) for i in range(1, nodes + 1)]
    modulus = extension.modulus
    target = [
        extension.element([(7 * j + 3) % modulus, (5 * j + 1) % 3, j % 3, (j + 2) % 3][:degree])
        for j in range(nodes * order)
    ]
    while target and not any(target[-1]):
        target.pop()
    assert len(target) <= nodes * order
    if len(target) < nodes * order:
        target.append(extension.element([1] + [0] * (degree - 1)))
    assert len(target) == nodes * order
    jets = [taylor_jets(extension, node, target, order) for node in node_list]
    recovered = hermite_reconstruct(extension, node_list, jets, order)
    assert recovered == tuple(target)


def test_hermite_mismatches_at_or_above_the_condition_count():
    extension = counting_extension(2, 2)
    node_list = [node_element(extension, i) for i in (1, 2)]
    target = [extension.element([3, 0]), extension.element([1, 2]), extension.element([8, 1]), extension.element([2, 2])]
    jets = [taylor_jets(extension, node, target, 2) for node in node_list]
    assert hermite_reconstruct(extension, node_list, jets, 2) == tuple(target)
    longer = target + [extension.element([5, 3])]
    jets_long = [taylor_jets(extension, node, longer, 2) for node in node_list]
    assert hermite_reconstruct(extension, node_list, jets_long, 2) != tuple(longer)


def test_search_refusals_do_not_return_witnesses():
    with pytest.raises(ExtensionRefusal):
        search_irreducible_polynomial(7)
    with pytest.raises(ExtensionRefusal):
        search_irreducible_polynomial(2, budgets=ExtensionBudgets(max_candidates=1))
    with pytest.raises(ExtensionRefusal):
        search_irreducible_polynomial(6, budgets=ExtensionBudgets(max_division_checks=3))
    with pytest.raises(ValueError):
        search_irreducible_polynomial(0)
    with pytest.raises(ValueError):
        CountingExtension(0, 1, (0,))


def test_precision_one_extension_is_integers_modulo_three():
    extension = counting_extension(1, 1)
    assert extension.modulus == 3
    assert extension.mul(extension.constant(2), extension.constant(2)) == extension.constant(1)
    assert extension.add(extension.constant(2), extension.constant(2)) == extension.constant(1)
    assert extension.decode_count(extension.nat_cast(2)) == 2
