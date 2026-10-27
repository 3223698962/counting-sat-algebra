"""Independent finite checks of the Appendix D modular primitives."""

import json
import random
from dataclasses import replace
from fractions import Fraction
from itertools import product
from math import comb, prod

import pytest

from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.modular_ring import (
    ModularBudgetRefusal, ModularRing, evaluate_polynomial, extension_ring_for_size, find_irreducible_polynomial,
    interpolate_polynomial, is_irreducible_mod2,
)
from counting_sat_algebra.modular_network import (
    BooleanTensorNetwork, ModularBudgets, TensorVertex, contract_network,
    count_cnf_modular, encode_cnf_network,
)
from counting_sat_algebra.modular_recovery import (
    PacketLayout, TaperedJet, cauchy_local_polynomial, embed_network,
    evaluate_deformation, plan_tapered_recovery, recover_network_modular,
    recover_tapered_jets,
)


def _parallel_network(*, nonlinear=False, precision=5, ell=2):
    ring = ModularRing(precision, ell)

    def row(bits):
        constant = 1 + (bits[0] * bits[1] if nonlinear else bits[0])
        return ring.from_base((constant, 2 - bits[1])[:ell])

    def column(bits):
        return ring.from_base((1 + bits[1], -bits[0])[:ell])

    return BooleanTensorNetwork(ring, 2, (
        TensorVertex('row', 0, (0, 1), ring, row),
        TensorVertex('column', 1, (0, 1), ring, column),
    ))


def _poly_product(left, right):
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] += x * y
    return out


def _rational_cauchy(a, row, column):
    L = 1 << a
    polynomial = [Fraction(1, prod(-q for q in range(-L, L) if q))]
    for q in range(-L, L):
        if q != column - row:
            polynomial = _poly_product(polynomial, (-q, 1))
    return polynomial


def _exact_shift(coefficients, t, alpha, length):
    return tuple(sum((coefficients[j] * comb(j, k) * t ** (j - k) * alpha ** k
                      for j in range(k, len(coefficients))), Fraction(0))
                 for k in range(length))


def _nonlinear_deformation_polynomials():
    # Sum over Y first: A interpolates to (1+u0*u1)+(2-u1)V.
    # The creation-time representative of -X0 modulo 32 is 31*X0.
    a, b = _rational_cauchy(1, 0, 0), _rational_cauchy(1, 0, 1)
    c, d = _rational_cauchy(1, 1, 0), _rational_cauchy(1, 1, 1)
    constant, v_coefficient = [Fraction(6)] + [Fraction(0)] * 6, [Fraction(74)] + [Fraction(0)] * 6
    for weight, left, right in ((3, a, c), (4, b, d), (2, a, d), (2, b, c)):
        for j, coefficient in enumerate(_poly_product(left, right)):
            constant[j] += weight * coefficient
    for weight, left, right in ((62, a, c), (31, b, d), (31, a, d), (31, b, c)):
        for j, coefficient in enumerate(_poly_product(left, right)):
            v_coefficient[j] += weight * coefficient
    for weight, polynomial in ((-3, c), (-4, d)):
        for j, coefficient in enumerate(polynomial):
            v_coefficient[j] += weight * coefficient
    return constant, v_coefficient


def _exact_samples(ring, plan, polynomials):
    samples = []
    for node in plan.nodes:
        local = ring.at_precision(node.precision)
        jets = [_exact_shift(p, node.t, 1 << (plan.a + 1), node.multiplicity) for p in polynomials]
        coefficients = []
        for j in range(node.multiplicity):
            base = [local.from_fraction(jet[j])[0] for jet in jets]
            coefficients.append(local.from_base(base))
        samples.append(TaperedJet(node.t, local, tuple(coefficients)))
    return samples


def test_irreducible_search_matches_independent_trial_division():
    def divides(polynomial, divisor):
        while polynomial.bit_length() >= divisor.bit_length():
            polynomial ^= divisor << (polynomial.bit_length() - divisor.bit_length())
        return polynomial == 0

    for degree in range(1, 8):
        for lower in range(1 << degree):
            candidate = (1 << degree) | lower
            irreducible = not any(divides(candidate, divisor)
                                  for d in range(1, degree // 2 + 1)
                                  for divisor in range(1 << d, 1 << (d + 1)))
            assert is_irreducible_mod2(candidate) == irreducible
        assert is_irreducible_mod2(find_irreducible_polynomial(degree))
    assert find_irreducible_polynomial(4) == 0b10011
    with pytest.raises(ModularBudgetRefusal, match='search not started'):
        find_irreducible_polynomial(1000)
    with pytest.raises(ModularBudgetRefusal, match='candidate budget'):
        find_irreducible_polynomial(4, max_candidates=1)
    with pytest.raises(ValueError, match='irreducible'):
        ModularRing(5, polynomial=0b10101)


def test_integer_extension_reduction_and_truncated_v_products():
    ring = ModularRing(5, 3, 0b111)  # U^2 = -U-1, over Z/32 rather than F_2.
    U, V = ring.element((0, 1)), ring.element((0, 0, 1))
    assert ring.mul(U, U) == ring.element((31, 31))
    assert ring.pow(U, 3) == ring.one
    assert ring.pow(V, 3) == ring.zero
    assert ring.pow(ring.add(ring.one, V), 2) == ring.from_base((1, 2, 1))
    assert ring.from_fraction(Fraction(-2, 3)) == ring.scalar(10)
    with pytest.raises(ValueError, match='odd denominator'):
        ring.from_fraction(Fraction(1, 2))
    with pytest.raises(ValueError, match='nonunit'):
        ring.inverse(ring.add(ring.scalar(2), V))
    with pytest.raises(ValueError, match='coordinates'):
        ring.mul((1,), ring.one)
    with pytest.raises(ValueError, match='base coefficient ring'):
        ring.base_coefficients(U)


def test_reference_extension_degree_and_preflight_limits():
    assert extension_ring_for_size(2, 5, 2) == ModularRing(5, 2, 0b10011)
    assert extension_ring_for_size(3, 7).extension_degree == 8
    with pytest.raises(ModularBudgetRefusal, match='search not started'):
        extension_ring_for_size(1 << 100, 5)
    with pytest.raises(ModularBudgetRefusal, match='search not started'):
        extension_ring_for_size(2, 5000)
    with pytest.raises(ValueError, match='at least two'):
        extension_ring_for_size(1, 5)


@pytest.mark.parametrize('precision', [1, 2, 9])
@pytest.mark.parametrize('ell', [1, 3])
def test_newton_unit_inversion_and_ring_laws(precision, ell):
    ring = ModularRing(precision, ell, 0b10011)
    rng = random.Random(100 + precision + ell)
    for _ in range(10):
        x, y, z = [ring.element(tuple(rng.randrange(ring.modulus) for _ in range(ring.size))) for _ in range(3)]
        if not any(value & 1 for value in x[:4]):
            x = ring.add(x, ring.one)
        assert ring.mul(x, ring.inverse(x)) == ring.one
        assert ring.mul(x, y) == ring.mul(y, x)
        assert ring.mul(ring.mul(x, y), z) == ring.mul(x, ring.mul(y, z))
        assert ring.mul(x, ring.add(y, z)) == ring.add(ring.mul(x, y), ring.mul(x, z))


def test_residue_field_interpolation_retains_all_u_and_v_coordinates():
    ring = ModularRing(7, 3, 0b10011)
    points = tuple(ring.residue_point(index) for index in range(16))
    for j, point in enumerate(points):
        for earlier in points[:j]:
            difference = ring.sub(point, earlier)
            assert ring.mul(difference, ring.inverse(difference)) == ring.one
    coefficients = tuple(ring.element(tuple((j + 3) * (i + 7) for i in range(ring.size))) for j in range(7))
    values = tuple(evaluate_polynomial(ring, coefficients, point) for point in points[:7])
    assert interpolate_polynomial(ring, points[:7], values) == coefficients
    with pytest.raises(ValueError, match='nonunit'):
        interpolate_polynomial(ring, (ring.scalar(0), ring.scalar(2)), (ring.one, ring.one))


def test_canonical_lifting_preserves_creation_precision_not_ring_operations():
    source, target = ModularRing(2), ModularRing(8)
    value = source.scalar(3)
    assert source.canonical_lift(source.add(value, value), target) == (2,)
    assert target.add(source.canonical_lift(value, target), source.canonical_lift(value, target)) == (6,)
    network = BooleanTensorNetwork(target, 1, (
        TensorVertex('a', 0, (0,), source, lambda bits: (-1,)),
        TensorVertex('b', 1, (0,), source, lambda bits: (6,)),
    ))
    # Each edge bit contributes 3*2, using the two separate canonical lifts.
    assert contract_network(network).value == (12,)
    assert contract_network(network.at_precision(12)).value == (12,)
    embedded = embed_network(network.at_precision(12), 0b10011)
    assert embedded.ring.base_coefficients(contract_network(embedded).value) == (12,)
    assert all(vertex.creation_ring.precision == 2 for vertex in embedded.vertices)
    with pytest.raises(ValueError, match='lower precision'):
        target.canonical_lift((1,), source)


def test_parallel_ports_isolated_scalars_and_empty_networks():
    network = _parallel_network()
    result = contract_network(network)
    assert result.value == (9, 4) and result.boolean_vectors_scanned == 4
    scalar = TensorVertex('scalar', 0, (), network.ring, lambda bits: (2, 1))
    with_scalar = replace(network, vertices=network.vertices + (scalar,))
    assert contract_network(with_scalar).value == (18, 17)
    empty = BooleanTensorNetwork(network.ring, 0, ())
    assert contract_network(empty).value == network.ring.one
    assert contract_network(empty).tensor_evaluations == 0
    zero = TensorVertex('zero', 1, (), network.ring, lambda bits: (0,))
    assert contract_network(replace(empty, vertices=(zero,))).value == network.ring.zero
    with pytest.raises(ValueError, match='distinct ports'):
        replace(network.vertices[0], ports=(0, 0))
    with pytest.raises(ValueError, match='each side'):
        replace(network, vertices=(network.vertices[0], replace(network.vertices[1], side=0)))
    with pytest.raises(ValueError, match='each side'):
        replace(network, edges=3)
    with pytest.raises(ValueError, match='creation precision'):
        network.at_precision(4)


def test_cnf_encoding_counts_free_bits_tautologies_duplicates_and_empty_clauses():
    examples = [
        CNF(0, []), CNF(4, []), CNF(0, [[]]),
        CNF(3, [[(0, False), (0, False), (1, True)], [(2, False), (2, True)]]),
        CNF(3, [[(0, False)], [(0, True)]]),
        CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]]),
    ]
    for cnf in examples:
        result = count_cnf_modular(cnf)
        assert result.counted and result.count == brute_count(cnf)
        assert result.theorem_1_1_complete is result.lean_proof_complete is False
    assert result.count == 196 and result.edges == 12 and result.vertices == 14
    assert result.free_variables == 2 and result.boolean_vectors_scanned == 4096
    assert result.tensor_evaluations == 57344


def test_trivalent_equality_tree_and_subdivision_preserve_unique_assignments():
    for degree in (1, 2, 3, 4, 5, 8):
        cnf = CNF(3, [[(0, False)]] * degree)
        encoded = encode_cnf_network(cnf, budgets=ModularBudgets(max_edges=32))
        assert not encoded.refused and encoded.network is not None
        assert encoded.edges == 2 * (degree + max(0, degree - 3))
        assert encoded.network.max_degree == (2 if degree < 3 else 3)
        assert sum(vertex.side == 1 for vertex in encoded.network.vertices) == encoded.edges // 2
        if degree <= 4:
            assert count_cnf_modular(cnf).count == 4
    # Encoding captures its clause data rather than the mutable CNF object.
    cnf = CNF(1, [[(0, False)]])
    encoded = encode_cnf_network(cnf)
    cnf.clauses[0].append((0, True))
    assert contract_network(encoded.network).value == (1,)


def test_random_small_cnf_counts_with_cube_oracle_disabled_in_backend(monkeypatch):
    rng = random.Random(520)
    formulas = [CNF(3, [[(rng.randrange(3), bool(rng.randrange(2))) for _ in range(2)] for _ in range(2)])
                for _ in range(12)]
    expected = [brute_count(cnf) for cnf in formulas]

    def forbidden(*args, **kwargs):
        raise AssertionError('a preexisting counter must not implement the modular backend')

    monkeypatch.setattr('counting_sat_algebra.cnf.brute_count', forbidden)
    monkeypatch.setattr('counting_sat_algebra.component_count.count_cnf_components', forbidden)
    monkeypatch.setattr('counting_sat_algebra.tensor_count.count_placed_formula_by_contraction', forbidden)
    for cnf, count in zip(formulas, expected):
        assert count_cnf_modular(cnf).count == count


@pytest.mark.parametrize('budgets', [
    ModularBudgets(max_edges=3), ModularBudgets(max_vertices=4),
    ModularBudgets(max_boolean_vectors=15), ModularBudgets(max_entry_evaluations=79),
    ModularBudgets(max_precision_bits=1), ModularBudgets(max_ring_coordinates=0),
])
def test_cnf_refusals_precede_network_allocation(monkeypatch, budgets):
    def forbidden(*args, **kwargs):
        raise AssertionError('network allocation started before refusal')

    monkeypatch.setattr('counting_sat_algebra.modular_network.BooleanTensorNetwork', forbidden)
    result = count_cnf_modular(CNF(1, [[(0, False)], [(0, False)]]), budgets=budgets)
    assert result.refused and result.count is None and not result.network_allocated
    assert result.boolean_vectors_scanned == result.tensor_evaluations == 0


def test_leaf_refusal_does_not_sample_implicit_entries():
    def forbidden(bits):
        raise AssertionError('entry evaluation started before refusal')

    network = _parallel_network()
    network = replace(network, vertices=tuple(replace(vertex, entry=forbidden) for vertex in network.vertices))
    result = contract_network(network, budgets=ModularBudgets(max_boolean_vectors=3))
    assert result.refused and result.value is None
    assert result.boolean_vectors_scanned == result.tensor_evaluations == 0


@pytest.mark.parametrize('a', [1, 2])
def test_cauchy_local_coefficients_against_rational_product(a):
    ring = ModularRing(9, 2, 0b111)
    L, alpha = 1 << a, 1 << (a + 1)
    for row, column in product(range(L), repeat=2):
        exact = _rational_cauchy(a, row, column)
        assert exact[0] == int(row == column)
        for t in range(-L, L):
            expected = _exact_shift(exact, t, alpha, 2 * L)
            actual = cauchy_local_polynomial(ring, a, t, row, column)
            assert actual == tuple(ring.from_fraction(value) for value in expected)
            if t != column - row:
                assert all(value[0] % 2 == 0 for value in actual)


def test_deformation_restores_zero_and_has_nonconstant_values():
    network = embed_network(_parallel_network(), 0b10011)
    for layout in (PacketLayout(1, ((0, 1),)), PacketLayout(2, ((None, 1, 0, None),)),
                   PacketLayout(1, ((0, None), (None, 1)))):
        result = evaluate_deformation(network, layout, 0, network.ring.zero)
        assert result.value == contract_network(network).value
        assert result.port_pairs_scanned == 16 and result.tensor_evaluations == 32
    result = evaluate_deformation(network, PacketLayout(1, ((0, 1),)), 1, network.ring.zero)
    assert network.ring.base_coefficients(result.value) == (26, 13)
    assert result.value != contract_network(network).value
    refused = evaluate_deformation(network, PacketLayout(1, ((0, 1),)), 1, network.ring.zero,
                                   budgets=ModularBudgets(max_boolean_vectors=15))
    assert refused.refused and refused.value is None and not refused.kernel_allocated
    assert refused.port_pairs_scanned == refused.tensor_evaluations == 0
    with pytest.raises(ValueError, match='cover every'):
        evaluate_deformation(network, PacketLayout(1, ((0, None),)), 1, network.ring.zero)


def test_tapered_recovery_matches_independent_full_degree_rational_polynomial():
    network = embed_network(_parallel_network(nonlinear=True), 0b10011)
    layout = PacketLayout(1, ((0, 1),))
    plan = plan_tapered_recovery(network, layout, 1)
    polynomials = _nonlinear_deformation_polynomials()
    assert polynomials[0][-1] == Fraction(11, 4) and polynomials[1][-1] == Fraction(155, 4)
    assert plan.degree_bound == 6 and plan.interpolation_points == 7
    samples = _exact_samples(network.ring, plan, polynomials)
    expected = network.ring.from_base((8, 5))
    for node in plan.nodes:
        local = network.at_precision(node.precision)
        points = tuple(local.ring.residue_point(j) for j in range(plan.interpolation_points))
        values = tuple(evaluate_deformation(local, layout, node.t, point).value for point in points)
        all_coefficients = interpolate_polynomial(local.ring, points, values)
        exact = [_exact_shift(p, node.t, 4, plan.interpolation_points) for p in polynomials]
        assert all_coefficients == tuple(
            local.ring.from_base([local.ring.from_fraction(jet[j])[0] for jet in exact])
            for j in range(plan.interpolation_points)
        )
    assert recover_tapered_jets(network.ring, plan, samples) == expected
    assert contract_network(network).value == expected
    recovered = recover_network_modular(network, layout, 1)
    assert recovered.recovered and recovered.value == expected
    assert recovered.port_pairs_scanned == 224 and recovered.tensor_evaluations == 448
    assert recovered.recursive_contraction is recovered.theorem_1_1_complete is False
    # Arbitrary lifts of the same residues give the same final answer.
    altered = [replace(sample, coefficients=tuple(tuple(value - sample.ring.modulus for value in coefficient)
                                                  for coefficient in sample.coefficients)) for sample in samples]
    assert recover_tapered_jets(network.ring, plan, altered) == expected
    with pytest.raises(ValueError, match='precision is below'):
        recover_tapered_jets(network.ring, plan, [replace(samples[0], ring=samples[0].ring.at_precision(5)), samples[1]])
    with pytest.raises(ValueError, match='missing tapered'):
        recover_tapered_jets(network.ring, plan, [replace(samples[0], coefficients=samples[0].coefficients[:1]), samples[1]])
    with pytest.raises(ValueError, match='exactly one jet'):
        recover_tapered_jets(network.ring, plan, samples[:1])


def test_multiple_valuation_levels_and_padded_packets_recover_all_v_coefficients():
    network = embed_network(_parallel_network(), 0b10011)
    layout = PacketLayout(2, ((0, None, 1, None),))
    plan = plan_tapered_recovery(network, layout, 2)
    assert plan.h == 3 and plan.degree_bound == 14 and plan.interpolation_points == 15
    assert [node.t for node in plan.nodes] == [-3, -2, -1, 1, 2, 3]
    assert [node.multiplicity for node in plan.nodes] == [6, 3, 6, 6, 3, 6]
    assert [node.precision for node in plan.nodes] == [149, 101, 149, 149, 101, 149]
    recovered = recover_network_modular(network, layout, 2)
    assert recovered.value == network.ring.from_base((9, 4))
    assert recovered.port_pairs_scanned == plan.port_pair_bound == 1440
    assert recovered.tensor_evaluations == plan.tensor_evaluation_bound == 2880
    with pytest.raises(ValueError, match='Lemma D.5 requires'):
        recover_network_modular(network, layout, 3)


@pytest.mark.parametrize('budgets', [
    ModularBudgets(max_edges=1), ModularBudgets(max_vertices=1),
    ModularBudgets(max_boolean_vectors=223), ModularBudgets(max_entry_evaluations=447),
    ModularBudgets(max_precision_bits=132), ModularBudgets(max_ring_coordinates=7),
    ModularBudgets(max_field_degree=3), ModularBudgets(max_packet_slots=1),
    ModularBudgets(max_interpolation_points=6), ModularBudgets(max_jet_slots=751),
    ModularBudgets(max_recovery_bits=840),
])
def test_recovery_refusals_precede_jets_kernels_and_entry_calls(budgets):
    def forbidden(bits):
        raise AssertionError('entry routine was called after a preflight refusal')

    network = embed_network(_parallel_network(), 0b10011)
    network = replace(network, vertices=tuple(replace(vertex, entry=forbidden) for vertex in network.vertices))
    result = recover_network_modular(network, PacketLayout(1, ((0, 1),)), 1, budgets=budgets)
    assert result.refused and result.value is None and not result.jets_allocated
    assert result.port_pairs_scanned == result.tensor_evaluations == 0


def test_recovery_requires_enough_distinct_residue_points():
    result = recover_network_modular(_parallel_network(), PacketLayout(1, ((0, 1),)), 1)
    assert result.refused and 'too few distinct' in result.reason and result.value is None
    assert result.port_pairs_scanned == 0 and not result.jets_allocated


def test_step20_cli_report_keeps_finite_scope_and_observed_results(tmp_path, capsys):
    from counting_sat_algebra import __version__
    from counting_sat_algebra.__main__ import main

    output = tmp_path / 'step20.json'
    assert main(['verify', '--step20', '--output', str(output)]) == 0
    report = json.loads(output.read_text(encoding='utf-8'))
    assert report['version'] == __version__ and report['steps'] == [20]
    assert report['all_checks_passed'] and report['failed_step'] is None
    assert report['theorem_1_1_complete'] is report['lean_proof_complete'] is False
    assert len(report['rows']) == 14 and all(row['ok'] for row in report['rows'])
    assert {row['label'] for row in report['rows']} == {'COUNTED', 'CHECKED', 'REFUSED'}
    rows = {row['name']: row for row in report['rows']}
    assert rows['appendix-d-cnf-network']['count'] == 196
    assert rows['appendix-d-tapered-precision-levels']['port_pairs_scanned'] == 1440
    assert all(row['recursive_contraction'] is False for row in report['rows'])
    assert '14/14 checks passed' in capsys.readouterr().out
