"""Independent conditional-expectation and packet geometry checks for D.4.2."""

import json
import random
from collections import Counter
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from counting_sat_algebra.cnf import CNF, brute_count
from counting_sat_algebra.modular_network import BooleanTensorNetwork, TensorVertex, encode_cnf_network
from counting_sat_algebra.modular_ring import ModularBudgetRefusal, ModularRing, extension_ring_for_size
from counting_sat_algebra.modular_recovery import (
    cauchy_local_polynomial, embed_network, recover_network_modular, valuation2,
)
from counting_sat_algebra.modular_packing import (
    PackingBudgets, appendix_d_integers, appendix_d_proposal, check_packing,
    child_topology, choose_collision_labels, pack_network, packing_hypotheses,
)


def _forbidden(bits):
    raise AssertionError('packing must not evaluate tensor entries')


def _matching(edges, *, isolates=0):
    ring = ModularRing(1)
    vertices = tuple(TensorVertex(f'{side}-{edge}', side, (edge,), ring, _forbidden)
                     for side in (0, 1) for edge in range(edges))
    vertices += tuple(TensorVertex(f'isolated-{j}', j % 2, (), ring, _forbidden) for j in range(isolates))
    return BooleanTensorNetwork(ring, edges, vertices)


def _multigraph(order=(0, 1, 2, 3)):
    ends = ((0, 2), (0, 2), (0, 3), (1, 2), (1, 3), (1, 3))
    ring = ModularRing(5)
    vertices = tuple(TensorVertex(str(v), 0 if v < 2 else 1,
                                  tuple(i for i, edge in enumerate(ends) if v in edge), ring, _forbidden)
                     for v in order)
    return BooleanTensorNetwork(ring, len(ends), vertices)


def _endpoints(network):
    endpoints = [[None, None] for _ in range(network.edges)]
    for vertex, tensor in enumerate(network.vertices):
        for edge in tensor.ports:
            endpoints[edge][tensor.side] = vertex
    return endpoints


def _exhaustive_expectation(network, assigned, K):
    """Enumerate all remaining vertex labels; never use the CE score formula."""
    missing = [v for v in range(len(network.vertices)) if v not in assigned]
    ends = _endpoints(network)
    total = 0
    for choice in product(range(K), repeat=len(missing)):
        labels = dict(assigned)
        labels.update(zip(missing, choice))
        loads = Counter((labels[row], labels[column]) for row, column in ends)
        total += sum(value * value for value in loads.values())
    return Fraction(total, K ** len(missing)) - Fraction(network.edges ** 2, K * K)


@pytest.mark.parametrize('r,a', [(2, 1), (3, 1), (2, 2)])
@pytest.mark.parametrize('order', [(0, 1, 2, 3), (2, 0, 3, 1), (3, 1, 0, 2)])
def test_every_candidate_expectation_matches_full_label_enumeration(r, a, order):
    network = _multigraph(order)
    result = choose_collision_labels(network, r, a, trace=True)
    K = r ** a
    assigned = {}
    assert result.initial_potential == _exhaustive_expectation(network, assigned, K)
    for step in result.choices:
        assert step.before == _exhaustive_expectation(network, assigned, K)
        assert sum(step.candidates, Fraction()) / K == step.before
        for label, expectation in enumerate(step.candidates):
            assert expectation == _exhaustive_expectation(network, {**assigned, step.vertex: label}, K)
        assert step.label == min(range(K), key=step.candidates.__getitem__)
        assert step.after == min(step.candidates) <= step.before
        assigned[step.vertex] = step.label
    assert result.final_potential == _exhaustive_expectation(network, assigned, K)
    assert result.label_trials == K * len(network.vertices)
    assert choose_collision_labels(network, r, a, trace=True) == result


def test_parallel_and_shared_incidence_corrections_are_present():
    result = choose_collision_labels(_multigraph(), 2, 1, trace=True)
    assert result.initial_potential == Fraction(23, 2)
    assert result.final_potential == 1 and result.labels == (0, 1, 0, 1)
    assert result.pair_loads == ((2, 1), (1, 2))
    assert result.initial_potential != Fraction(6 * 3, 4)  # The independent-edge value is wrong here.


def test_random_graphs_and_interleaved_sides_match_pair_probability_oracle():
    rng = random.Random(621)
    ring = ModularRing(2)
    for _ in range(12):
        ends = tuple((rng.randrange(3), 3 + rng.randrange(3)) for _ in range(8))
        order = list(range(6))
        rng.shuffle(order)
        network = BooleanTensorNetwork(ring, len(ends), tuple(
            TensorVertex(str(v), int(v >= 3), tuple(i for i, edge in enumerate(ends) if v in edge), ring, _forbidden)
            for v in order
        ))
        result = choose_collision_labels(network, 2, 1, trace=True)
        assigned = {}
        for step in result.choices:
            for label, expected in enumerate(step.candidates):
                actual = _exhaustive_expectation(network, {**assigned, step.vertex: label}, 2)
                assert expected == actual
            assigned[step.vertex] = step.label


def test_d4_integer_schedule_ceilings_and_zero_guard():
    for s in range(2, 80):
        proposal = appendix_d_proposal(s, 2 * s)
        lam = (s - 1).bit_length()
        divisor = max(k for k in range(1, lam + 1) if k * k <= lam)
        assert proposal.H * divisor >= s > (proposal.H - 1) * divisor
        assert proposal.delta ** 10 >= s ** 9 > (proposal.delta - 1) ** 10
        assert proposal.integers.u == proposal.integers.a == 0
        assert proposal.hypotheses is None and 'not formed' in proposal.reason
    assert appendix_d_proposal(16, 8).leaf
    active = appendix_d_proposal(16, 20)
    assert (active.H, active.delta, active.integers.lambda_s) == (8, 13, 4)
    assert not active.leaf
    for lam in (131071, 131072, 262144, 655360):
        for depth in (0, 1, 2, 3, 10):
            integers = appendix_d_integers(lam, depth)
            terms = [lam // (131072 * (1 << j)) for j in range(depth + 1)]
            assert integers.u == terms[-1] and integers.a == 2 * sum(terms)
            assert integers.previous_a == 2 * sum(terms[:-1])
    witness = appendix_d_integers(131072)
    assert (witness.u, witness.a, witness.previous_a) == (1, 2, 0)
    assert witness.r == 64 and witness.epsilon == Fraction(1, 131072)
    assert appendix_d_integers(655360, 10**12).a == 16
    with pytest.raises(ModularBudgetRefusal, match='intermediate bits'):
        appendix_d_proposal(1 << 500, 20)


def test_d6_collision_guard_uses_the_square_of_the_pair_class_count():
    h = packing_hypotheses(128, 2, 1, 4)
    assert (h.L, h.labels, h.D_arr) == (2, 2, 4)
    assert h.admitted and not h.failures
    assert not packing_hypotheses(127, 2, 1, 4).collision
    assert not packing_hypotheses(32, 2, 1, 4).collision
    large = packing_hypotheses(8192, 64, 1, 4)
    assert not large.collision and large.density and large.packet_degree
    assert 2 * large.D_arr ** 2 * large.delta == 134217728
    assert 'D_arr^2' in large.failures[0]


def test_d6_admitted_finite_packing_and_exact_padding():
    result = pack_network(_matching(128), 2, 1, delta=4, trace=True)
    assert result.packed and result.check.ok and result.plan.hypotheses.admitted
    packed = result.packing
    assert packed.H_pad == 32 and len(packed.layout.packets) == 64
    assert packed.labeling.initial_potential == 96 and packed.labeling.final_potential == 0
    assert result.check.row_group_ports == result.check.column_group_ports == ((64, 64),)
    assert result.check.max_child_edges == 128 and result.check.max_packet_degree == 4
    assert not result.check.child_shrink
    assert result.label_trials == 512 and result.score_updates <= result.plan.score_update_bound
    assert result.plan.padded_slot_bound == 168 and result.plan.packet_bound == 84
    assert result.tensor_evaluations == 0 and not result.recursive_contraction


def test_r64_actual_conclusions_and_child_shrink_without_claiming_d6():
    network = _matching(8192)
    refused = pack_network(network, 64, 1, delta=4)
    assert refused.refused and refused.packing is None and not refused.labels_allocated
    result = pack_network(network, 64, 1, delta=4, require_hypotheses=False)
    assert result.packed and result.check.ok and not result.plan.hypotheses.admitted
    assert result.packing.H_pad == 2 and len(result.packing.layout.packets) == 4096
    assert result.packing.labeling.initial_potential == 8190 and result.packing.labeling.final_potential == 0
    assert all(value == 2 for row in result.packing.labeling.pair_loads for value in row)
    assert result.label_trials == 1048576 and result.score_updates == 1572864
    assert result.check.row_group_ports == result.check.column_group_ports == ((128,) * 64,)
    assert result.check.max_child_edges == 256 and result.check.child_shrink
    for group in (0, 63):
        child = child_topology(result.packing, 0, group)
        assert child.edges == 256 <= network.edges // 8
        assert child.max_degree == 4 and len(child.original_vertices) == 256
    assert result.tensor_evaluations == 0
    assert result.theorem_1_1_complete is result.lean_proof_complete is False


def test_multicoordinate_equal_and_unequal_orbits_cover_every_edge_once():
    result = pack_network(_matching(32, isolates=2), 2, 2, require_hypotheses=False)
    assert result.packed and result.check.ok
    packed = result.packing
    assert packed.H_pad == 4 and len(packed.layout.packets) == 16
    assert sorted(edge for packet in packed.layout.packets for edge in packet if edge is not None) == list(range(32))
    equal_counts = {sum(low == high for low, high in symbols) for symbols in packed.packet_symbols}
    assert equal_counts == {0, 1, 2}
    for coordinate in range(2):
        children = [child_topology(packed, coordinate, g) for g in range(2)]
        assert sum(child.edges for child in children) == 64
        assert sorted(v for child in children for v in child.original_vertices) == list(range(66))
        for child in children:
            assert child.max_degree <= 8
            assert sorted(e for ports in child.original_ports for e in ports) == list(range(child.edges))
            assert sorted(e for ports in child.packet_ports for e in ports) == list(range(child.edges))
            for vertex, ports in zip(child.original_vertices, child.original_ports):
                assert tuple(child.port_origins[e] for e in ports) == tuple((vertex, edge) for edge in packed.source.vertices[vertex].ports)


def test_packet_orientation_implies_locality_at_nonzero_cauchy_nodes():
    packed = pack_network(_matching(32), 2, 2, require_hypotheses=False).packing
    endpoints = _endpoints(packed.source)
    ring = ModularRing(9)
    for packet in packed.layout.packets:
        for x, row_edge in enumerate(packet):
            if row_edge is None:
                continue
            for y, col_edge in enumerate(packet):
                if col_edge is None:
                    continue
                for t in (-3, -2, -1, 1, 2, 3):
                    coordinate = valuation2(t)
                    if packed.symbol(endpoints[row_edge][0], coordinate) != packed.symbol(endpoints[col_edge][1], coordinate):
                        assert y - x != t
                        assert cauchy_local_polynomial(ring, 2, t, x, y)[0] == ring.zero


def test_automatic_layout_connects_to_modular_recovery_on_a_connected_cnf_network():
    cnf = CNF(2, [[(0, False)], [(0, False)]])
    encoded = encode_cnf_network(cnf)
    extension = extension_ring_for_size(2, encoded.network.ring.precision)
    network = embed_network(encoded.network, extension.polynomial)
    packed = pack_network(network, 2, 1, require_hypotheses=False)
    assert packed.packed and network.edges == 4
    recovered = recover_network_modular(network, packed.packing.layout, 1)
    assert recovered.recovered and extension.base_coefficients(recovered.value) == (1,)
    assert (recovered.value[0] << encoded.free_variables) == brute_count(cnf) == 2
    assert recovered.port_pairs_scanned == 6656 and recovered.tensor_evaluations == 33280


@pytest.mark.parametrize('budgets', [
    PackingBudgets(max_edges=127), PackingBudgets(max_vertices=255),
    PackingBudgets(max_labels=1), PackingBudgets(max_pair_classes=3),
    PackingBudgets(max_packet_slots=1), PackingBudgets(max_label_trials=511),
    PackingBudgets(max_score_updates=1023), PackingBudgets(max_padded_slots=167),
    PackingBudgets(max_packets=83), PackingBudgets(max_trace_entries=511),
    PackingBudgets(max_parameter_bits=7),
])
def test_all_resource_refusals_precede_label_search_and_packet_allocation(monkeypatch, budgets):
    def forbidden(*args, **kwargs):
        raise AssertionError('label or packet construction began after a preflight refusal')

    monkeypatch.setattr('counting_sat_algebra.modular_packing._choose_labels', forbidden)
    monkeypatch.setattr('counting_sat_algebra.modular_packing._make_packets', forbidden)
    result = pack_network(_matching(128), 2, 1, delta=4, trace=True, budgets=budgets)
    assert result.refused and result.packing is None
    assert not result.labels_allocated and not result.packets_allocated
    assert result.label_trials == result.score_updates == result.tensor_evaluations == 0


def test_default_parameter_refusal_and_failed_actual_conclusions_are_distinct():
    tiny = _matching(1)
    default = pack_network(tiny, 2, 1)
    assert default.refused and default.label_trials == 0 and default.packing is None
    finite = pack_network(tiny, 2, 1, require_hypotheses=False)
    assert finite.refused and finite.label_trials == 4 and finite.packing is None
    assert finite.labels_allocated and finite.packets_allocated
    assert 'packet count' in finite.reason and finite.tensor_evaluations == 0


def test_checker_rejects_forged_geometry_potentials_and_admission_flags():
    packed = pack_network(_matching(32), 2, 2, require_hypotheses=False).packing
    assert check_packing(packed).ok
    assert not check_packing(replace(packed, hypotheses=replace(packed.hypotheses, collision=True))).ok
    wrong = replace(packed.labeling, initial_potential=Fraction(-1))
    assert not check_packing(replace(packed, labeling=wrong)).ok
    wrong = replace(packed.labeling, pair_loads=((32,),))
    assert not check_packing(replace(packed, labeling=wrong)).ok
    symbols = list(packed.packet_symbols)
    symbols[0] = ((1, 1), (1, 1))
    tampered = replace(packed, packet_symbols=tuple(symbols))
    assert not check_packing(tampered).ok
    with pytest.raises(ValueError, match='invalid packing'):
        child_topology(tampered, 0, 0)


def test_child_topology_budgets_and_empty_scalar_vertices():
    packed = pack_network(_matching(4), 2, 1, require_hypotheses=False).packing
    with pytest.raises(ModularBudgetRefusal, match='child edges'):
        child_topology(packed, 0, 0, budgets=PackingBudgets(max_child_edges=1))
    with pytest.raises(ModularBudgetRefusal, match='child vertices'):
        child_topology(packed, 0, 0, budgets=PackingBudgets(max_vertices=1))
    with pytest.raises(ValueError, match='outside'):
        child_topology(packed, 1, 0)
    empty = pack_network(_matching(0, isolates=2), 2, 1, require_hypotheses=False)
    assert empty.packed and not empty.packing.layout.packets
    child = child_topology(empty.packing, 0, 0)
    assert child.edges == 0 and child.original_vertices == (0, 1)
    assert child.original_ports == ((), ())


def test_step21_report_records_actual_bounds_and_keeps_theorem_flags_false(tmp_path, capsys):
    from counting_sat_algebra import __version__
    from counting_sat_algebra.__main__ import main

    output = tmp_path / 'step21.json'
    assert main(['verify', '--step21', '--output', str(output)]) == 0
    report = json.loads(output.read_text(encoding='utf-8'))
    assert report['version'] == __version__ and report['steps'] == [21]
    assert report['all_checks_passed'] and report['failed_step'] is None
    assert report['theorem_1_1_complete'] is report['lean_proof_complete'] is False
    assert len(report['rows']) == 14 and all(row['ok'] for row in report['rows'])
    assert {row['label'] for row in report['rows']} == {'COUNTED', 'CHECKED', 'REFUSED'}
    rows = {row['name']: row for row in report['rows']}
    assert rows['appendix-d-r64-actual-packet-and-child-bounds']['max_child_edges'] == 256
    assert rows['appendix-d-r64-actual-packet-and-child-bounds']['D6_admitted'] is False
    assert rows['appendix-d-automatic-layout-cnf-recovery']['count'] == 2
    assert all(row['recursive_contraction'] is False for row in report['rows'])
    assert '14/14 checks passed' in capsys.readouterr().out
