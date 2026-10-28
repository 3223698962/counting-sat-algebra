"""Step 21 finite checks: collision packing and prospective child topologies.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction

from .protocol import _row, stamp_completion
from ..cnf import CNF, brute_count
from ..modular_ring import ModularBudgetRefusal, ModularRing, extension_ring_for_size
from ..modular_network import BooleanTensorNetwork, TensorVertex
from ..modular_recovery import embed_network, recover_network_modular
from ..modular_packing import (
    PackingBudgets,
    STEP21_SCOPE,
    appendix_d_integers,
    appendix_d_proposal,
    child_topology,
    choose_collision_labels,
    pack_network,
    packing_hypotheses,
)


def run_packet_packing_suite() -> list[dict]:
    print('=== (W) Appendix D collision labels, packets and child topology ===', flush=True)
    rows = []
    proposal = appendix_d_proposal(16, 20)
    integers = proposal.integers
    rows.append(
        _row(
            'appendix-d-active-schedule-guard',
            0,
            20,
            (integers.lambda_s, integers.u, integers.a, proposal.H, proposal.delta)
            == (4, 0, 0, 8, 13)
            and not proposal.leaf
            and proposal.hypotheses is None,
            'CHECKED',
            s=16,
            lambda_s=4,
            u=0,
            a=0,
            H=8,
            delta=13,
            L_formed=False,
            D_arr_formed=False,
            reason=proposal.reason,
        )
    )
    witness = appendix_d_integers(131072)
    deeper = appendix_d_integers(262144, 1)
    rows.append(
        _row(
            'appendix-d-logarithmic-schedule-witness',
            0,
            0,
            (witness.u, witness.a, witness.previous_a) == (1, 2, 0)
            and (deeper.u, deeper.a, deeper.previous_a) == (1, 6, 4)
            and witness.r == 64
            and witness.epsilon == Fraction(1, 131072),
            'CHECKED',
            lambda_s=131072,
            depth=0,
            u=1,
            a=2,
            previous_a=0,
            second_lambda=262144,
            second_depth=1,
            second_u=1,
            second_a=6,
            second_previous_a=4,
            r=64,
            epsilon='1/131072',
            size_s_constructed=False,
        )
    )
    admitted = packing_hypotheses(128, 2, 1, 4)
    too_small = packing_hypotheses(32, 2, 1, 4)
    paper_r = packing_hypotheses(8192, 64, 1, 4)
    rows.append(
        _row(
            'appendix-d-squared-collision-guard',
            0,
            0,
            admitted.admitted
            and not too_small.collision
            and not paper_r.collision
            and paper_r.density
            and paper_r.packet_degree,
            'CHECKED',
            finite_r2_threshold=128,
            r64_D_arr=4096,
            r64_collision_threshold=134217728,
            r64_M=8192,
            r64_D6_admitted=False,
            note='The collision condition contains D_arr squared. These finite tuples are not D.4 schedule members.',
        )
    )

    def no_entry(bits):
        raise AssertionError('collision packing must not evaluate tensor entries')

    ring = ModularRing(1)

    def matching(M):
        return BooleanTensorNetwork(
            ring,
            M,
            tuple(
                TensorVertex(f'{side}-{edge}', side, (edge,), ring, no_entry)
                for side in (0, 1)
                for edge in range(M)
            ),
        )

    ends = ((0, 2), (0, 2), (0, 3), (1, 2), (1, 3), (1, 3))
    small = BooleanTensorNetwork(
        ring,
        len(ends),
        tuple(
            TensorVertex(
                str(v),
                int(v >= 2),
                tuple(i for i, edge in enumerate(ends) if v in edge),
                ring,
                no_entry,
            )
            for v in range(4)
        ),
    )
    labeling = choose_collision_labels(small, 2, 1, trace=True)

    def pair_oracle(assigned):
        value = Fraction(-(len(ends) ** 2), 4)
        for first in ends:
            for second in ends:
                probability = Fraction(1)
                for left, right in zip(first, second):
                    if left == right:
                        factor = Fraction(1)
                    elif left in assigned and right in assigned:
                        factor = Fraction(int(assigned[left] == assigned[right]))
                    else:
                        factor = Fraction(1, 2)
                    probability *= factor
                value += probability
        return value

    assigned = {}
    ce_ok = labeling.initial_potential == Fraction(23, 2) and labeling.final_potential == 1
    for choice in labeling.choices:
        ce_ok &= choice.before == pair_oracle(assigned)
        ce_ok &= all(
            value == pair_oracle({**assigned, choice.vertex: label})
            for label, value in enumerate(choice.candidates)
        )
        ce_ok &= choice.after == min(choice.candidates) <= choice.before
        assigned[choice.vertex] = choice.label
    rows.append(
        _row(
            'appendix-d-shared-parallel-collision-expectation',
            0,
            small.edges,
            ce_ok,
            'CHECKED',
            initial_potential=str(labeling.initial_potential),
            final_potential=str(labeling.final_potential),
            labels=list(labeling.labels),
            pair_loads=[list(row) for row in labeling.pair_loads],
            candidate_expectations_checked=True,
            label_trials=labeling.label_trials,
            score_updates=labeling.score_updates,
            tensor_evaluations=0,
        )
    )

    finite = pack_network(matching(128), 2, 1, delta=4)
    finite_ok = (
        finite.packed
        and finite.check.ok
        and finite.plan.hypotheses.admitted
        and finite.packing.H_pad == 32
        and len(finite.packing.layout.packets) == 64
        and finite.packing.labeling.final_potential == 0
        and finite.check.max_child_edges == 128
        and finite.check.max_packet_degree == 4
        and not finite.check.child_shrink
    )
    rows.append(
        _row(
            'appendix-d-finite-D6-packing',
            0,
            128,
            finite_ok,
            'CHECKED',
            r=2,
            a=1,
            delta=4,
            D6_admitted=True,
            D4_schedule_member=False,
            H_pad=finite.packing.H_pad if finite.packing else None,
            packets=len(finite.packing.layout.packets) if finite.packing else None,
            max_group_ports=64,
            max_child_edges=128,
            max_packet_degree=4,
            child_shrink=False,
            initial_potential='96',
            final_potential='0',
            label_trials=finite.label_trials,
            score_updates=finite.score_updates,
            tensor_evaluations=finite.tensor_evaluations,
        )
    )
    multi = pack_network(matching(32), 2, 2, require_hypotheses=False)
    multi_ok = (
        multi.packed
        and multi.check.ok
        and multi.packing.H_pad == 4
        and len(multi.packing.layout.packets) == 16
        and multi.check.max_packet_degree == 8
        and {sum(low == high for low, high in symbols) for symbols in multi.packing.packet_symbols}
        == {0, 1, 2}
    )
    rows.append(
        _row(
            'appendix-d-multicoordinate-packet-orbits',
            0,
            32,
            multi_ok,
            'CHECKED',
            r=2,
            a=2,
            H_pad=4,
            packets=16,
            max_packet_degree=8,
            equal_coordinate_counts=[0, 1, 2],
            D6_admitted=False,
            tensor_evaluations=multi.tensor_evaluations,
        )
    )

    large = matching(8192)
    strict = pack_network(large, 64, 1, delta=4)
    rows.append(
        _row(
            'appendix-d-r64-D6-refusal',
            0,
            large.edges,
            strict.refused
            and strict.packing is None
            and not strict.labels_allocated
            and strict.label_trials == 0,
            'REFUSED',
            r=64,
            a=1,
            delta=4,
            D6_admitted=False,
            packing=None,
            labels_allocated=False,
            packets_allocated=False,
            label_trials=0,
            reason=strict.reason,
        )
    )
    actual = pack_network(large, 64, 1, delta=4, require_hypotheses=False)
    child = child_topology(actual.packing, 0, 0) if actual.packed else None
    actual_ok = (
        actual.packed
        and actual.check.ok
        and not actual.plan.hypotheses.admitted
        and actual.packing.H_pad == 2
        and len(actual.packing.layout.packets) == 4096
        and actual.packing.labeling.initial_potential == 8190
        and actual.packing.labeling.final_potential == 0
        and all(load == 2 for row in actual.packing.labeling.pair_loads for load in row)
        and actual.check.row_group_ports == actual.check.column_group_ports == ((128,) * 64,)
        and child is not None
        and child.edges == 256
        and child.max_degree == 4
        and actual.check.child_shrink
        and actual.label_trials == 1048576
        and actual.score_updates == 1572864
    )
    rows.append(
        _row(
            'appendix-d-r64-actual-packet-and-child-bounds',
            0,
            large.edges,
            actual_ok,
            'CHECKED',
            r=64,
            a=1,
            delta=4,
            D6_admitted=False,
            D4_schedule_member=False,
            H_pad=2,
            packets=4096,
            initial_potential='8190',
            final_potential='0',
            pair_class_load=2,
            max_group_ports=128,
            max_child_edges=256,
            max_packet_degree=4,
            child_shrink=True,
            label_trials=actual.label_trials,
            score_updates=actual.score_updates,
            tensor_evaluations=actual.tensor_evaluations,
            child_topology_only=True,
        )
    )

    cnf = CNF(2, [[(0, False)], [(0, False)]])
    from ..modular_network import encode_cnf_network

    encoded = encode_cnf_network(cnf)
    field = extension_ring_for_size(2, encoded.network.ring.precision)
    network = embed_network(encoded.network, field.polynomial)
    automatic = pack_network(network, 2, 1, require_hypotheses=False)
    recovered = (
        recover_network_modular(network, automatic.packing.layout, 1) if automatic.packed else None
    )
    count = (
        (recovered.value[0] << encoded.free_variables)
        if recovered is not None and recovered.recovered
        else None
    )
    rows.append(
        _row(
            'appendix-d-automatic-layout-cnf-recovery',
            cnf.n,
            cnf.m,
            count == brute_count(cnf) == 2
            and recovered.port_pairs_scanned == 6656
            and recovered.tensor_evaluations == 33280,
            'COUNTED',
            count=count,
            brute=2,
            network_edges=4,
            free_variables=1,
            port_pairs_scanned=recovered.port_pairs_scanned if recovered else 0,
            tensor_evaluations=recovered.tensor_evaluations if recovered else 0,
            backend='Step 20 local recovery with automatically constructed packets',
            D6_admitted=False,
        )
    )

    for name, limits, trace in (
        ('label-trials', PackingBudgets(max_label_trials=1), False),
        ('trace', PackingBudgets(max_trace_entries=1), True),
        ('padding', PackingBudgets(max_padded_slots=167), False),
        ('score-updates', PackingBudgets(max_score_updates=1023), False),
    ):
        refused = pack_network(matching(128), 2, 1, delta=4, budgets=limits, trace=trace)
        rows.append(
            _row(
                f'appendix-d-packing-{name}-budget',
                0,
                128,
                refused.refused
                and refused.packing is None
                and not refused.labels_allocated
                and not refused.packets_allocated
                and refused.label_trials == refused.score_updates == 0,
                'REFUSED',
                packing=None,
                label_trials=0,
                score_updates=0,
                tensor_evaluations=0,
                labels_allocated=False,
                packets_allocated=False,
                reason=refused.reason,
            )
        )
    child_refused = False
    child_reason = 'packing did not complete'
    if actual.packed:
        try:
            child_topology(actual.packing, 0, 0, budgets=PackingBudgets(max_child_edges=255))
        except ModularBudgetRefusal as error:
            child_refused = True
            child_reason = str(error)
    rows.append(
        _row(
            'appendix-d-child-topology-budget',
            0,
            8192,
            child_refused,
            'REFUSED',
            child_edges=256,
            child_edge_budget=255,
            child_topology_allocated=False,
            packing_already_constructed=True,
            tensor_evaluations=0,
            reason=child_reason,
        )
    )
    for row in rows:
        stamp_completion(row)
        print(f'  [W] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP21_SCOPE}', flush=True)
    return rows
