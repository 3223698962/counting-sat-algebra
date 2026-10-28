"""Step 20 finite checks: Appendix D rings, networks, and local recovery.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from dataclasses import asdict

from .protocol import _row, stamp_completion
from ..cnf import CNF, brute_count
from ..modular_ring import (
    ModularRing,
    evaluate_polynomial,
    extension_ring_for_size,
    interpolate_polynomial,
)
from ..modular_network import (
    BooleanTensorNetwork,
    ModularBudgets,
    STEP20_SCOPE,
    TensorVertex,
    contract_network,
    count_cnf_modular,
)
from ..modular_recovery import (
    PacketLayout,
    cauchy_local_polynomial,
    embed_network,
    evaluate_deformation,
    recover_network_modular,
)


def run_appendix_d_suite() -> list[dict]:
    print('=== (V) Appendix D finite modular networks and local recovery ===', flush=True)
    rows = []
    ring = extension_ring_for_size(2, 5, 2)
    unit = ring.element((1, 1, 0, 0, 2))
    V = ring.from_base((0, 1))
    points = tuple(ring.residue_point(j) for j in range(7))
    coefficients = tuple(ring.element((j, j + 1, 0, 0, 2 * j)) for j in range(7))
    values = tuple(evaluate_polynomial(ring, coefficients, point) for point in points)
    ring_ok = (
        ring.extension_degree == 4
        and ring.polynomial == 0b10011
        and ring.mul(unit, ring.inverse(unit)) == ring.one
        and ring.mul(V, V) == ring.zero
        and interpolate_polynomial(ring, points, values) == coefficients
    )
    rows.append(
        _row(
            'appendix-d-coefficient-ring',
            0,
            0,
            ring_ok,
            'CHECKED',
            s=2,
            precision=5,
            ell=2,
            extension_degree=4,
            polynomial_bits=bin(ring.polynomial),
            interpolation_points=7,
            all_coordinates_checked=True,
        )
    )

    demo = CNF(8, [[(0, False), (1, False), (2, False)], [(3, False), (4, False), (5, False)]])
    counted = count_cnf_modular(demo)
    rows.append(
        _row(
            'appendix-d-cnf-network',
            demo.n,
            demo.m,
            counted.counted
            and counted.count == brute_count(demo) == 196
            and counted.edges == 12
            and counted.vertices == 14
            and counted.max_degree == 3
            and counted.free_variables == 2
            and counted.boolean_vectors_scanned == 4096
            and counted.tensor_evaluations == 57344,
            'COUNTED',
            **asdict(counted),
            clause_ports=6,
            backend='direct modular leaf',
        )
    )
    tree = CNF(3, [[(0, False)]] * 4)
    tree_count = count_cnf_modular(tree)
    rows.append(
        _row(
            'appendix-d-equality-tree',
            tree.n,
            tree.m,
            tree_count.counted
            and tree_count.count == brute_count(tree) == 4
            and tree_count.edges == 10
            and tree_count.vertices == 11
            and tree_count.max_degree == 3,
            'COUNTED',
            **asdict(tree_count),
            clause_ports=4,
            trivalent_equality_vertices=2,
        )
    )
    empty, unsat = count_cnf_modular(CNF(5, [])), count_cnf_modular(CNF(5, [[]]))
    rows.append(
        _row(
            'appendix-d-trivial-formulas',
            5,
            0,
            empty.counted
            and empty.count == 32
            and unsat.counted
            and unsat.count == 0
            and empty.edges == unsat.edges == 0
            and unsat.network_allocated is False,
            'COUNTED',
            counts=[empty.count, unsat.count],
            edges=0,
            tensor_evaluations=0,
            empty_network_assignments=empty.boolean_vectors_scanned,
        )
    )

    creation, later = ModularRing(2), ModularRing(8)
    lifted = BooleanTensorNetwork(
        later,
        1,
        (
            TensorVertex('row', 0, (0,), creation, lambda bits: (-1,)),
            TensorVertex('column', 1, (0,), creation, lambda bits: (6,)),
        ),
    )
    lifted_result = contract_network(lifted)
    rows.append(
        _row(
            'appendix-d-creation-precision',
            0,
            1,
            lifted_result.value == (12,),
            'CHECKED',
            creation_precision=2,
            call_precision=8,
            value=list(lifted_result.value or ()),
            lifted_entry_values=[3, 2],
            boolean_vectors_scanned=lifted_result.boolean_vectors_scanned,
            note='Each entry is reduced at precision 2 and lifted separately; the two edge terms sum to 12.',
        )
    )

    local_ring = ring.at_precision(9)
    kernel = cauchy_local_polynomial(local_ring, 1, 1, 0, 1)
    expected_kernel = tuple(local_ring.scalar(value) for value in (-3, -22, -48, -32))
    identity_ok = all(
        cauchy_local_polynomial(local_ring, 1, 0, x, y)[0] == local_ring.scalar(int(x == y))
        for x in range(2)
        for y in range(2)
    )
    rows.append(
        _row(
            'appendix-d-integral-cauchy-coordinates',
            0,
            0,
            kernel == expected_kernel and identity_ok,
            'CHECKED',
            a=1,
            L=2,
            alpha=4,
            t=1,
            row=0,
            column=1,
            exact_scalar_coefficients=[-3, -22, -48, -32],
            precision=9,
            identity_at_zero_checked=True,
        )
    )

    base = ModularRing(5, 2)
    network = BooleanTensorNetwork(
        base,
        2,
        (
            TensorVertex('row', 0, (0, 1), base, lambda bits: (1 + bits[0], 2 - bits[1])),
            TensorVertex('column', 1, (0, 1), base, lambda bits: (1 + bits[1], -bits[0])),
        ),
    )
    extension = embed_network(network, ring.polynomial)
    layout = PacketLayout(1, ((0, 1),))
    recovered = recover_network_modular(extension, layout, 1)
    at_one = evaluate_deformation(extension, layout, 1, ring.zero)
    recovered_ok = (
        recovered.recovered
        and recovered.value == contract_network(extension).value == ring.from_base((9, 4))
        and at_one.value == ring.from_base((26, 13))
        and recovered.port_pairs_scanned == 224
        and recovered.tensor_evaluations == 448
        and recovered.plan is not None
        and recovered.plan.degree_bound == 6
        and [node.precision for node in recovered.plan.nodes] == [133, 133]
    )
    rows.append(
        _row(
            'appendix-d-local-recovery',
            0,
            2,
            recovered_ok,
            'CHECKED',
            coefficients=list(ring.base_coefficients(recovered.value))
            if recovered.value is not None
            else None,
            value_at_one=[26, 13],
            precision=5,
            ell=2,
            plan=asdict(recovered.plan) if recovered.plan is not None else None,
            port_pairs_scanned=recovered.port_pairs_scanned,
            tensor_evaluations=recovered.tensor_evaluations,
            recursive_contraction=False,
        )
    )

    nonlinear = BooleanTensorNetwork(
        base,
        2,
        (
            TensorVertex('row', 0, (0, 1), base, lambda bits: (1 + bits[0] * bits[1], 2 - bits[1])),
            network.vertices[1],
        ),
    )
    nonlinear = embed_network(nonlinear, ring.polynomial)
    full_degree = recover_network_modular(nonlinear, layout, 1)
    rows.append(
        _row(
            'appendix-d-nonlinear-row-recovery',
            0,
            2,
            full_degree.recovered
            and full_degree.value == contract_network(nonlinear).value == ring.from_base((8, 5)),
            'CHECKED',
            coefficients=[8, 5],
            precision=5,
            ell=2,
            port_pairs_scanned=full_degree.port_pairs_scanned,
            tensor_evaluations=full_degree.tensor_evaluations,
            note='The degree-6 rational deformation is checked coefficientwise by the independent pytest oracle.',
        )
    )

    tapered = recover_network_modular(extension, PacketLayout(2, ((0, None, 1, None),)), 2)
    tapered_ok = (
        tapered.recovered
        and tapered.value == ring.from_base((9, 4))
        and tapered.plan is not None
        and [node.multiplicity for node in tapered.plan.nodes] == [6, 3, 6, 6, 3, 6]
        and [node.precision for node in tapered.plan.nodes] == [149, 101, 149, 149, 101, 149]
        and tapered.port_pairs_scanned == 1440
        and tapered.tensor_evaluations == 2880
    )
    rows.append(
        _row(
            'appendix-d-tapered-precision-levels',
            0,
            2,
            tapered_ok,
            'CHECKED',
            coefficients=[9, 4],
            precision=5,
            ell=2,
            plan=asdict(tapered.plan) if tapered.plan is not None else None,
            port_pairs_scanned=tapered.port_pairs_scanned,
            tensor_evaluations=tapered.tensor_evaluations,
            recursive_contraction=False,
        )
    )

    refused = count_cnf_modular(demo, budgets=ModularBudgets(max_boolean_vectors=4095))
    rows.append(
        _row(
            'appendix-d-network-vector-budget',
            demo.n,
            demo.m,
            refused.refused
            and refused.count is None
            and not refused.network_allocated
            and refused.boolean_vectors_scanned == refused.tensor_evaluations == 0,
            'REFUSED',
            **asdict(refused),
        )
    )
    for name, budgets in (
        ('precision', ModularBudgets(max_precision_bits=132)),
        ('port-pairs', ModularBudgets(max_boolean_vectors=223)),
        ('jet-storage', ModularBudgets(max_jet_slots=751)),
    ):
        refused_recovery = recover_network_modular(extension, layout, 1, budgets=budgets)
        rows.append(
            _row(
                f'appendix-d-recovery-{name}-budget',
                0,
                2,
                refused_recovery.refused
                and refused_recovery.value is None
                and not refused_recovery.jets_allocated
                and refused_recovery.port_pairs_scanned == refused_recovery.tensor_evaluations == 0,
                'REFUSED',
                value=None,
                reason=refused_recovery.reason,
                port_pairs_scanned=0,
                tensor_evaluations=0,
                jets_allocated=False,
            )
        )
    short_field = recover_network_modular(network, layout, 1)
    rows.append(
        _row(
            'appendix-d-insufficient-residue-points',
            0,
            2,
            short_field.refused
            and short_field.value is None
            and 'too few distinct' in short_field.reason
            and short_field.port_pairs_scanned == 0
            and not short_field.jets_allocated,
            'REFUSED',
            value=None,
            reason=short_field.reason,
            field_size=2,
            required_points=7,
            port_pairs_scanned=0,
            jets_allocated=False,
        )
    )
    for row in rows:
        stamp_completion(row)
        print(f'  [V] {row["name"]}: {row["label"]}', flush=True)
    print(f'  {STEP20_SCOPE}', flush=True)
    return rows
