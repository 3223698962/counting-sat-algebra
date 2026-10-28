"""Finite Cauchy deformation and tapered Hermite recovery (D.4.3–D.4.4).

The evaluation of (D.8) still streams all independent row/column port pairs.
It does not use D.4.5 mode tuples or the child networks of (D.16).  All nodes,
precisions, interpolation storage and pair scans are planned before evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial
from typing import Sequence

from .modular_network import (
    BooleanTensorNetwork, ModularBudgets, TensorVertex, _side_product,
    bounded_power_of_two, check_limit, check_ring_budget,
)
from .modular_ring import (
    Element, ModularBudgetRefusal, ModularRing, interpolate_polynomial,
)


def valuation2(value: int) -> int:
    if value == 0:
        raise ValueError("zero has no finite two-adic valuation")
    value = abs(value)
    return (value & -value).bit_length() - 1


@dataclass(frozen=True)
class PacketLayout:
    """Same edge slot at both ends, with None for a padded, unused slot.

    This is a supplied finite layout.  It has no symbol labels and does not
    assert Lemma D.3's collision or group-load bounds.
    """

    a: int
    packets: tuple[tuple[int | None, ...], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "packets", tuple(tuple(packet) for packet in self.packets))
        if not isinstance(self.a, int) or self.a < 1:
            raise ValueError("packet exponent a must be a positive integer")
        for packet in self.packets:
            length = len(packet)
            if length == 0 or length & (length - 1) or length.bit_length() - 1 != self.a:
                raise ValueError("each packet must have exactly 2^a slots")

    def locations(self, edges: int) -> dict[int, tuple[int, int]]:
        locations = {}
        for packet_index, packet in enumerate(self.packets):
            for slot, edge in enumerate(packet):
                if edge is None:
                    continue
                if not isinstance(edge, int) or not 0 <= edge < edges or edge in locations:
                    raise ValueError("each real edge must occupy one distinct packet slot")
                locations[edge] = (packet_index, slot)
        if len(locations) != edges:
            raise ValueError("packet layout does not cover every network edge")
        return locations


def embed_network(network: BooleanTensorNetwork, polynomial: int) -> BooleanTensorNetwork:
    """Embed base-ring input tensors into the D.4.1 coefficient extension."""
    if polynomial == network.ring.polynomial:
        return network
    if network.ring.extension_degree != 1:
        raise ValueError("only base-ring networks can change coefficient extension")
    target = ModularRing(network.ring.precision, network.ring.ell, polynomial)
    vertices = []
    for vertex in network.vertices:
        creation = target.at_precision(vertex.creation_ring.precision)

        def entry(bits, original=vertex, ring=creation):
            source = original.creation_ring.element(original.entry(bits))
            return ring.from_base(original.creation_ring.base_coefficients(source))

        vertices.append(TensorVertex(vertex.name, vertex.side, vertex.ports, creation, entry))
    return BooleanTensorNetwork(target, network.edges, tuple(vertices))


def _cauchy_factors(a: int, t: int, delta: int):
    L, alpha = 1 << a, 1 << (a + 1)
    if not -L <= t < L or not -(L - 1) <= delta < L:
        raise ValueError("Cauchy node or slot difference outside the grid")
    denominator = (-1) ** (L - 1) * factorial(L) * factorial(L - 1)
    powers = -valuation2(denominator)
    factors = []
    for q in range(-L, L):
        if q == delta:
            continue
        content = a + 1 if q == t else valuation2(t - q)
        powers += content
        factors.append(((t - q) // (1 << content), alpha // (1 << content)))
    if powers < 0:
        raise ArithmeticError("Cauchy local coefficients are not two-adically integral")
    odd_part = denominator // (1 << valuation2(denominator))
    return tuple(factors), powers, odd_part


def _cauchy_entry(ring: ModularRing, a: int, t: int, delta: int, zeta: Element) -> Element:
    factors, powers, odd_part = _cauchy_factors(a, t, delta)
    value = ring.scalar((1 << powers) * pow(odd_part, -1, ring.modulus))
    for constant, slope in factors:
        value = ring.mul(value, ring.add(ring.scalar(constant), ring.scale(zeta, slope)))
    return value


def cauchy_local_polynomial(
    ring: ModularRing, a: int, t: int, row: int, column: int,
    *, budgets: ModularBudgets | None = None,
) -> tuple[Element, ...]:
    """Coefficients of W_xy(t+2^(a+1) z), canceling powers of two first."""
    limits = budgets or ModularBudgets()
    if not isinstance(a, int) or a < 1:
        raise ValueError("a must be a positive integer")
    L = bounded_power_of_two(a, limits.max_packet_slots, "packet slots")
    check_ring_budget(ring, limits)
    check_limit(2 * L * ring.size, limits.max_jet_slots, "Cauchy coefficient slots")
    if not 0 <= row < L or not 0 <= column < L:
        raise ValueError("row and column must be real packet slot indices")
    factors, powers, odd_part = _cauchy_factors(a, t, column - row)
    coefficients = [(1 << powers) * pow(odd_part, -1, ring.modulus) % ring.modulus]
    for constant, slope in factors:
        updated = [0] * (len(coefficients) + 1)
        for j, value in enumerate(coefficients):
            updated[j] += constant * value
            updated[j + 1] += slope * value
        coefficients = [value % ring.modulus for value in updated]
    return tuple(ring.scalar(value) for value in coefficients)


@dataclass(frozen=True)
class DeformationResult:
    value: Element | None
    refused: bool
    reason: str
    port_pairs_scanned: int = 0
    tensor_evaluations: int = 0
    kernel_allocated: bool = False


def _deformation_value(network, layout, t, zeta, locations, stats):
    ring, M = network.ring, network.edges
    kernel = [[ring.zero] * M for _ in range(M)]
    for row in range(M):
        packet, x = locations[row]
        for column in range(M):
            other_packet, y = locations[column]
            if packet == other_packet:
                kernel[row][column] = _cauchy_entry(ring, layout.a, t, y - x, zeta)
    total = ring.zero
    for mask in range(1 << (2 * M)):
        rows = tuple((mask >> edge) & 1 for edge in range(M))
        columns = tuple((mask >> (M + edge)) & 1 for edge in range(M))
        value = ring.mul(_side_product(network, rows, 0, stats), _side_product(network, columns, 1, stats))
        for row in range(M):
            form = ring.zero
            for column, bit in enumerate(columns):
                if bit:
                    form = ring.add(form, kernel[row][column])
            value = ring.mul(value, form if rows[row] else ring.sub(ring.one, form))
        total = ring.add(total, value)
        stats['port_pairs_scanned'] += 1
    return total


def evaluate_deformation(
    network: BooleanTensorNetwork, layout: PacketLayout, t: int, zeta: Element,
    *, budgets: ModularBudgets | None = None,
) -> DeformationResult:
    """Evaluate (D.8) by a finite port-pair scan; no tensor is sampled on refusal."""
    limits = budgets or ModularBudgets()
    locations = layout.locations(network.edges)
    network.ring._check(zeta)
    try:
        L = bounded_power_of_two(layout.a, limits.max_packet_slots, "packet slots")
        if not -L <= t < L:
            raise ValueError("node outside the Cauchy grid")
        check_ring_budget(network.ring, limits)
        check_limit(network.edges, limits.max_edges, "network edges")
        check_limit(len(network.vertices), limits.max_vertices, "network vertices")
        pairs = bounded_power_of_two(2 * network.edges, limits.max_boolean_vectors, "port pairs")
        check_limit(pairs * len(network.vertices), limits.max_entry_evaluations, "tensor evaluations")
        check_limit(network.edges ** 2 * network.ring.size, limits.max_jet_slots, "kernel coordinate slots")
    except ModularBudgetRefusal as error:
        return DeformationResult(None, True, f"{error}; no kernel allocated or port pairs scanned")
    stats = dict(port_pairs_scanned=0, tensor_evaluations=0)
    value = _deformation_value(network, layout, t, zeta, locations, stats)
    return DeformationResult(value, False, "finite (D.8) port-pair sum", **stats, kernel_allocated=True)


@dataclass(frozen=True)
class TaperedNode:
    t: int
    valuation: int
    multiplicity: int
    precision: int


@dataclass(frozen=True)
class TaperedPlan:
    edges: int
    a: int
    u: int
    precision: int
    h: int
    degree_bound: int
    interpolation_points: int
    nodes: tuple[TaperedNode, ...]
    port_pair_bound: int
    tensor_evaluation_bound: int
    jet_coordinate_bound: int
    recovery_bit_bound: int


def plan_tapered_recovery(
    network: BooleanTensorNetwork, layout: PacketLayout, u: int,
    *, budgets: ModularBudgets | None = None,
) -> TaperedPlan:
    """Lemma D.5's exact finite integers, plus preflight resource upper bounds.

    The coordinate bound covers retained jets and interpolation/kernel arrays,
    excluding Python object overhead and arithmetic temporaries.  The rational
    bound uses a product of all recovery-weight denominators; it is deliberately
    conservative and is not the full algorithm's bit-complexity bound.
    """
    limits = budgets or ModularBudgets()
    ring, M, a = network.ring, network.edges, layout.a
    layout.locations(M)
    if not isinstance(u, int) or not 1 <= u <= a or M < u:
        raise ValueError("Lemma D.5 requires 1 <= u <= a and M >= u")
    check_ring_budget(ring, limits)
    check_limit(M, limits.max_edges, "network edges")
    check_limit(len(network.vertices), limits.max_vertices, "network vertices")
    L = bounded_power_of_two(a, limits.max_packet_slots, "packet slots")
    degree = M * (2 * L - 1)
    points = degree + 1
    check_limit(points, limits.max_interpolation_points, "interpolation points")
    if points > 1 << ring.extension_degree:
        raise ModularBudgetRefusal("residue field has too few distinct interpolation points")
    h = (3 * M + u - 1) // u
    nodes = tuple(
        TaperedNode(t, valuation2(t), h * (u - valuation2(t)),
                    ring.precision + 8 * h * (L >> valuation2(t)) + 8 * M * (a + 1))
        for t in range(-L, L) if t and valuation2(t) < u
    )
    largest_precision = max(node.precision for node in nodes)
    check_limit(largest_precision, limits.max_precision_bits, "child precision bits")
    pairs = bounded_power_of_two(
        2 * M, limits.max_boolean_vectors // (len(nodes) * points), "port pairs per evaluation",
    )
    pair_bound = len(nodes) * points * pairs
    entry_bound = pair_bound * len(network.vertices)
    check_limit(entry_bound, limits.max_entry_evaluations, "tensor evaluations")
    total_q = sum(node.multiplicity for node in nodes)
    largest_q = max(node.multiplicity for node in nodes)
    slots = ring.size * (total_q + 8 * points + 2 * largest_q + M * M + 10)
    check_limit(slots, limits.max_jet_slots, "jet coordinate slots")
    bit_bound = largest_precision + (total_q + largest_q) * (a + 2) * (total_q + 1) + total_q.bit_length() + 2
    check_limit(bit_bound, limits.max_recovery_bits, "rational recovery bits")
    return TaperedPlan(M, a, u, ring.precision, h, degree, points, nodes, pair_bound, entry_bound, slots, bit_bound)


@dataclass(frozen=True)
class TaperedJet:
    t: int
    ring: ModularRing
    coefficients: tuple[Element, ...]


def recover_tapered_jets(
    ring: ModularRing, plan: TaperedPlan, samples: Sequence[TaperedJet],
) -> Element:
    """Equation (D.9), using modular quotients and a rational final combination.

    The caller supplies jets of the same polynomial satisfying Lemma D.5's
    integrality and degree hypotheses.  Shape and precision are checked here;
    arbitrary unrelated jets cannot certify those polynomial hypotheses.
    """
    if ring.precision != plan.precision:
        raise ValueError("recovery ring precision does not match the plan")
    by_node = {sample.t: sample for sample in samples}
    if len(by_node) != len(samples) or set(by_node) != {node.t for node in plan.nodes}:
        raise ValueError("need exactly one jet at each active tapered node")
    for node in plan.nodes:
        sample = by_node[node.t]
        if (sample.ring.ell, sample.ring.polynomial) != (ring.ell, ring.polynomial):
            raise ValueError("jet coefficient basis does not match the recovery ring")
        if sample.ring.precision < node.precision:
            raise ValueError("jet precision is below the Lemma D.5 requirement")
        if len(sample.coefficients) < node.multiplicity:
            raise ValueError("missing tapered jet coefficients")
        for value in sample.coefficients:
            sample.ring._check(value)
    total = [Fraction(0) for _ in range(ring.size)]
    alpha = 1 << (plan.a + 1)
    for node in plan.nodes:
        t, q = node.t, node.multiplicity
        local = ring.at_precision(node.precision)
        modulus = local.modulus
        reciprocal = [1] + [0] * (q - 1)
        E0 = Fraction(1)
        for other in plan.nodes:
            if other.t == t:
                continue
            E0 *= Fraction(-other.t, t - other.t) ** other.multiplicity
            ratio = local.from_fraction(Fraction(alpha, t - other.t))[0]
            factor = [((-1) ** j * comb(other.multiplicity + j - 1, j) * pow(ratio, j, modulus)) % modulus
                      for j in range(q)]
            updated = [0] * q
            for j, left in enumerate(reciprocal):
                for k in range(q - j):
                    updated[j + k] += left * factor[k]
            reciprocal = [value % modulus for value in updated]
        coefficients = by_node[t].coefficients
        weight = E0
        for j in range(q):
            quotient = [0] * ring.size
            for k in range(j + 1):
                for coordinate, value in enumerate(coefficients[k]):
                    quotient[coordinate] += reciprocal[j - k] * value
            for coordinate, value in enumerate(quotient):
                total[coordinate] += weight * (value % modulus)
            weight *= Fraction(-t, alpha)
    if any(value.denominator % 2 == 0 for value in total):
        raise ArithmeticError("reconstructed value is not integral; supplied jets violate the hypotheses")
    return ring.element(tuple(value.numerator * pow(value.denominator, -1, ring.modulus) for value in total))


@dataclass(frozen=True)
class ModularRecoveryResult:
    value: Element | None
    refused: bool
    reason: str
    plan: TaperedPlan | None = None
    port_pairs_scanned: int = 0
    tensor_evaluations: int = 0
    jets_allocated: bool = False
    recursive_contraction: bool = False
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False

    @property
    def recovered(self) -> bool:
        return self.value is not None and not self.refused


def recover_network_modular(
    network: BooleanTensorNetwork, layout: PacketLayout, u: int,
    *, budgets: ModularBudgets | None = None,
) -> ModularRecoveryResult:
    """Finite D.8 evaluations, residue-field interpolation, then D.9 recovery."""
    try:
        plan = plan_tapered_recovery(network, layout, u, budgets=budgets)
    except ModularBudgetRefusal as error:
        return ModularRecoveryResult(None, True, f"{error}; no jets allocated or port pairs scanned")
    locations = layout.locations(network.edges)
    samples = []
    stats = dict(port_pairs_scanned=0, tensor_evaluations=0)
    for node in plan.nodes:
        local = network.at_precision(node.precision)
        points = tuple(local.ring.residue_point(j) for j in range(plan.interpolation_points))
        values = tuple(_deformation_value(local, layout, node.t, point, locations, stats) for point in points)
        coefficients = interpolate_polynomial(local.ring, points, values)
        samples.append(TaperedJet(node.t, local.ring, coefficients[:node.multiplicity]))
    value = recover_tapered_jets(network.ring, plan, samples)
    return ModularRecoveryResult(
        value, False, "finite Cauchy deformation and tapered Hermite recovery", plan,
        **stats, jets_allocated=True,
    )
