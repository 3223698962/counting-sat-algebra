"""Finite Boolean tensor networks and the Appendix D.2 CNF encoding.

Contraction here is the direct, streamed leaf routine of Appendix D.1.
It is exponential in the number of edges, with a preflight budget.  The
mode expansion and recursive child contractions from D.4.5–D.4.7 are not
implemented. Finite collision layouts are constructed in modular_packing.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Callable, Sequence

from .cnf import CNF
from .modular_ring import Element, ModularBudgetRefusal, ModularRing


STEP20_SCOPE = (
    "Finite Appendix D.1–D.2 networks and D.4.1, D.4.3–D.4.4 local arithmetic "
    "and recovery. Direct edge and port-pair enumeration is budgeted. "
    "Step 20 itself does not construct collision labels or packet topology. "
    "Weighted modes, recursive child contractions, and Lemma D.2's "
    "sparse-system encoding remain open. No Theorem D.1 time/space bound or "
    "Theorem 1.1 implementation; no general 2^{o(n)} algorithm."
)


@dataclass(frozen=True)
class ModularBudgets:
    max_edges: int = 24
    max_vertices: int = 100_000
    max_boolean_vectors: int = 262_144
    max_entry_evaluations: int = 2_000_000
    max_precision_bits: int = 4096
    max_ring_coordinates: int = 4096
    max_field_degree: int = 32
    max_packet_slots: int = 16
    max_interpolation_points: int = 128
    max_jet_slots: int = 200_000
    max_recovery_bits: int = 2_000_000

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ValueError(f"{field.name} must be a nonnegative integer")


def check_limit(value: int, limit: int, name: str) -> None:
    if value > limit:
        raise ModularBudgetRefusal(f"{name} {value} exceeds budget {limit}")


def bounded_power_of_two(exponent: int, limit: int, name: str) -> int:
    if limit < 1 or exponent >= limit.bit_length():
        raise ModularBudgetRefusal(f"{name} 2^{exponent} exceeds budget {limit}")
    return 1 << exponent


def check_ring_budget(ring: ModularRing, budgets: ModularBudgets) -> None:
    check_limit(ring.precision, budgets.max_precision_bits, "precision bits")
    check_limit(ring.extension_degree, budgets.max_field_degree, "extension degree")
    check_limit(ring.size, budgets.max_ring_coordinates, "ring coordinates")


@dataclass(frozen=True)
class TensorVertex:
    """An immutable description evaluated at its own creation precision.

    ``entry`` must be a deterministic function of the ordered Boolean ports,
    returning V-major coordinates.  Its captured description must be immutable.
    It is never passed the later precision of the calling network.
    """

    name: str
    side: int
    ports: tuple[int, ...]
    creation_ring: ModularRing
    entry: Callable[[tuple[int, ...]], Sequence[int]]

    def __post_init__(self) -> None:
        object.__setattr__(self, "ports", tuple(self.ports))
        if self.side not in (0, 1):
            raise ValueError("vertex side must be zero or one")
        if any(not isinstance(port, int) or port < 0 for port in self.ports):
            raise ValueError("port indices must be nonnegative integers")
        if len(set(self.ports)) != len(self.ports):
            raise ValueError("parallel edges need distinct ports")
        if not callable(self.entry):
            raise ValueError("tensor entry must be callable")

    def evaluate(self, bits: tuple[int, ...], ring: ModularRing) -> Element:
        if len(bits) != len(self.ports) or any(bit not in (0, 1) for bit in bits):
            raise ValueError("tensor entry requires one bit per ordered port")
        original = self.creation_ring.element(self.entry(bits))
        return self.creation_ring.canonical_lift(original, ring)


@dataclass(frozen=True)
class BooleanTensorNetwork:
    ring: ModularRing
    edges: int
    vertices: tuple[TensorVertex, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "vertices", tuple(self.vertices))
        if not isinstance(self.edges, int) or self.edges < 0:
            raise ValueError("edge count must be a nonnegative integer")
        incidences: dict[int, list[int]] = {}
        names = set()
        for vertex in self.vertices:
            if vertex.name in names:
                raise ValueError("vertex names must be distinct")
            names.add(vertex.name)
            source = vertex.creation_ring
            if (source.ell, source.polynomial) != (self.ring.ell, self.ring.polynomial):
                raise ValueError("all tensors must use the network coefficient basis")
            if source.precision > self.ring.precision:
                raise ValueError("network precision is below a tensor's creation precision")
            for edge in vertex.ports:
                if edge >= self.edges:
                    raise ValueError("port is outside the network edge range")
                incidences.setdefault(edge, []).append(vertex.side)
        if len(incidences) != self.edges or any(sorted(sides) != [0, 1] for sides in incidences.values()):
            raise ValueError("every edge must have exactly one endpoint on each side")

    @property
    def max_degree(self) -> int:
        return max((len(vertex.ports) for vertex in self.vertices), default=0)

    def at_precision(self, precision: int) -> BooleanTensorNetwork:
        return BooleanTensorNetwork(self.ring.at_precision(precision), self.edges, self.vertices)


@dataclass(frozen=True)
class NetworkContractionResult:
    value: Element | None
    refused: bool
    reason: str
    boolean_vectors_scanned: int = 0
    tensor_evaluations: int = 0
    vector_bound: int | None = None

    @property
    def contracted(self) -> bool:
        return self.value is not None and not self.refused


def _network_budget(network: BooleanTensorNetwork, budgets: ModularBudgets) -> int:
    check_ring_budget(network.ring, budgets)
    check_limit(network.edges, budgets.max_edges, "network edges")
    check_limit(len(network.vertices), budgets.max_vertices, "network vertices")
    vectors = bounded_power_of_two(network.edges, budgets.max_boolean_vectors, "edge assignments")
    check_limit(vectors * len(network.vertices), budgets.max_entry_evaluations, "tensor evaluations")
    return vectors


def _side_product(network: BooleanTensorNetwork, bits: tuple[int, ...], side: int, stats: dict) -> Element:
    value = network.ring.one
    for vertex in network.vertices:
        if vertex.side == side:
            local = tuple(bits[edge] for edge in vertex.ports)
            value = network.ring.mul(value, vertex.evaluate(local, network.ring))
            stats['tensor_evaluations'] += 1
    return value


def contract_network(
    network: BooleanTensorNetwork, *, budgets: ModularBudgets | None = None,
) -> NetworkContractionResult:
    """Stream the leaf sum over Boolean edges; refuse before entry evaluation."""
    limits = budgets or ModularBudgets()
    try:
        vectors = _network_budget(network, limits)
    except ModularBudgetRefusal as error:
        return NetworkContractionResult(None, True, f"{error}; edge scan not started")
    total = network.ring.zero
    stats = dict(boolean_vectors_scanned=0, tensor_evaluations=0)
    for mask in range(vectors):
        bits = tuple((mask >> edge) & 1 for edge in range(network.edges))
        product = network.ring.mul(_side_product(network, bits, 0, stats), _side_product(network, bits, 1, stats))
        total = network.ring.add(total, product)
        stats['boolean_vectors_scanned'] += 1
    return NetworkContractionResult(
        total, False, "direct Appendix D leaf sum", **stats, vector_bound=vectors,
    )


@dataclass(frozen=True)
class CnfNetworkEncoding:
    network: BooleanTensorNetwork | None
    free_variables: int
    clause_ports: int
    edges: int
    vertices: int
    refused: bool = False
    reason: str = ""
    constant_count: int | None = None


@dataclass(frozen=True)
class ModularCountResult:
    count: int | None
    refused: bool
    reason: str
    precision: int
    edges: int = 0
    vertices: int = 0
    free_variables: int = 0
    max_degree: int = 0
    boolean_vectors_scanned: int = 0
    tensor_evaluations: int = 0
    network_allocated: bool = False
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False

    @property
    def counted(self) -> bool:
        return self.count is not None and not self.refused


def _equality(bits: tuple[int, ...]) -> tuple[int]:
    return (int(not bits or all(bit == bits[0] for bit in bits)),)


def _clause_entry(negations: tuple[bool, ...]):
    def entry(bits: tuple[int, ...]) -> tuple[int]:
        return (int(any(bit != neg for bit, neg in zip(bits, negations))),)
    return entry


def _encode_cnf(cnf: CNF, limits: ModularBudgets, *, for_count: bool) -> CnfNetworkEncoding:
    clauses = []
    occurrences: dict[int, list[tuple[int, int]]] = {}
    for raw in cnf.clauses:
        signs = {}
        tautology = False
        for variable, negation in raw:
            if variable in signs and signs[variable] != negation:
                tautology = True
            signs[variable] = negation
        if tautology:
            continue
        clause = tuple(sorted(signs.items()))
        if not clause:
            return CnfNetworkEncoding(None, 0, 0, 0, 0, constant_count=0, reason="empty clause")
        index = len(clauses)
        clauses.append(clause)
        for port, (variable, _) in enumerate(clause):
            occurrences.setdefault(variable, []).append((index, port))
    free = cnf.n - len(occurrences)
    ports = sum(map(len, clauses))
    old_edges = ports + sum(max(0, len(hits) - 3) for hits in occurrences.values())
    old_vertices = len(clauses) + sum(max(1, len(hits) - 2) for hits in occurrences.values())
    edges, vertices = 2 * old_edges, old_vertices + old_edges
    try:
        check_limit(cnf.n + 1, limits.max_precision_bits, "precision bits")
        check_limit(1, limits.max_ring_coordinates, "ring coordinates")
        check_limit(1, limits.max_field_degree, "extension degree")
        check_limit(edges, limits.max_edges, "network edges")
        check_limit(vertices, limits.max_vertices, "network vertices")
        if for_count:
            vectors = bounded_power_of_two(edges, limits.max_boolean_vectors, "edge assignments")
            check_limit(vectors * vertices, limits.max_entry_evaluations, "tensor evaluations")
    except ModularBudgetRefusal as error:
        return CnfNetworkEncoding(None, free, ports, edges, vertices, True, f"{error}; no network allocated")

    ring = ModularRing(cnf.n + 1)
    raw_ports = [[-1] * len(clause) for clause in clauses]
    entries = [_clause_entry(tuple(neg for _, neg in clause)) for clause in clauses]
    endpoints: list[tuple[int, int]] = []

    def equality_vertex() -> int:
        index = len(raw_ports)
        raw_ports.append([])
        entries.append(_equality)
        return index

    def connect(left: tuple[int, int | None], right: tuple[int, int | None]) -> None:
        edge = len(endpoints)
        endpoints.append((left[0], right[0]))
        for vertex, port in (left, right):
            if port is None:
                raw_ports[vertex].append(edge)
            else:
                raw_ports[vertex][port] = edge

    for hits in occurrences.values():
        degree = len(hits)
        chain = [equality_vertex() for _ in range(max(1, degree - 2))]
        if degree <= 3:
            for hit in hits:
                connect(hit, (chain[0], None))
        else:
            for hit in hits[:2]:
                connect(hit, (chain[0], None))
            for j in range(1, len(chain) - 1):
                connect(hits[j + 1], (chain[j], None))
            for hit in hits[-2:]:
                connect(hit, (chain[-1], None))
            for left, right in zip(chain, chain[1:]):
                connect((left, None), (right, None))

    tensors = []
    for vertex, (ordered, entry) in enumerate(zip(raw_ports, entries)):
        subdivided = tuple(2 * edge + int(endpoints[edge][0] != vertex) for edge in ordered)
        tensors.append(TensorVertex(f"original-{vertex}", 0, subdivided, ring, entry))
    for edge in range(old_edges):
        tensors.append(TensorVertex(f"subdivision-{edge}", 1, (2 * edge, 2 * edge + 1), ring, _equality))
    network = BooleanTensorNetwork(ring, edges, tuple(tensors))
    if len(network.vertices) != vertices or len(endpoints) != old_edges:
        raise ArithmeticError("CNF network construction disagreed with its plan")
    return CnfNetworkEncoding(network, free, ports, edges, vertices)


def encode_cnf_network(cnf: CNF, *, budgets: ModularBudgets | None = None) -> CnfNetworkEncoding:
    """D.2: trivalent equality trees, then subdivide every wire exactly once."""
    return _encode_cnf(cnf, budgets or ModularBudgets(), for_count=False)


def count_cnf_modular(cnf: CNF, *, budgets: ModularBudgets | None = None) -> ModularCountResult:
    """Count a CNF via the D.2 network and the budgeted modular leaf routine."""
    limits = budgets or ModularBudgets()
    encoded = _encode_cnf(cnf, limits, for_count=True)
    metadata = dict(
        precision=cnf.n + 1, edges=encoded.edges, vertices=encoded.vertices,
        free_variables=encoded.free_variables,
    )
    if encoded.constant_count is not None:
        return ModularCountResult(encoded.constant_count, False, encoded.reason, **metadata)
    if encoded.refused:
        return ModularCountResult(None, True, encoded.reason, **metadata)
    network = encoded.network
    if network is None:
        raise ArithmeticError("successful encoding has no network")
    result = contract_network(network, budgets=limits)
    if result.refused or result.value is None:
        raise ArithmeticError("accepted CNF plan failed the identical leaf preflight")
    residue = (result.value[0] << encoded.free_variables) % network.ring.modulus
    if residue > 1 << cnf.n:
        raise ArithmeticError("CNF network residue is outside the possible count range")
    return ModularCountResult(
        residue, False, "D.2 encoding with direct modular leaf contraction", **metadata,
        max_degree=network.max_degree, boolean_vectors_scanned=result.boolean_vectors_scanned,
        tensor_evaluations=result.tensor_evaluations, network_allocated=True,
    )
