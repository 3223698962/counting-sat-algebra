"""Appendix D.4.2 integer schedules and finite collision packing.

The conditional expectation is computed from exact incidence histograms.
No tensor entry is read.  Packet orientations and the prospective child
topologies are checked, but no mode expansion or recursive contraction is
performed here.  In particular, finite packing is not Theorem D.1.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, fields
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import isqrt

from .modular_network import BooleanTensorNetwork, bounded_power_of_two, check_limit
from .modular_recovery import PacketLayout
from .modular_ring import ModularBudgetRefusal


STEP21_SCOPE = (
    "Finite Appendix D.4.2 integer schedules, exact collision conditional expectations, "
    "packet orientations and child topologies. The D.6 collision guard is M >= "
    "2*D_arr^2*Delta. Checked finite conclusions do not imply those hypotheses. "
    "D.4.5 modes, recursive child contractions, D.3 incidence encoding and the "
    "Theorem D.1 time/space bounds remain open. No Theorem 1.1 implementation "
    "or general 2^{o(n)} algorithm."
)


@dataclass(frozen=True)
class PackingBudgets:
    max_edges: int = 65_536
    max_vertices: int = 131_072
    max_labels: int = 4096
    max_pair_classes: int = 65_536
    max_packet_slots: int = 256
    max_label_trials: int = 5_000_000
    max_score_updates: int = 20_000_000
    max_padded_slots: int = 2_000_000
    max_packets: int = 1_000_000
    max_trace_entries: int = 4096
    max_child_edges: int = 131_072
    max_parameter_bits: int = 4096

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ValueError(f"{field.name} must be a nonnegative integer")


def _integer(value: int, minimum: int, name: str) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")


def _bounded_power(base: int, exponent: int, limit: int, name: str) -> int:
    value = 1
    for _ in range(exponent):
        if value > limit // base:
            raise ModularBudgetRefusal(f"{name} exceeds budget {limit}")
        value *= base
    return value


@dataclass(frozen=True)
class AppendixDIntegers:
    lambda_s: int
    depth: int
    u: int
    a: int
    previous_a: int
    r: int = 64
    epsilon: Fraction = Fraction(1, 128 * 1024)


def appendix_d_integers(lambda_s: int, depth: int = 0) -> AppendixDIntegers:
    """D.4 without forming s, L or D_arr; also usable for logarithmic witnesses."""
    _integer(lambda_s, 1, "lambda_s")
    _integer(depth, 0, "depth")
    initial = lambda_s // (128 * 1024)
    # Sum_j floor(t/2^j) = 2*t-popcount(t). Subtract the tail without a loop
    # over either the depth or the (potentially large) logarithmic witness.
    u = initial >> depth if depth < initial.bit_length() else 0
    tail = initial >> (depth + 1) if depth + 1 < initial.bit_length() else 0
    a = 2 * (2 * initial - initial.bit_count() - 2 * tail + tail.bit_count())
    return AppendixDIntegers(lambda_s, depth, u, a, a - 2 * u)


def _ceil_root(value: int, degree: int) -> int:
    root = 1 << ((value.bit_length() + degree - 1) // degree)
    while True:
        smaller = ((degree - 1) * root + value // root ** (degree - 1)) // degree
        if smaller >= root:
            break
        root = smaller
    return root + int(root ** degree != value)


@dataclass(frozen=True)
class PackingHypotheses:
    edges: int
    r: int
    a: int
    delta: int
    L: int
    labels: int
    D_arr: int
    collision: bool
    density: bool
    packet_degree: bool

    @property
    def admitted(self) -> bool:
        return self.collision and self.density and self.packet_degree

    @property
    def failures(self) -> tuple[str, ...]:
        return tuple(name for ok, name in (
            (self.collision, "M >= 2*D_arr^2*Delta"),
            (self.density, "M >= L*D_arr"),
            (self.packet_degree, "2*L <= Delta"),
        ) if not ok)


def packing_hypotheses(
    edges: int, r: int, a: int, delta: int, *, budgets: PackingBudgets | None = None,
) -> PackingHypotheses:
    _integer(edges, 0, "edges")
    _integer(r, 2, "r")
    _integer(a, 1, "a")
    _integer(delta, 1, "delta")
    limits = budgets or PackingBudgets()
    check_limit(max(edges.bit_length(), delta.bit_length()), limits.max_parameter_bits, "parameter bits")
    L = bounded_power_of_two(a, limits.max_packet_slots, "packet slots")
    labels = _bounded_power(r, a, limits.max_labels, "vertex labels")
    D_arr = labels * labels
    check_limit(D_arr, limits.max_pair_classes, "pair-label classes")
    return PackingHypotheses(edges, r, a, delta, L, labels, D_arr,
                             edges >= 2 * D_arr * D_arr * delta,
                             edges >= L * D_arr, 2 * L <= delta)


@dataclass(frozen=True)
class AppendixDProposal:
    s: int
    edges: int
    depth: int
    integers: AppendixDIntegers
    H: int
    delta: int
    leaf: bool
    hypotheses: PackingHypotheses | None
    reason: str


def appendix_d_proposal(
    s: int, edges: int, depth: int = 0, *, budgets: PackingBudgets | None = None,
) -> AppendixDProposal:
    """H=ceil(s/floor(sqrt(lambda))) and Delta=ceil(s^(9/10)), exactly.

    This describes a schedule; it neither executes the exceptional-size leaf
    fallback nor asserts the graph/entry hypotheses of Theorem D.1.
    """
    _integer(s, 2, "s")
    _integer(edges, 0, "edges")
    limits = budgets or PackingBudgets()
    check_limit(9 * s.bit_length(), limits.max_parameter_bits, "root parameter intermediate bits")
    integers = appendix_d_integers((s - 1).bit_length(), depth)
    divisor = isqrt(integers.lambda_s)
    H = (s + divisor - 1) // divisor
    delta = _ceil_root(s ** 9, 10)
    if edges <= H:
        return AppendixDProposal(s, edges, depth, integers, H, delta, True, None, "direct-leaf threshold")
    if integers.u < 1 or integers.u > edges:
        return AppendixDProposal(s, edges, depth, integers, H, delta, False, None,
                                 "active-call guard 1 <= u <= M fails; L and D_arr are not formed")
    hypotheses = packing_hypotheses(edges, integers.r, integers.a, delta, budgets=limits)
    return AppendixDProposal(s, edges, depth, integers, H, delta, False, hypotheses,
                             "D.6 admitted" if hypotheses.admitted else "D.6 numerical hypotheses fail")


@dataclass(frozen=True)
class _Incidences:
    endpoints: tuple[tuple[int, int], ...]
    neighbors: tuple[dict[int, int], ...]
    degrees: tuple[int, ...]


def _incidences(network: BooleanTensorNetwork) -> _Incidences:
    ends = [[-1, -1] for _ in range(network.edges)]
    for index, vertex in enumerate(network.vertices):
        for edge in vertex.ports:
            ends[edge][vertex.side] = index
    neighbors: list[dict[int, int]] = [{} for _ in network.vertices]
    for row, column in ends:
        neighbors[row][column] = neighbors[row].get(column, 0) + 1
        neighbors[column][row] = neighbors[column].get(row, 0) + 1
    return _Incidences(tuple((row, column) for row, column in ends), tuple(neighbors),
                       tuple(len(vertex.ports) for vertex in network.vertices))


def _initial_numerator(data: _Incidences, network: BooleanTensorNetwork, K: int) -> int:
    parallel_pairs = sum(m * m for i, vertex in enumerate(network.vertices) if vertex.side == 0
                         for m in data.neighbors[i].values())
    shared_only = sum(degree * degree for degree in data.degrees) - 2 * parallel_pairs
    # This is K^2 E[sum_gamma (X_gamma-M/K^2)^2]. Ordered pairs are used.
    return parallel_pairs * (K * K - 1) + shared_only * (K - 1)


@dataclass(frozen=True)
class CollisionChoice:
    vertex: int
    label: int
    before: Fraction
    candidates: tuple[Fraction, ...]
    after: Fraction


@dataclass(frozen=True)
class CollisionLabels:
    labels: tuple[int, ...]
    initial_potential: Fraction
    final_potential: Fraction
    pair_loads: tuple[tuple[int, ...], ...]
    choices: tuple[CollisionChoice, ...]
    label_trials: int
    score_updates: int


def _label_plan(network, K, data, limits, trace):
    trials = K * len(network.vertices)
    # Each edge contributes to an assigned-neighbor histogram at most once.
    updates = trials + K * network.edges + sum(degree * degree for degree in data.degrees)
    check_limit(trials, limits.max_label_trials, "label trials")
    check_limit(updates, limits.max_score_updates, "score updates")
    if trace:
        check_limit(trials, limits.max_trace_entries, "conditional-expectation trace entries")
    return trials, updates


def _choose_labels(network, K, data, trace):
    """Exact CE update, using completed loads and assigned-neighbor histograms.

    For a row vertex v with degree d and z unassigned incident ports, write
    C_ab for completed pair loads, R_a for row-only loads, and d_b for v's
    assigned-neighbor histogram. Up to a candidate-independent constant the
    score for a is

      K sum_b C_ab d_b + z sum_b C_ab + d R_a
        + (K-1) sum_{unassigned neighbors w} m_vw d_w(a).

    The last term handles shared endpoints, including parallel incidences.
    If c_a are these scores, the change in K^2 E[Phi] is
    2(K c_a - sum_b c_b). Column scores are the transposed formula.
    """
    labels: list[int | None] = [None] * len(network.vertices)
    histograms: list[dict[int, int]] = [{} for _ in network.vertices]
    completed = [[0] * K for _ in range(K)]
    totals = [[0] * K, [0] * K]
    pending = [[0] * K, [0] * K]
    numerator = initial = _initial_numerator(data, network, K)
    choices = []
    trials = updates = 0
    for vertex, tensor in enumerate(network.vertices):
        side = tensor.side
        histogram = histograms[vertex]
        degree = data.degrees[vertex]
        unassigned = degree - sum(histogram.values())
        scores = [unassigned * totals[side][label] + degree * pending[side][label] for label in range(K)]
        updates += K
        for other_label, multiplicity in histogram.items():
            for label in range(K):
                load = completed[label][other_label] if side == 0 else completed[other_label][label]
                scores[label] += K * load * multiplicity
                updates += 1
        for neighbor, multiplicity in data.neighbors[vertex].items():
            if labels[neighbor] is None:
                for label, weight in histograms[neighbor].items():
                    scores[label] += (K - 1) * multiplicity * weight
                    updates += 1
        chosen = min(range(K), key=scores.__getitem__)
        total_scores = sum(scores)
        after = numerator + 2 * (K * scores[chosen] - total_scores)
        if after > numerator or after < 0:
            raise ArithmeticError("collision conditional expectation increased or became negative")
        if trace:
            choices.append(CollisionChoice(
                vertex, chosen, Fraction(numerator, K * K),
                tuple(Fraction(numerator + 2 * (K * score - total_scores), K * K) for score in scores),
                Fraction(after, K * K),
            ))
        numerator = after
        trials += K
        labels[vertex] = chosen
        for neighbor, multiplicity in data.neighbors[vertex].items():
            other_label = labels[neighbor]
            if other_label is None:
                pending[side][chosen] += multiplicity
                histograms[neighbor][chosen] = histograms[neighbor].get(chosen, 0) + multiplicity
            else:
                row, column = (chosen, other_label) if side == 0 else (other_label, chosen)
                completed[row][column] += multiplicity
                totals[0][row] += multiplicity
                totals[1][column] += multiplicity
                pending[1 - side][other_label] -= multiplicity
        histograms[vertex].clear()
    final = Fraction(sum(load * load for row in completed for load in row), 1) - Fraction(network.edges ** 2, K * K)
    if final != Fraction(numerator, K * K):
        raise ArithmeticError("final collision loads disagree with conditional expectation")
    return CollisionLabels(tuple(int(label) for label in labels), Fraction(initial, K * K), final,
                           tuple(tuple(row) for row in completed), tuple(choices), trials, updates)


def choose_collision_labels(
    network: BooleanTensorNetwork, r: int, a: int, *, trace: bool = False,
    budgets: PackingBudgets | None = None,
) -> CollisionLabels:
    """Choose labels by successive minimization of the exact conditional potential."""
    limits = budgets or PackingBudgets()
    check_limit(network.edges, limits.max_edges, "network edges")
    check_limit(len(network.vertices), limits.max_vertices, "network vertices")
    parameters = packing_hypotheses(network.edges, r, a, max(1, network.max_degree), budgets=limits)
    data = _incidences(network)
    _label_plan(network, parameters.labels, data, limits, trace)
    return _choose_labels(network, parameters.labels, data, trace)


@dataclass(frozen=True)
class PackingPlan:
    hypotheses: PackingHypotheses
    initial_potential: Fraction
    label_trial_bound: int
    score_update_bound: int
    padded_slot_bound: int
    packet_bound: int


@dataclass(frozen=True)
class PackedNetwork:
    source: BooleanTensorNetwork
    hypotheses: PackingHypotheses
    labeling: CollisionLabels
    layout: PacketLayout
    packet_symbols: tuple[tuple[tuple[int, int], ...], ...]
    H_pad: int

    def symbol(self, vertex: int, coordinate: int) -> int:
        return self.labeling.labels[vertex] // (self.hypotheses.r ** coordinate) % self.hypotheses.r


@dataclass(frozen=True)
class PackingCheck:
    problems: tuple[str, ...]
    row_group_ports: tuple[tuple[int, ...], ...]
    column_group_ports: tuple[tuple[int, ...], ...]
    max_child_edges: int
    max_packet_degree: int
    child_shrink: bool

    @property
    def ok(self) -> bool:
        return not self.problems


def check_packing(packed: PackedNetwork) -> PackingCheck:
    """Check the actual finite certificate, separately from the D.6 hypotheses."""
    network, h = packed.source, packed.hypotheses
    empty = PackingCheck(("invalid packet certificate",), (), (), 0, 0, False)
    if h.edges != network.edges or h.L != 1 << h.a or h.labels != h.r ** h.a or h.D_arr != h.labels ** 2:
        return empty
    if (h.collision, h.density, h.packet_degree) != (
        h.edges >= 2 * h.D_arr * h.D_arr * h.delta, h.edges >= h.L * h.D_arr, 2 * h.L <= h.delta,
    ):
        return PackingCheck(("D.6 flags disagree with their exact inequalities",), (), (), 0, 0, False)
    if packed.layout.a != h.a or len(packed.packet_symbols) != len(packed.layout.packets):
        return empty
    if len(packed.labeling.labels) != len(network.vertices) or any(not 0 <= label < h.labels for label in packed.labeling.labels):
        return empty
    try:
        packed.layout.locations(network.edges)
    except ValueError as error:
        return PackingCheck((str(error),), (), (), 0, 0, False)
    data = _incidences(network)
    problems = []
    row_counts = [[0] * h.r for _ in range(h.a)]
    col_counts = [[0] * h.r for _ in range(h.a)]
    for vertex, tensor in enumerate(network.vertices):
        counts = row_counts if tensor.side == 0 else col_counts
        for coordinate in range(h.a):
            counts[coordinate][packed.symbol(vertex, coordinate)] += len(tensor.ports)
    if any(count * h.r > 3 * network.edges for group in row_counts + col_counts for count in group):
        problems.append("real group ports exceed 3M/r")
    if len(packed.layout.packets) * h.L > 3 * network.edges:
        problems.append("packet count exceeds 3M/L")
    if packed.H_pad < 0 or packed.H_pad % h.L or len(packed.layout.packets) * h.L != packed.H_pad * h.D_arr:
        problems.append("padding does not match the uniform pair-class construction")
    actual_loads = Counter((packed.labeling.labels[row], packed.labeling.labels[column]) for row, column in data.endpoints)
    energy = Fraction(sum(load * load for load in actual_loads.values())) - Fraction(network.edges ** 2, h.D_arr)
    initial = Fraction(_initial_numerator(data, network, h.labels), h.D_arr)
    if (initial != packed.labeling.initial_potential or energy != packed.labeling.final_potential
            or energy > initial or initial > 2 * network.edges * h.delta):
        problems.append("collision potential certificate fails")
    loads = packed.labeling.pair_loads
    if (len(loads) != h.labels or any(len(row) != h.labels for row in loads)
            or any(loads[row][column] != actual_loads.get((row, column), 0)
                   for row in range(h.labels) for column in range(h.labels))):
        problems.append("pair-load table disagrees with the labeled graph")
    if max(actual_loads.values(), default=0) > packed.H_pad:
        problems.append("a pair-label class exceeds its padding")
    max_packet_degree = 0
    for packet, symbols in zip(packed.layout.packets, packed.packet_symbols):
        if len(symbols) != h.a or any(len(pair) != 2 or not 0 <= pair[0] <= pair[1] < h.r for pair in symbols):
            problems.append("invalid packet symbol pairs")
            continue
        degrees = [Counter() for _ in range(h.a)]
        for slot, edge in enumerate(packet):
            if edge is None:
                continue
            row, column = data.endpoints[edge]
            for coordinate, (low, high) in enumerate(symbols):
                expected = (high, low) if slot >> coordinate & 1 else (low, high)
                actual = (packed.symbol(row, coordinate), packed.symbol(column, coordinate))
                if actual != expected:
                    problems.append("packet slot orientation disagrees with endpoint labels")
                degrees[coordinate].update(actual)
        max_packet_degree = max(max_packet_degree, max((max(counts.values(), default=0) for counts in degrees), default=0))
    if max(network.max_degree, max_packet_degree) > h.delta:
        problems.append("prospective child degree exceeds Delta")
    child_edges = max((row_counts[i][g] + col_counts[i][g] for i in range(h.a) for g in range(h.r)), default=0)
    if child_edges * h.r > 6 * network.edges:
        problems.append("prospective child edges exceed 6M/r")
    return PackingCheck(tuple(dict.fromkeys(problems)), tuple(map(tuple, row_counts)), tuple(map(tuple, col_counts)),
                        child_edges, max_packet_degree, 8 * child_edges <= network.edges)


@dataclass(frozen=True)
class PackingResult:
    packing: PackedNetwork | None
    refused: bool
    reason: str
    plan: PackingPlan | None = None
    check: PackingCheck | None = None
    label_trials: int = 0
    score_updates: int = 0
    labels_allocated: bool = False
    packets_allocated: bool = False
    tensor_evaluations: int = 0
    recursive_contraction: bool = False
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False

    @property
    def packed(self) -> bool:
        return self.packing is not None and not self.refused


def _make_packets(network, h, data, labeling):
    max_load = max((load for row in labeling.pair_loads for load in row), default=0)
    H_pad = (max_load + h.L - 1) // h.L * h.L
    if H_pad == 0:
        return PackedNetwork(network, h, labeling, PacketLayout(h.a, ()), (), 0)
    classes: dict[tuple[int, int], list[int]] = {}
    for edge, (row, column) in enumerate(data.endpoints):
        key = (labeling.labels[row], labeling.labels[column])
        classes.setdefault(key, []).append(edge)
    packets, all_symbols = [], []
    pairs = tuple(combinations_with_replacement(range(h.r), 2))
    powers = tuple(h.r ** i for i in range(h.a))
    for symbols in product(pairs, repeat=h.a):
        equal = tuple(i for i, (low, high) in enumerate(symbols) if low == high)
        copies = 1 << len(equal)
        slot_data = []
        for slot in range(h.L):
            row = sum((high if slot >> i & 1 else low) * powers[i] for i, (low, high) in enumerate(symbols))
            column = sum((low if slot >> i & 1 else high) * powers[i] for i, (low, high) in enumerate(symbols))
            copy = sum(((slot >> i) & 1) << j for j, i in enumerate(equal))
            slot_data.append((classes.get((row, column), ()), copy))
        for block in range(H_pad // copies):
            packet = tuple(edges[block * copies + copy] if block * copies + copy < len(edges) else None
                           for edges, copy in slot_data)
            packets.append(packet)
            all_symbols.append(symbols)
    return PackedNetwork(network, h, labeling, PacketLayout(h.a, tuple(packets)), tuple(all_symbols), H_pad)


def pack_network(
    network: BooleanTensorNetwork, r: int, a: int, *, delta: int | None = None,
    require_hypotheses: bool = True, trace: bool = False, budgets: PackingBudgets | None = None,
) -> PackingResult:
    """Label and pack a network, returning no certificate on any refusal.

    Numerical D.6 admission is the default. ``require_hypotheses=False`` is
    an explicit finite experiment: the actual packet and group conclusions
    must still pass. Every resource budget is checked before labeling starts.
    """
    limits = budgets or PackingBudgets()
    try:
        check_limit(network.edges, limits.max_edges, "network edges")
        check_limit(len(network.vertices), limits.max_vertices, "network vertices")
        _integer(a, 1, "a")
        L = bounded_power_of_two(a, limits.max_packet_slots, "packet slots")
        bound = max(network.max_degree, 2 * L) if delta is None else delta
        _integer(bound, 1, "delta")
        if bound < network.max_degree:
            raise ValueError("Delta is below the actual vertex degree")
        h = packing_hypotheses(network.edges, r, a, bound, budgets=limits)
        if require_hypotheses and not h.admitted:
            return PackingResult(None, True, "D.6 requires " + "; ".join(h.failures) + "; labeling not started")
        data = _incidences(network)
        trials, updates = _label_plan(network, h.labels, data, limits, trace)
        initial = Fraction(_initial_numerator(data, network, h.labels), h.D_arr)
        ceiling = (initial.numerator + initial.denominator - 1) // initial.denominator
        root = isqrt(ceiling)
        root += int(root * root != ceiling)
        load_bound = min(network.edges, (network.edges + h.D_arr - 1) // h.D_arr + root)
        pad_bound = (load_bound + h.L - 1) // h.L * h.L
        slots = pad_bound * h.D_arr
        packets = slots // h.L
        check_limit(slots, limits.max_padded_slots, "padded slots")
        check_limit(packets, limits.max_packets, "packets")
        plan = PackingPlan(h, initial, trials, updates, slots, packets)
    except ModularBudgetRefusal as error:
        return PackingResult(None, True, f"{error}; labeling and packet allocation not started")
    labeling = _choose_labels(network, h.labels, data, trace)
    packed = _make_packets(network, h, data, labeling)
    checked = check_packing(packed)
    if (len(packed.layout.packets) > plan.packet_bound or labeling.score_updates > plan.score_update_bound
            or labeling.label_trials > plan.label_trial_bound):
        raise ArithmeticError("packing exceeded its preflight upper bounds")
    if not checked.ok:
        return PackingResult(None, True, "; ".join(checked.problems), plan, checked,
                             labeling.label_trials, labeling.score_updates, True, True)
    return PackingResult(packed, False, "finite packet certificate checked", plan, checked,
                         labeling.label_trials, labeling.score_updates, True, True)


@dataclass(frozen=True)
class ChildTopology:
    coordinate: int
    group: int
    original_vertices: tuple[int, ...]
    original_ports: tuple[tuple[int, ...], ...]
    packet_ports: tuple[tuple[int, ...], ...]
    port_origins: tuple[tuple[int, int], ...]

    @property
    def edges(self) -> int:
        return len(self.port_origins)

    @property
    def max_degree(self) -> int:
        return max(map(len, self.original_ports + self.packet_ports), default=0)


def child_topology(
    packed: PackedNetwork, coordinate: int, group: int, *, budgets: PackingBudgets | None = None,
) -> ChildTopology:
    """D.4.6 port geometry only; no packet tensor or contraction is evaluated."""
    h, network = packed.hypotheses, packed.source
    _integer(coordinate, 0, "coordinate")
    _integer(group, 0, "group")
    if coordinate >= h.a or group >= h.r:
        raise ValueError("coordinate or group is outside the label alphabet")
    checked = check_packing(packed)
    if not checked.ok:
        raise ValueError("invalid packing: " + "; ".join(checked.problems))
    limits = budgets or PackingBudgets()
    edge_count = checked.row_group_ports[coordinate][group] + checked.column_group_ports[coordinate][group]
    check_limit(edge_count, limits.max_child_edges, "child edges")
    check_limit(len(packed.layout.packets), limits.max_packets, "child packet vertices")
    vertex_count = sum(packed.symbol(vertex, coordinate) == group for vertex in range(len(network.vertices)))
    check_limit(vertex_count + len(packed.layout.packets), limits.max_vertices, "child vertices")
    locations = packed.layout.locations(network.edges)
    vertices, ordered, origins = [], [], []
    packet_ports: list[list[int]] = [[] for _ in packed.layout.packets]
    for vertex, tensor in enumerate(network.vertices):
        if packed.symbol(vertex, coordinate) != group:
            continue
        vertices.append(vertex)
        ports = []
        for original_edge in tensor.ports:
            edge = len(origins)
            origins.append((vertex, original_edge))
            ports.append(edge)
            packet_ports[locations[original_edge][0]].append(edge)
        ordered.append(tuple(ports))
    child = ChildTopology(coordinate, group, tuple(vertices), tuple(ordered), tuple(map(tuple, packet_ports)), tuple(origins))
    if child.edges != edge_count or child.max_degree > h.delta:
        raise ArithmeticError("child topology violated its checked edge or degree bound")
    return child
