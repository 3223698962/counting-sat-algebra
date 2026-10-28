"""Finite Appendix D.4.5 weighted modes for one packet.

The mode table, its exponent tuples, and one crossing coefficient are
budgeted. A refusal returns no mode table, no exponent list, and no count.
Packet-tensor evaluation and recursive child contraction stay out of this
module. The weighted-support gate below copies an inequality shape from a
different algorithm layer; it is not that layer's polynomial identity.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from fractions import Fraction
from typing import Sequence

from .modular_network import BooleanTensorNetwork, bounded_power_of_two, check_limit
from .modular_recovery import (
    PacketLayout,
    cauchy_local_polynomial,
    recover_network_modular,
    valuation2,
)
from .modular_ring import (
    ModularBudgetRefusal,
    ModularRing,
    evaluate_polynomial,
)


STEP22_SCOPE = (
    "Finite Appendix D.4.5 weighted modes for one explicit packet. "
    "Mode tables, geometric crossing residues and exponent tuples are "
    "admitted only inside an explicit budget. A refusal returns no table, "
    "no tuple list and no count. D.4.6 packet-tensor coefficients, D.4.7 "
    "recursive contraction, D.3 incidence encoding and the Theorem D.1 "
    "time/space bounds remain open. No Theorem 1.1 implementation and no "
    "general 2^{o(n)} algorithm."
)


@dataclass(frozen=True)
class ModeBudgets:
    max_modes: int = 4096
    max_tuples: int = 65_536
    max_partial_choices: int = 1_000_000
    max_precision_bits: int = 64
    max_packet_slots: int = 16
    max_real_columns: int = 64
    max_support_records: int = 64

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ValueError(f"{field.name} must be a nonnegative integer")


def _integer(value: int, minimum: int, name: str) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")


@dataclass(frozen=True)
class ModeRequest:
    """One packet, one active node, and the occupied column slots.

    ``row_symbol`` and ``column_symbol`` are the endpoint symbols when bit
    ``v2(t)`` of the slot is clear. They differ exactly on a crossing packet.
    """

    a: int
    t: int
    edges: int
    precision: int
    row_symbol: int
    column_symbol: int
    real_columns: tuple[int, ...]
    packet: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "real_columns", tuple(self.real_columns))


@dataclass(frozen=True)
class WeightedMode:
    packet: int
    coordinate: int
    tile_origin: int
    weight: int
    source_group: int
    target_group: int
    real_columns: tuple[int, ...]


@dataclass(frozen=True)
class FrozenModeDescription:
    """Immutable mode data. This is not a D.4.6 packet tensor.

    Upstream packet records also carry labels, arity, degree, a cache, rows
    and residue origins. Those fields are not constructed here.
    """

    packet: int
    coordinate: int
    row_symbol: int
    column_symbol: int
    creation_precision: int
    exponent_budget: int
    modes: tuple[WeightedMode, ...]
    tensor_evaluated: bool = False

    def __post_init__(self) -> None:
        if self.tensor_evaluated:
            raise ValueError("a frozen mode description does not evaluate a packet tensor")


@dataclass(frozen=True)
class SupportRecord:
    """Declared degrees for the weighted-support refusal gate.

    The inequality ``deg_x <= deg_u + deg_z`` is the shape of an admissibility
    predicate from a different formalization. It does not encode Appendix D
    packet modes, and satisfying it is not a polynomial identity.
    """

    deg_u: int
    deg_z: int
    deg_x: int


@dataclass(frozen=True)
class ModeWitness:
    refused: bool
    reason: str
    count: int | None = None
    base_value: tuple[int, ...] | None = None
    description: FrozenModeDescription | None = None
    tuples: tuple[tuple[int, ...], ...] | None = None
    partial_choices: int = 0
    modes_allocated: bool = False
    tuples_allocated: bool = False
    support_stored: bool = False
    recovery_called: bool = False
    port_pairs_scanned: int = 0
    tensor_evaluations: int = 0
    recovery_tensor_evaluations: int = 0
    coefficient_valuations: tuple[int, ...] = ()
    recursive_contraction: bool = False
    theorem_1_1_complete: bool = False
    lean_proof_complete: bool = False


def _refuse(reason: str, **details) -> ModeWitness:
    return ModeWitness(True, reason, count=None, **details)


def _layout(request: ModeRequest, budgets: ModeBudgets):
    _integer(request.a, 1, "a")
    _integer(request.edges, 0, "edges")
    _integer(request.precision, 1, "precision")
    _integer(request.row_symbol, 0, "row_symbol")
    _integer(request.column_symbol, 0, "column_symbol")
    _integer(request.packet, 0, "packet")
    if not isinstance(request.t, int) or isinstance(request.t, bool) or request.t == 0:
        raise ValueError("t must be a nonzero integer")
    length = bounded_power_of_two(request.a, budgets.max_packet_slots, "packet slots")
    if not -length <= request.t < length:
        raise ValueError("t must lie in the Cauchy grid excluding zero")
    coordinate = valuation2(request.t)
    if coordinate >= request.a:
        raise ValueError("mode expansion requires v2(t) < a")
    check_limit(request.precision, budgets.max_precision_bits, "precision bits")
    columns = request.real_columns
    check_limit(len(columns), budgets.max_real_columns, "real columns")
    if list(columns) != sorted(set(columns)) or any(slot < 0 or slot >= length for slot in columns):
        raise ValueError("real columns must be strictly increasing slot indices")
    block = length >> coordinate
    tile_width = 1 << (coordinate + 1)
    crossing = request.row_symbol != request.column_symbol
    tiles = tile_width if crossing else 0
    return length, coordinate, block, tile_width, crossing, tiles


def _weights(tiles: int, precision: int) -> tuple[int, ...]:
    return tuple(weight for _tile in range(tiles) for weight in range(1, precision))


def _tuple_product(weights: tuple[int, ...], precision: int, limit: int) -> int | None:
    bound = 1
    for weight in weights:
        choices = (precision - 1) // weight + 1
        if bound > limit // choices:
            return None
        bound *= choices
    return bound


def _enumerate(weights: tuple[int, ...], precision: int) -> tuple[tuple[tuple[int, ...], ...], int]:
    found: list[tuple[int, ...]] = []
    chosen = [0] * len(weights)
    visits = 0

    def rec(index: int, used: int) -> None:
        nonlocal visits
        visits += 1
        if index == len(weights):
            found.append(tuple(chosen))
            return
        weight = weights[index]
        span = (precision - 1 - used) // weight
        for beta in range(span + 1):
            chosen[index] = beta
            rec(index + 1, used + weight * beta)

    rec(0, 0)
    return tuple(found), visits


def _modes(request: ModeRequest, coordinate: int, tile_width: int) -> tuple[WeightedMode, ...]:
    modes = []
    for origin in range(tile_width):
        bit = (origin >> coordinate) & 1
        source = request.column_symbol if bit == 0 else request.row_symbol
        target = request.row_symbol if bit == 0 else request.column_symbol
        occupied = tuple(slot for slot in request.real_columns if slot % tile_width == origin)
        for weight in range(1, request.precision):
            modes.append(WeightedMode(
                request.packet, coordinate, origin, weight, source, target, occupied,
            ))
    return tuple(modes)


def build_mode_table(
    request: ModeRequest, *, budgets: ModeBudgets | None = None,
) -> ModeWitness:
    """Allocate a mode table only after the weight, tuple and visit bounds fit.

    Exponent tuples are those with weighted degree strictly below ``precision``.
    The returned count is always ``None``: a tuple list is not a solution count.
    """
    limits = budgets or ModeBudgets()
    try:
        _length, coordinate, block, tile_width, crossing, tiles = _layout(request, limits)
        weights = _weights(tiles, request.precision)
        check_limit(len(weights), limits.max_modes, "weighted modes")
        if crossing and tiles * block > 6 * request.edges:
            return _refuse(
                "one weight class exceeds 6*M/D; mode table not allocated",
            )
        product = _tuple_product(weights, request.precision, limits.max_tuples)
        if product is None or product > limits.max_tuples:
            return _refuse("exponent-tuple bound exceeds budget; mode table not allocated")
        if product * len(weights) + 1 > limits.max_partial_choices:
            return _refuse("partial-choice bound exceeds budget; mode table not allocated")
    except ModularBudgetRefusal as error:
        return _refuse(f"{error}; mode table not allocated")
    modes = () if not crossing else _modes(request, coordinate, tile_width)
    tuples, visits = _enumerate(weights, request.precision)
    if visits > len(weights) * len(tuples) + 1 or len(tuples) > product:
        raise ArithmeticError("exponent enumeration exceeded its preflight bound")
    description = FrozenModeDescription(
        request.packet, coordinate, request.row_symbol, request.column_symbol,
        request.precision, request.precision, modes,
    )
    return ModeWitness(
        False,
        "finite weighted-mode table checked" if crossing else "local packet has an empty mode table",
        description=description,
        tuples=tuples,
        partial_choices=visits,
        modes_allocated=True,
        tuples_allocated=True,
    )


def _fraction_power(value: Fraction, exponent: int) -> Fraction:
    result = Fraction(1)
    for _ in range(exponent):
        result *= value
    return result


def _series_terms(a: int, t: int, row: int, column: int, zeta: int, precision: int):
    length = 1 << a
    node = Fraction(t + (1 << (a + 1)) * zeta)
    coordinate = valuation2(t)
    tile_width = 1 << (coordinate + 1)
    origin = column & (tile_width - 1)
    gap = Fraction(column // tile_width)
    base_gap = node - origin + row
    numerator = node
    derivative = Fraction(1)
    for point in range(-length, length):
        if point == 0:
            continue
        numerator *= node - point
        derivative *= -point
    leading = numerator / derivative / base_gap
    ratio = Fraction(tile_width) / base_gap
    total = Fraction(0)
    valuations = []
    for index in range(precision - 1):
        coefficient = leading * _fraction_power(ratio, index)
        if coefficient.denominator % 2 == 0:
            raise ArithmeticError("crossing coefficient is not two-adically integral")
        if coefficient.numerator == 0:
            valuation = precision
        else:
            valuation = valuation2(coefficient.numerator) - valuation2(coefficient.denominator)
        if valuation < index + 1:
            raise ArithmeticError("crossing coefficient is not divisible by 2^weight")
        valuations.append(valuation)
        total += coefficient * _fraction_power(gap, index)
    if total.denominator % 2 == 0:
        raise ArithmeticError("truncated crossing series is not two-adically integral")
    modulus = 1 << precision
    residue = total.numerator * pow(total.denominator, -1, modulus) % modulus
    return residue, tuple(valuations)


def recover_crossing_residue(
    request: ModeRequest, row: int, column: int, zeta: int,
    *, budgets: ModeBudgets | None = None,
) -> ModeWitness:
    """Recover one crossing Cauchy entry from the truncated series (D.12).

    The direct value is the Step 20 local polynomial at the same node. The
    count is their common residue. A budget refusal does not evaluate it.
    """
    limits = budgets or ModeBudgets()
    _integer(row, 0, "row")
    _integer(column, 0, "column")
    _integer(zeta, 0, "zeta")
    witness = build_mode_table(request, budgets=limits)
    if witness.refused or witness.description is None:
        return witness
    coordinate = witness.description.coordinate
    length = 1 << request.a
    if row >= length or column >= length:
        raise ValueError("row and column must be packet slot indices")
    if ((row >> coordinate) & 1) != ((column >> coordinate) & 1):
        raise ValueError("local slots are kept unchanged and are not expanded")
    if request.row_symbol == request.column_symbol:
        raise ValueError("equal packet symbols have no crossing expansion")
    ring = ModularRing(request.precision)
    direct = evaluate_polynomial(
        ring,
        cauchy_local_polynomial(ring, request.a, request.t, row, column),
        ring.scalar(zeta),
    )
    residue, valuations = _series_terms(
        request.a, request.t, row, column, zeta, request.precision,
    )
    if direct[0] != residue:
        raise ArithmeticError("geometric series disagrees with the Cauchy entry")
    return ModeWitness(
        False,
        "crossing Cauchy entry recovered from the weighted geometric series",
        count=residue,
        description=witness.description,
        tuples=witness.tuples,
        partial_choices=witness.partial_choices,
        modes_allocated=True,
        tuples_allocated=True,
        coefficient_valuations=valuations,
    )


def support_admissible(record: SupportRecord) -> bool:
    """True when the declared X degree fits the U and Z degrees."""
    _integer(record.deg_u, 0, "deg_u")
    _integer(record.deg_z, 0, "deg_z")
    _integer(record.deg_x, 0, "deg_x")
    return record.deg_x <= record.deg_u + record.deg_z


def retained_degree_ok(record: SupportRecord, mu: int, degree_cap: int) -> bool:
    """Admissible degrees inside the open rectangle meet the retained bound.

    Outside that rectangle the predicate is vacuously true. It does not store
    a polynomial.
    """
    _integer(mu, 1, "mu")
    _integer(degree_cap, 1, "degree_cap")
    if not support_admissible(record):
        return False
    if record.deg_u < mu and record.deg_z < degree_cap:
        return record.deg_x <= mu + degree_cap - 2
    return True


def admit_support(
    record: SupportRecord, *, budgets: ModeBudgets | None = None,
) -> ModeWitness:
    """Store one admissible degree record, or refuse before storing it."""
    limits = budgets or ModeBudgets()
    if not support_admissible(record):
        return _refuse("inadmissible weighted support; degree record not stored")
    try:
        check_limit(1, limits.max_support_records, "support records")
    except ModularBudgetRefusal as error:
        return _refuse(f"{error}; degree record not stored")
    return ModeWitness(
        False,
        "weighted-support degrees admitted",
        support_stored=True,
    )


def admit_exponent_tuple(
    weights: Sequence[int], exponents: Sequence[int], precision: int,
    *, budgets: ModeBudgets | None = None,
) -> ModeWitness:
    """Keep one tuple only when its weighted degree is strictly below precision."""
    limits = budgets or ModeBudgets()
    _integer(precision, 1, "precision")
    if len(weights) != len(exponents):
        raise ValueError("weights and exponents must have the same length")
    for weight in weights:
        _integer(weight, 1, "weight")
    for exponent in exponents:
        _integer(exponent, 0, "exponent")
    weighted = sum(weight * exponent for weight, exponent in zip(weights, exponents))
    if weighted >= precision:
        return _refuse("exponent tuple meets or exceeds the precision budget; tuple not stored")
    try:
        check_limit(1, limits.max_tuples, "exponent tuples")
    except ModularBudgetRefusal as error:
        return _refuse(f"{error}; tuple not stored")
    return ModeWitness(
        False,
        "exponent tuple inside the precision budget",
        tuples=(tuple(exponents),),
        tuples_allocated=True,
    )


def recover_network_with_modes(
    network: BooleanTensorNetwork, layout: PacketLayout, u: int, request: ModeRequest,
    *, mode_budgets: ModeBudgets | None = None, recovery_budgets=None,
) -> ModeWitness:
    """Admit modes first. Step 20 recovery runs only after the table exists.

    A mode-budget refusal does not call recovery and returns no count. The
    recovery value is the existing tapered count; modes do not replace it and
    do not evaluate a packet tensor.
    """
    if layout.a != request.a:
        raise ValueError("mode request and packet layout use different packet lengths")
    witness = build_mode_table(request, budgets=mode_budgets)
    if witness.refused:
        return witness
    recovered = recover_network_modular(network, layout, u, budgets=recovery_budgets)
    if not recovered.recovered or recovered.value is None:
        return _refuse(
            recovered.reason,
            recovery_called=True,
            port_pairs_scanned=recovered.port_pairs_scanned,
            recovery_tensor_evaluations=recovered.tensor_evaluations,
        )
    base = network.ring.base_coefficients(recovered.value)
    return ModeWitness(
        False,
        "modes admitted; Step 20 recovery returned the base value",
        count=base[0],
        base_value=base,
        description=witness.description,
        tuples=witness.tuples,
        partial_choices=witness.partial_choices,
        modes_allocated=True,
        tuples_allocated=True,
        recovery_called=True,
        port_pairs_scanned=recovered.port_pairs_scanned,
        recovery_tensor_evaluations=recovered.tensor_evaluations,
    )
