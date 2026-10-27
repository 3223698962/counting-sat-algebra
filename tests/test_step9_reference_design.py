"""Reference Lemma 4.5 expectations and Definition 3.1's rational kernel."""

from __future__ import annotations

import itertools
import random
from dataclasses import replace
from fractions import Fraction
from types import SimpleNamespace

import pytest

from counting_sat_algebra.cnf import CNF
from counting_sat_algebra.conditional_balance import choose_schedules
from counting_sat_algebra.linear_algebra import RationalRowSpan
from counting_sat_algebra.replacement_matrix import (
    build_auxiliary_slice, design_problems_at, step3_examples, vector_in_row_span,
)
from counting_sat_algebra.scheduler import (
    _finish_chosen_placement, _search_schedules, check_balance,
    schedule_placement, scopes_from_cnf, step7_examples,
)


def _potential(n, scopes, r, Q, variables, factors):
    """Direct reference Phi and Xi, independent of the probability formulas."""
    counts = [sum(v in scope for scope in scopes) for v in range(n)]
    M = sum(counts)
    norm = M * (max(counts) + max(map(len, scopes)))
    labels = {}
    for e, scope in enumerate(scopes):
        for v in scope:
            key = variables[v], factors[e]
            labels[key] = labels.get(key, 0) + 1
    phi = Fraction(sum(c * c for c in labels.values())) - Fraction(M * M, r ** (2 * Q))
    xi = sum(
        (sum(sigma[t] == g for sigma in variables) - Fraction(n, r)) ** 2
        for t in range(Q) for g in range(r)
    )
    return phi / norm + xi / (Q * n)


@pytest.mark.parametrize("r,Q", [(2, 1), (2, 2), (3, 1)])
def test_conditional_expectation_matches_all_completions(r, Q):
    n, scopes = 3, ((0, 1), (1, 2), ())
    chosen = choose_schedules(n, scopes, r, Q)
    alphabet = tuple(itertools.product(range(r), repeat=Q))
    slots = [('v', v) for v in range(n)] + [('f', e) for e, scope in enumerate(scopes) if scope]
    prefix = {}
    path_index = 0

    def expectation(fixed):
        pending = [slot for slot in slots if slot not in fixed]
        total = Fraction(0)
        for symbols in itertools.product(alphabet, repeat=len(pending)):
            filled = dict(fixed)
            filled.update(zip(pending, symbols))
            variables = tuple(filled['v', v] for v in range(n))
            factors = tuple(filled['f', e] if scope else alphabet[0] for e, scope in enumerate(scopes))
            total += _potential(n, scopes, r, Q, variables, factors)
        return total / len(alphabet) ** len(pending)

    assert expectation(prefix) == chosen.initial_potential
    for kind, index in [('v', v) for v in range(n)] + [('f', e) for e in range(len(scopes))]:
        symbol = chosen.variables[index] if kind == 'v' else chosen.factors[index]
        if kind == 'v' or scopes[index]:
            candidates = [(expectation({**prefix, (kind, index): candidate}), candidate) for candidate in alphabet]
            assert min(candidates)[1] == symbol
            prefix[kind, index] = symbol
        else:
            assert symbol == alphabet[0]
        path_index += 1
        assert expectation(prefix) == chosen.expectation_path[path_index]
    assert chosen.final_potential <= chosen.initial_potential <= 2


def test_cached_row_span_matches_independent_dense_elimination():
    rng = random.Random(19)
    for width in range(6):
        for height in range(6):
            rows = [[Fraction(rng.randrange(-2, 3), rng.randrange(1, 4)) for _ in range(width)] for _ in range(height)]
            span = RationalRowSpan(rows, width)
            for _ in range(12):
                target = [Fraction(rng.randrange(-2, 3), 2) for _ in range(width)]
                assert span.contains(target) == vector_in_row_span(rows, target)
    span = RationalRowSpan([[Fraction(1), Fraction(1)], [Fraction(2), Fraction(2)]], 2)
    assert span.rank == 1
    with pytest.raises(ValueError):
        RationalRowSpan([[Fraction(1)]], 2)
    with pytest.raises(ValueError):
        span.contains([Fraction(0)])


def test_linear_design_agrees_with_small_boolean_checks_and_rejects_damage():
    for example in step3_examples() + step7_examples(heavy=False):
        placed = example.placed
        for t in range(1, placed.auxiliary.Q + 1):
            assert design_problems_at(placed, t) == design_problems_at(placed, t, boolean_check=True) == []
    placed = step3_examples()[0].placed
    broken = replace(placed, measurement_rows=())
    assert any("outside the measurement span" in p for p in design_problems_at(broken, 1))
    assert design_problems_at(broken, 1) == design_problems_at(broken, 1, boolean_check=True)
    malformed = replace(placed, measurement_rows=((Fraction(1),),))
    assert any("length" in p for p in design_problems_at(malformed, 1))


def test_admitted_1024_units_finish_without_a_boolean_walk(monkeypatch):
    auxiliary = build_auxiliary_slice(2, 2, 1)

    def boom(*args, **kwargs):
        raise AssertionError("the Boolean design cube was enumerated")

    monkeypatch.setattr('counting_sat_algebra.replacement_matrix.itertools', SimpleNamespace(product=boom))
    cnf = CNF(1024, [[(v, bool(v % 2))] for v in range(1024)])
    outcome = schedule_placement(cnf, auxiliary=auxiliary)
    assert outcome.balance.admitted and outcome.balance.degree_gap_ok
    assert outcome.accepted and outcome.placed is not None
    assert outcome.balancer_built and outcome.design_checked and not outcome.caller_supplied
    assert outcome.design_method == 'linear'
    assert outcome.multiplicity == 64 and outcome.multiplicity_bound == 128
    assert outcome.max_group_size == outcome.outside_group_size == 512
    assert outcome.final_potential <= outcome.initial_potential <= 2
    assert outcome.design_problems == outcome.restoration_problems == ()


def test_degree_gap_and_allocation_gates_are_kept():
    cnf = CNF(64, [[(v, False)] for v in range(64)])
    outcome = schedule_placement(cnf, r=2, Q=1, R=1)
    assert outcome.design_checked and outcome.placed is not None
    assert not outcome.accepted and outcome.balance.degree_gap_ok is False
    assert 'degree gap' in outcome.reason
    with pytest.raises(ValueError, match='design_method'):
        schedule_placement(cnf, r=2, Q=1, R=1, design_method='unknown')


def test_placement_entry_budget_is_checked_before_slice_allocation(monkeypatch):
    cnf = CNF(8, [[(v, False)] for v in range(8)])
    balance = check_balance(cnf.n, scopes_from_cnf(cnf), r=2, Q=2, R=1)
    search = _search_schedules(cnf, balance)

    def boom(*args, **kwargs):
        raise AssertionError('a slice was allocated before the placement budget gate')

    monkeypatch.setattr('counting_sat_algebra.scheduler.BALANCER_PLACEMENT_ENTRY_BUDGET', 0)
    monkeypatch.setattr('counting_sat_algebra.scheduler.build_auxiliary_slice', boom)
    outcome = _finish_chosen_placement(cnf, balance, search, None)
    assert not outcome.built_slice and outcome.placed is None and not outcome.accepted
    assert outcome.balancer_attempted and not outcome.design_checked
    assert outcome.initial_potential == search.initial_potential
    assert 'coefficient entries' in outcome.reason and 'place_formula was not called' in outcome.reason
