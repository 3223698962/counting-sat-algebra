"""Lemma 4.8's schedule and the budget-feasible finite member selected from it."""

from fractions import Fraction

import pytest

from counting_sat_algebra import __version__, verify
from counting_sat_algebra.cnf import CNF
from counting_sat_algebra.reference_parameters import (
    ChoiceBudgets, check_reference_parameters, choose_asymptotic_parameters_for_cnf,
    lemma_4_8_integers, lemma_4_8_proposal,
)
from counting_sat_algebra.reference_verify import run_asymptotic_family_suite
from counting_sat_algebra.replacement_matrix import build_auxiliary_slice
from counting_sat_algebra.scheduler import scopes_from_cnf


def _units(n, alternating=True):
    return CNF(n, [[(v, bool(v % 2) if alternating else False)] for v in range(n)])


def test_lemma_4_8_integers_fail_the_preliminary_guard_before_dimension():
    nu, r, Q, R = lemma_4_8_integers(1 << 19, beta=1)
    assert (nu, r, Q, R) == (19, 2, 215, 2900)
    assert r >= 2 and Q <= R
    # The guard stops before D = r**(2Q). The integer schedule itself has no dimension field.
    assert lemma_4_8_integers(14, 1) == (4, 0, 0, 16)


def test_closed_form_at_the_known_unit_point_does_not_form_d():
    cnf = _units(16384)
    proposal = lemma_4_8_proposal(cnf.n, scopes_from_cnf(cnf))
    assert (proposal.log_n, proposal.log_log_n, proposal.beta) == (14, 4, 1)
    assert (proposal.r, proposal.Q, proposal.R) == (0, 0, 16)
    assert proposal.dimension_formed is False and proposal.dimension is None
    assert proposal.tau is None and proposal.epsilon is None and proposal.interval is None
    assert proposal.admitted is False and proposal.lemma_4_7_admitted is False
    assert proposal.reference is None
    assert any('r=0 < 2' in reason for reason in proposal.reasons)


def test_chooser_matches_lemma_4_7_and_improves_the_hand_tau():
    cnf = _units(16384)
    scopes = scopes_from_cnf(cnf)
    first = choose_asymptotic_parameters_for_cnf(cnf)
    second = choose_asymptotic_parameters_for_cnf(cnf)
    assert first.source == second.source == 'finite-budget-member'
    assert first.examined == second.examined
    selected = first.selected
    assert selected is not None and selected.placeable
    assert (selected.r, selected.Q, selected.R, selected.tau, selected.A0) == (2, 3, 2, Fraction(5, 16), 0)
    assert selected.tail_probability == Fraction(5, 32)
    assert selected.reference == check_reference_parameters(
        cnf.n, scopes, r=2, Q=3, R=2, tau=Fraction(5, 16), A0=0,
    )
    hand = check_reference_parameters(cnf.n, scopes, r=2, Q=3, R=2, tau=Fraction(1, 2))
    assert hand.admitted and hand.tail_probability == Fraction(5, 32)
    assert selected.tau < Fraction(1, 2)
    assert 2 * selected.tau * cnf.n == 10240
    assert selected.meets_r_at_least_32 is False and selected.meets_six_r_at_most_q is False
    paper = check_reference_parameters(cnf.n, scopes, r=2, Q=3, R=2, tau=Fraction(1, 28))
    assert paper.admitted is False


def test_smaller_unit_formula_fails_a_numerical_hypothesis():
    cnf = _units(8192, alternating=False)
    choice = choose_asymptotic_parameters_for_cnf(cnf)
    assert choice.source == 'none' and choice.selected is None and choice.admitted is False
    assert choice.placeable is False
    assert any('binomial tail 7/16' in item for item in choice.obstacles)
    assert any('2(b+w)D^2 = 16384 > M = 8192' in item for item in choice.obstacles)
    assert 'no budget-bounded tuple' in choice.reasons[0]


def test_trial_budget_refuses_before_schedule_search(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('schedule search or a dense slice started')

    monkeypatch.setattr('counting_sat_algebra.sparse_placement.choose_schedules', boom)
    monkeypatch.setattr('counting_sat_algebra.sparse_placement.build_structured_auxiliary', boom)
    monkeypatch.setattr('counting_sat_algebra.replacement_matrix.build_auxiliary_slice', boom)
    cnf = _units(16384)
    choice = choose_asymptotic_parameters_for_cnf(cnf, budgets=ChoiceBudgets(max_trials=1))
    assert choice.admitted and choice.placeable is False
    assert choice.selected is not None and choice.selected.trials == 262144
    assert choice.reasons and 'before any schedule search' in choice.reasons[0]


def test_dimension_budget_reports_no_admitted_tuple_without_a_matrix(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('a matrix or schedule search started')

    monkeypatch.setattr('counting_sat_algebra.replacement_matrix.build_auxiliary_slice', boom)
    monkeypatch.setattr('counting_sat_algebra.sparse_placement.schedule_sparse_placement', boom)
    cnf = _units(16384)
    choice = choose_asymptotic_parameters_for_cnf(cnf, budgets=ChoiceBudgets(max_dimension=32))
    assert choice.selected is None and choice.admitted is False
    assert any('D=64 exceeds the choice dimension budget 32' in item for item in choice.obstacles)


def test_dense_side_length_cap_still_refuses_full_m0():
    with pytest.raises(ValueError, match=r'2\*\*384') as exc:
        build_auxiliary_slice(2, 192, 32)
    assert 'not built' in str(exc.value)


def test_step16_cli_is_wired(tmp_path, monkeypatch):
    row = dict(name='stub', n=1, m=1, brute=-1, ok=True, label='CHECKED')
    monkeypatch.setattr(verify, 'REFERENCE_SUITES', ((16, lambda: [row]),))
    output = tmp_path / 'step16.json'
    from counting_sat_algebra.__main__ import main
    assert main(['verify', '--step16', '--output', str(output)]) == 0
    import json
    report = json.loads(output.read_text(encoding='utf-8'))
    assert report['version'] == __version__ == '0.18.0'
    assert report['steps'] == [16] and report['all_checks_passed'] is True
    assert report['theorem_1_1_complete'] is report['lean_proof_complete'] is False
    assert 'Lemma 4.8' in report['scope']


def test_executed_family_checks_places_and_counts_or_refuses():
    rows = {row['name']: row for row in run_asymptotic_family_suite()}
    assert rows['lemma48-integers-lambda-2**19']['ok'] and rows['lemma48-integers-lambda-2**19']['label'] == 'CHECKED'
    assert rows['lemma48-integers-lambda-2**19']['dimension_formed'] is False
    assert rows['lemma48-units-n16384']['ok'] and rows['lemma48-units-n16384']['label'] == 'CHECKED'
    assert rows['lemma48-units-n16384']['numerical_hypotheses_admitted'] is False
    assert rows['finite-member-units-n8192']['ok'] and rows['finite-member-units-n8192']['label'] == 'CHECKED'
    assert rows['finite-member-units-n8192']['numerical_hypotheses_admitted'] is False
    chosen = rows['finite-member-units-n16384']
    assert chosen['ok'] and chosen['label'] == 'CHECKED'
    assert chosen['tau'] == '5/16' and chosen['measurement_allowance'] == 10240
    assert chosen['meets_r_at_least_32'] is False and chosen['meets_six_r_at_most_q'] is False
    assert rows['finite-member-trials-budget']['label'] == 'REFUSED' and rows['finite-member-trials-budget']['placement_started'] is False
    assert rows['dense-full-m0-still-refused']['label'] == 'REFUSED' and rows['dense-full-m0-still-refused']['matrix_allocated'] is False
    placed = rows['chosen-placement-n16384']
    assert placed['ok'] and placed['label'] == 'CHECKED'
    assert placed['storage']['polynomial_entries'] == 59392
    assert placed['storage']['measurement_entries'] == 14336
    assert placed['storage']['dense_polynomial_slots'] == 536870912
    assert placed['independent_measurements'] == 4096
    common = rows['chosen-n16384-common-budget']
    assert common['ok'] and common['label'] == 'REFUSED' and common['count'] is None and common['boolean_vectors'] == 0
    support = rows['chosen-n16384-support']
    assert support['ok'] and support['label'] == 'COUNTED' and support['count'] == 1
    assert support['boolean_vectors'] == 96 and support['local_evaluations'] == 13 and support['memo_hits'] == 10227
    tensor = rows['chosen-n16384-tensor']
    assert tensor['ok'] and tensor['label'] == 'COUNTED' and tensor['count'] == 1
    assert tensor['boolean_vectors'] == 176128 and tensor['jet_products'] == 262144
    assert tensor['peak_table_entries'] == 72 and tensor['factor_evaluations'] == 98304
    assert all(row['ok'] for row in rows.values())
