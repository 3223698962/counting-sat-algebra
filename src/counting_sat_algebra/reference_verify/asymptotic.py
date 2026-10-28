"""Step 16 finite checks: Lemma 4.8's integer schedule and one budget-feasible member.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction
from dataclasses import asdict

from .protocol import _row, clean_refusal
from ..cnf import CNF
from ..common_count import count_placed_formula
from ..replacement_matrix import build_auxiliary_slice
from ..reference_parameters import (
    ChoiceBudgets,
    choose_asymptotic_parameters_for_cnf,
    lemma_4_8_integers,
)
from ..sparse_placement import schedule_sparse_placement
from ..support_count import count_placed_formula_by_support
from ..tensor_count import ContractionBudgets, count_placed_formula_by_contraction


def _alternating_units(n: int) -> CNF:
    return CNF(n, [[(v, bool(v % 2))] for v in range(n)])


def run_asymptotic_family_suite() -> list[dict]:
    """Choose Lemma 4.8's schedule and run a budget-feasible finite member.

    The closed form is checked exactly. On these n it fails before D is formed.
    The executed tuple is whatever ``check_reference_parameters`` accepts inside
    the sparse-placement budgets. Counting uses the existing backends and returns
    no partial count when a budget fails. This is not an asymptotic runtime proof.
    """
    print(
        '=== (R) Asymptotic parameter family: Lemma 4.8 and a finite admitted member ===',
        flush=True,
    )
    rows = []
    large_n = 16384
    nu, r, Q, R = lemma_4_8_integers(1 << 19, 1)
    guard_fails = r >= 2 and Q <= R
    rows.append(
        _row(
            'lemma48-integers-lambda-2**19',
            0,
            0,
            guard_fails,
            'CHECKED',
            log_n=1 << 19,
            log_log_n=nu,
            r=r,
            Q=Q,
            R=R,
            preliminary_guard_passed=False,
            dimension_formed=False,
            matrix_allocated=False,
            numerical_hypotheses_admitted=False,
            note='Q <= R, so tau and D are not formed; this lambda is not an executed instance',
        )
    )
    print(
        f'  [R] lambda=2**19: r={r}, Q={Q}, R={R}; preliminary guard {"fails before D" if guard_fails else "FAIL"}',
        flush=True,
    )

    cnf = _alternating_units(large_n)
    choice = choose_asymptotic_parameters_for_cnf(cnf)
    proposal = choice.proposal
    proposal_ok = (
        proposal.log_n == 14
        and proposal.log_log_n == 4
        and proposal.beta == 1
        and proposal.r == 0
        and proposal.Q == 0
        and proposal.R == 16
        and proposal.dimension_formed is False
        and proposal.dimension is None
        and proposal.tau is None
        and proposal.admitted is False
        and proposal.lemma_4_7_admitted is False
        and any('r=0 < 2' in reason for reason in proposal.reasons)
    )
    rows.append(
        _row(
            'lemma48-units-n16384',
            cnf.n,
            cnf.m,
            proposal_ok,
            'CHECKED',
            log_n=proposal.log_n,
            log_log_n=proposal.log_log_n,
            beta=proposal.beta,
            r=proposal.r,
            Q=proposal.Q,
            R=proposal.R,
            tau=None,
            dimension=None,
            dimension_formed=proposal.dimension_formed,
            matrix_allocated=False,
            numerical_hypotheses_admitted=proposal.admitted,
            reasons=list(proposal.reasons),
            placement_started=False,
        )
    )
    print(
        f'  [R] n=16384 Lemma 4.8: r={proposal.r}, Q={proposal.Q}, R={proposal.R}; D not formed; hypotheses rejected',
        flush=True,
    )

    smaller = _alternating_units(8192)
    smaller_choice = choose_asymptotic_parameters_for_cnf(smaller)
    smaller_ok = (
        smaller_choice.source == 'none'
        and smaller_choice.selected is None
        and smaller_choice.admitted is False
        and any('binomial tail 7/16' in item for item in smaller_choice.obstacles)
        and any('2(b+w)D^2 = 16384' in item for item in smaller_choice.obstacles)
    )
    rows.append(
        _row(
            'finite-member-units-n8192',
            smaller.n,
            smaller.m,
            smaller_ok,
            'CHECKED',
            source=smaller_choice.source,
            numerical_hypotheses_admitted=False,
            examined=smaller_choice.examined,
            admitted_candidates=smaller_choice.admitted_candidates,
            obstacles=list(smaller_choice.obstacles),
            reasons=list(smaller_choice.reasons),
            placement_started=False,
            matrix_allocated=False,
        )
    )
    print(
        f'  [R] n=8192: no Lemma 4.7 tuple inside the budgets; {"CHECKED failure" if smaller_ok else "FAIL"}',
        flush=True,
    )

    selected = choice.selected
    choice_ok = (
        proposal_ok
        and choice.source == 'finite-budget-member'
        and selected is not None
        and selected.placeable
        and choice.admitted
        and selected.reference.admitted
        and (selected.r, selected.Q, selected.R, selected.tau, selected.A0, selected.dimension)
        == (2, 3, 2, Fraction(5, 16), 0, 64)
        and selected.tail_probability == Fraction(5, 32)
        and selected.tail_limit == Fraction(5, 32)
        and selected.meets_r_at_least_32 is False
        and selected.meets_six_r_at_most_q is False
        and 2 * selected.tau * cnf.n == 10240
    )
    rows.append(
        _row(
            'finite-member-units-n16384',
            cnf.n,
            cnf.m,
            bool(choice_ok),
            'CHECKED',
            source=choice.source,
            numerical_hypotheses_admitted=bool(selected and selected.reference.admitted),
            r=None if selected is None else selected.r,
            Q=None if selected is None else selected.Q,
            R=None if selected is None else selected.R,
            tau=None if selected is None else str(selected.tau),
            hand_tau='1/2',
            measurement_allowance=None if selected is None else int(2 * selected.tau * cnf.n),
            hand_measurement_allowance=16384,
            interval=None if selected is None else list(selected.interval),
            dimension=None if selected is None else selected.dimension,
            tail=None if selected is None else str(selected.tail_probability),
            tail_limit=None if selected is None else str(selected.tail_limit),
            alphabet=None if selected is None else selected.alphabet,
            trials=None if selected is None else selected.trials,
            examined=choice.examined,
            admitted_candidates=choice.admitted_candidates,
            meets_r_at_least_32=None if selected is None else selected.meets_r_at_least_32,
            meets_six_r_at_most_q=None if selected is None else selected.meets_six_r_at_most_q,
            count_checked=False,
            note='Least tau meeting the exact tail. Not the Lemma 4.8 closed form, and not R>=32, 6R<=Q.',
        )
    )
    print(
        f'  [R] n=16384 finite member: r=2, Q=3, R=2, tau=5/16; {"CHECKED" if choice_ok else "FAIL"}',
        flush=True,
    )

    tight = choose_asymptotic_parameters_for_cnf(cnf, budgets=ChoiceBudgets(max_trials=1))
    tight_ok = (
        tight.admitted
        and tight.selected is not None
        and tight.placeable is False
        and tight.selected.r == 2
        and tight.selected.trials > 1
    )
    rows.append(
        _row(
            'finite-member-trials-budget',
            cnf.n,
            cnf.m,
            tight_ok,
            'REFUSED',
            source=tight.source,
            numerical_hypotheses_admitted=tight.admitted,
            placeable=tight.placeable,
            trials=None if tight.selected is None else tight.selected.trials,
            trial_budget=1,
            placement_started=False,
            matrix_allocated=False,
            reasons=list(tight.reasons),
        )
    )
    print(
        f'  [R] trial budget 1: {"REFUSED before schedule search" if tight_ok else "FAIL"}',
        flush=True,
    )

    try:
        build_auxiliary_slice(2, 192, 32)
        dense_refused = False
    except ValueError as exc:
        dense_refused = '2**384' in str(exc) and 'not built' in str(exc)
    rows.append(
        _row(
            'dense-full-m0-still-refused',
            0,
            0,
            dense_refused,
            'REFUSED',
            r=2,
            Q=192,
            R=32,
            matrix_allocated=False,
            note='The dense side-length cap is unchanged. Step 16 does not allocate full M0.',
        )
    )
    print(
        f'  [R] dense M0 at R=32, Q=192: {"REFUSED before allocation" if dense_refused else "FAIL"}',
        flush=True,
    )

    if not choice_ok or selected is None:
        rows.append(
            _row(
                'chosen-placement-n16384',
                cnf.n,
                cnf.m,
                False,
                'CHECKED',
                reason='chooser did not admit a placeable tuple',
            )
        )
        return rows

    outcome = schedule_sparse_placement(
        cnf, r=selected.r, Q=selected.Q, R=selected.R, tau=selected.tau, A0=selected.A0
    )
    storage = outcome.placed.statistics if outcome.placed is not None else None
    rank_sum = (
        sum(outcome.design.group_ranks[0])
        if outcome.design and outcome.design.group_ranks
        else None
    )
    placed_ok = (
        outcome.accepted
        and outcome.parameters.admitted
        and outcome.parameters.tau == Fraction(5, 16)
        and storage is not None
        and storage.polynomial_entries == 59392
        and storage.measurement_entries == 14336
        and storage.dense_polynomial_slots == 536870912
        and outcome.placed.copies == 256
        and outcome.design.max_group_size == 8192
        and outcome.design.group_ranks == ((2048, 2048),) * 3
        and rank_sum is not None
        and rank_sum <= 2 * outcome.parameters.tau * cnf.n
    )
    rows.append(
        _row(
            'chosen-placement-n16384',
            cnf.n,
            cnf.m,
            bool(placed_ok),
            'CHECKED',
            numerical_hypotheses_admitted=outcome.parameters.admitted,
            tau=str(outcome.parameters.tau),
            storage=asdict(storage) if storage else None,
            multiplicity=outcome.placed.copies if outcome.placed else None,
            group_ranks=[list(ranks) for ranks in outcome.design.group_ranks]
            if outcome.design
            else [],
            measurement_allowance=int(2 * outcome.parameters.tau * cnf.n),
            independent_measurements=rank_sum,
            count_checked=False,
            reason=outcome.reason,
        )
    )
    print(
        f'  [R] n=16384 sparse placement at the chosen tau: {"CHECKED" if placed_ok else "FAIL"}',
        flush=True,
    )
    if not placed_ok or outcome.placed is None:
        return rows

    counted = count_placed_formula(outcome.placed)
    common_refused = clean_refusal(counted.status, counted.count, counted.boolean_vectors)
    rows.append(
        _row(
            'chosen-n16384-common-budget',
            cnf.n,
            cnf.m,
            common_refused,
            'REFUSED',
            count=None,
            q=counted.plan.q,
            J=counted.plan.J,
            boolean_vector_bound_bits=counted.plan.boolean_vector_bound.bit_length(),
            boolean_vectors=counted.boolean_vectors,
            reason=counted.reason,
        )
    )
    print(
        f'  [R] common backend: {"REFUSED before enumeration" if common_refused else "FAIL"}',
        flush=True,
    )

    support = count_placed_formula_by_support(outcome.placed, vector_budget=49152)
    support_ok = (
        support.counted
        and support.count == 1
        and support.plan is not None
        and support.plan.components == 10240
        and support.plan.largest_component == 4
        and support.plan.max_component_degree == 4
        and support.plan.boolean_vector_bound == 49152
        and support.boolean_vectors == 96
        and support.local_evaluations == 13
        and support.memo_hits == 10227
    )
    rows.append(
        _row(
            'chosen-n16384-support',
            cnf.n,
            cnf.m,
            support_ok,
            'COUNTED',
            count=support.count,
            expected_count=1,
            count_validation='one assignment for independent alternating-polarity unit clauses',
            same_chosen_placement=True,
            memo_enabled=True,
            plan=asdict(support.plan) if support.plan else None,
            boolean_vectors=support.boolean_vectors,
            local_evaluations=support.local_evaluations,
            memo_hits=support.memo_hits,
            note='Disjoint replacement supports are not a general connected-instance bound.',
        )
    )
    print(
        f'  [R] support backend: #F={support.count}; vectors={support.boolean_vectors}; {"COUNTED" if support_ok else "FAIL"}',
        flush=True,
    )

    tensor = count_placed_formula_by_contraction(outcome.placed, budgets=ContractionBudgets())
    if tensor.counted:
        tensor_ok = (
            tensor.count == 1
            and tensor.plan is not None
            and tensor.plan.components == 10240
            and tensor.plan.largest_component == 4
            and tensor.plan.max_bucket_variables == 4
            and tensor.boolean_vectors == 176128
            and tensor.factor_evaluations == 98304
            and tensor.jet_products == 262144
            and tensor.peak_table_entries == 72
            and tensor.plan.planning_steps == 268288
            and tensor.plan.peak_table_coefficients == 360
        )
        rows.append(
            _row(
                'chosen-n16384-tensor',
                cnf.n,
                cnf.m,
                bool(tensor_ok),
                'COUNTED',
                count=tensor.count,
                expected_count=1,
                count_validation='one assignment for independent alternating-polarity unit clauses',
                same_chosen_placement=True,
                memo_enabled=False,
                plan=asdict(tensor.plan) if tensor.plan else None,
                boolean_vectors=tensor.boolean_vectors,
                factor_evaluations=tensor.factor_evaluations,
                jet_products=tensor.jet_products,
                peak_table_entries=tensor.peak_table_entries,
                note='Finite Taylor elimination on these supports. Not Appendix D.',
            )
        )
        print(
            f'  [R] tensor backend: #F={tensor.count}; vectors={tensor.boolean_vectors}; {"COUNTED" if tensor_ok else "FAIL"}',
            flush=True,
        )
    else:
        tensor_refused = clean_refusal(tensor.status, tensor.count, tensor.boolean_vectors)
        rows.append(
            _row(
                'chosen-n16384-tensor',
                cnf.n,
                cnf.m,
                tensor_refused,
                'REFUSED',
                count=None,
                boolean_vectors=tensor.boolean_vectors,
                reason=tensor.reason,
                plan=asdict(tensor.plan) if tensor.plan else None,
                note='Refused before numerical tables. Not a partial count.',
            )
        )
        print(
            f'  [R] tensor backend: {"REFUSED before numerical tables" if tensor_refused else "FAIL"}',
            flush=True,
        )
    return rows
