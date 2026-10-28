"""Step 19 finite checks: Section 6 equation selection and disjoint sparsification.

Labels stay COUNTED, CHECKED, or REFUSED. A failed check is FAIL.
This module does not claim Theorem 1.1, a Lean proof, or recursive contraction.
"""

from __future__ import annotations

from fractions import Fraction

from .protocol import _row, completion_fields
from ..cnf import CNF, brute_count
from ..sparse_systems import SparseSystem, count_common_zeros
from ..equation_selection import (
    STEP19_SCOPE,
    Section6Budgets,
    choose_branching_M,
    construct_initial_matrix,
    count_via_disjoint_sparsification,
    count_violations_common_zeros,
    greedy_select_equations,
    kernel_support_check,
    largest_s,
    lemma_6_1_parameters,
    proposition_6_4_schedule,
)


def run_section6_suite() -> list[dict]:
    """Finite Section 6. Not Appendix D, not Proposition 6.4's runtime, not Theorem 1.1."""
    print('=== (U) Section 6 equation selection and disjoint sparsification ===', flush=True)
    rows = []
    parameters = lemma_6_1_parameters(3, 3, 1)
    built = construct_initial_matrix(3, 3, 1)
    matrix = built.matrix
    kernel = (
        None
        if matrix is None
        else kernel_support_check(
            matrix,
            numerator=parameters.m,
            denominator=parameters.P,
            max_boolean_vectors=64,
        )
    )
    kernel_ok = (
        built.constructed
        and matrix is not None
        and kernel is not None
        and kernel.ok
        and parameters.mu == 2
        and parameters.P == 80_000
        and parameters.D0 == 2_560_000
        and parameters.q_G == 1
        and parameters.h == 0
        and parameters.local_search_bits == 0
        and matrix.nrows == 2
        and matrix.nrows <= 3 * 3
        and matrix.max_row_weight == 3
        and matrix.max_row_weight <= parameters.D0
        and kernel.min_nonzero_weight == 2
        and kernel.nonzero_kernel_vectors == 3
        and kernel.vectors_scanned == 8
        and kernel.min_nonzero_weight * parameters.P >= parameters.m
        and built.rows_allocated
        and built.search_started is False
    )
    rows.append(
        _row(
            'lemma-6.1-tiny-kernel',
            3,
            3,
            kernel_ok,
            'CHECKED',
            mu=parameters.mu,
            P=parameters.P,
            D0=parameters.D0,
            rows=matrix.nrows if matrix else 0,
            max_row_weight=matrix.max_row_weight if matrix else None,
            min_nonzero_weight=None if kernel is None else kernel.min_nonzero_weight,
            nonzero_kernel_vectors=None if kernel is None else kernel.nonzero_kernel_vectors,
            vectors_scanned=0 if kernel is None else kernel.vectors_scanned,
            row_bound=9,
            count_checked=False,
            appendix_d_implemented=False,
            rho_bound_claimed=False,
            **completion_fields(),
            note='Even-parity A(0) on m=3. Every nonzero kernel vector has weight 2 >= 3/80000.',
        )
    )
    print(
        f'  [U] Lemma 6.1 m=3: min weight {None if kernel is None else kernel.min_nonzero_weight}; '
        f'{"CHECKED" if kernel_ok else "FAIL"}',
        flush=True,
    )

    search = construct_initial_matrix(2048, 2048, 1)
    search_bits = 0 if search.parameters is None else search.parameters.local_search_bits
    search_ok = (
        search.refused
        and search.matrix is None
        and search.rows_allocated is False
        and search.search_started is False
        and search_bits == 524_288
        and 'no local matrix allocated' in search.reason
    )
    rows.append(
        _row(
            'lemma-6.1-search-budget',
            2048,
            2048,
            search_ok,
            'REFUSED',
            count=None,
            local_search_bits=search_bits,
            rows_allocated=False,
            search_started=False,
            vectors_scanned=0,
            appendix_d_implemented=False,
            rho_bound_claimed=False,
            **completion_fields(),
            reason=search.reason,
        )
    )
    print(
        f'  [U] Lemma 6.1 search: {"REFUSED before allocation" if search_ok else "FAIL"}',
        flush=True,
    )

    tiny = CNF(2, [[(0, False)], [(1, False)]])
    greedy = greedy_select_equations(tiny, list_length=1, B=1)
    final = (
        None
        if greedy.equations is None
        else count_common_zeros(
            SparseSystem(tiny.n, 2, 1, greedy.equations),
        )
    )
    brute = brute_count(tiny)
    greedy_ok = (
        greedy.counted
        and greedy.count == brute == 1
        and final is not None
        and final.counted
        and final.count == 1
        and greedy.rounds_completed == 6
        and greedy.counts_performed == 24
        and greedy.candidates_enumerated == 24
        and greedy.boolean_vectors_scanned == 96
        and greedy.matrix_allocated
        and final.boolean_vectors == 4
    )
    rows.append(
        _row(
            'lemma-6.2-greedy-tiny',
            tiny.n,
            tiny.m,
            greedy_ok,
            'COUNTED',
            brute=brute,
            count=greedy.count,
            rounds=greedy.rounds_completed,
            list_length=1,
            counts_performed=greedy.counts_performed,
            boolean_vectors_scanned=greedy.boolean_vectors_scanned,
            common_zero_count=None if final is None else final.count,
            common_zero_vectors=0 if final is None else final.boolean_vectors,
            appendix_d_implemented=False,
            rho_bound_claimed=False,
            **completion_fields(),
            note='List length 1 fits. Final N(A) is 1, matching brute #F and count_common_zeros.',
        )
    )
    print(
        f'  [U] greedy n=2: N(A)={greedy.count}; counts={greedy.counts_performed}; '
        f'{"COUNTED" if greedy_ok else "FAIL"}',
        flush=True,
    )

    candidate = greedy_select_equations(
        tiny,
        list_length=1,
        B=1,
        budgets=Section6Budgets(max_candidate_lists=1),
    )
    candidate_ok = (
        candidate.refused
        and candidate.count is None
        and candidate.candidates_enumerated == 0
        and candidate.counts_performed == 0
        and candidate.boolean_vectors_scanned == 0
        and candidate.matrix_allocated is False
        and candidate.equations is None
        and 'no candidates enumerated' in candidate.reason
    )
    rows.append(
        _row(
            'lemma-6.2-candidate-budget',
            tiny.n,
            tiny.m,
            candidate_ok,
            'REFUSED',
            count=None,
            candidate_budget=1,
            candidates_enumerated=0,
            boolean_vectors_scanned=0,
            matrix_allocated=False,
            appendix_d_implemented=False,
            **completion_fields(),
            reason=candidate.reason,
        )
    )
    print(
        f'  [U] candidate budget: {"REFUSED before enumeration" if candidate_ok else "FAIL"}',
        flush=True,
    )

    cube = greedy_select_equations(
        tiny,
        list_length=1,
        B=1,
        budgets=Section6Budgets(max_boolean_vectors=1),
    )
    cube_ok = (
        cube.refused
        and cube.count is None
        and cube.candidates_enumerated == 0
        and cube.boolean_vectors_scanned == 0
        and cube.matrix_allocated is False
        and 'no assignments scanned' in cube.reason
    )
    rows.append(
        _row(
            'lemma-6.2-count-budget',
            tiny.n,
            tiny.m,
            cube_ok,
            'REFUSED',
            count=None,
            boolean_vector_budget=1,
            candidates_enumerated=0,
            boolean_vectors_scanned=0,
            matrix_allocated=False,
            appendix_d_implemented=False,
            **completion_fields(),
            reason=cube.reason,
        )
    )
    print(f'  [U] count budget: {"REFUSED before the cube" if cube_ok else "FAIL"}', flush=True)

    paper_P = lemma_6_1_parameters(tiny.n, tiny.m, 1).P
    paper = greedy_select_equations(tiny, list_length=paper_P, B=1)
    paper_ok = (
        paper_P == 80_000
        and paper.refused
        and paper.count is None
        and paper.candidates_enumerated == 0
        and paper.boolean_vectors_scanned == 0
        and paper.matrix_allocated is False
        and paper.equations is None
    )
    rows.append(
        _row(
            'lemma-6.2-paper-scale-P',
            tiny.n,
            tiny.m,
            paper_ok,
            'REFUSED',
            count=None,
            P=paper_P,
            candidates_enumerated=0,
            boolean_vectors_scanned=0,
            matrix_allocated=False,
            appendix_d_implemented=False,
            rho_bound_claimed=False,
            **completion_fields(),
            reason=paper.reason,
        )
    )
    print(f'  [U] paper-scale P={paper_P}: {"REFUSED" if paper_ok else "FAIL"}', flush=True)

    sparse = CNF(
        3,
        [
            [(0, False), (1, False)],
            [(0, False), (2, False)],
            [(1, False), (2, False)],
        ],
    )
    partitioned = count_via_disjoint_sparsification(
        sparse,
        k=2,
        M=2,
        leaf_counter=count_violations_common_zeros,
        verify_disjoint=True,
    )
    sparse_brute = brute_count(sparse)
    sparse_ok = (
        partitioned.counted
        and partitioned.count == sparse_brute == 4
        and partitioned.leaf_counts == (3, 1)
        and partitioned.leaves == 2
        and partitioned.nodes == 3
        and partitioned.depth == 1
        and partitioned.same_declared_n
        and partitioned.max_occurrence == 1
        and partitioned.occurrence_bound == 6
        and partitioned.max_occurrence <= partitioned.occurrence_bound
        and partitioned.tree is not None
        and all(leaf.n == sparse.n for leaf in partitioned.tree.leaves)
    )
    rows.append(
        _row(
            'disjoint-sparsification-n3',
            sparse.n,
            sparse.m,
            sparse_ok,
            'COUNTED',
            brute=sparse_brute,
            count=partitioned.count,
            leaf_counts=list(partitioned.leaf_counts),
            leaves=partitioned.leaves,
            nodes=partitioned.nodes,
            depth=partitioned.depth,
            occurrence_bound=partitioned.occurrence_bound,
            max_occurrence=partitioned.max_occurrence,
            same_declared_n=partitioned.same_declared_n,
            boolean_vectors_scanned=partitioned.boolean_vectors_scanned,
            step5_copying=False,
            appendix_d_implemented=False,
            rho_bound_claimed=False,
            **completion_fields(),
            note='Leaf common-zero counts sum to #F. Declared n stays 3. Step 5 copying is not this tree.',
        )
    )
    print(
        f'  [U] disjoint n=3: #F={partitioned.count}; leaves={list(partitioned.leaf_counts)}; '
        f'{"COUNTED" if sparse_ok else "FAIL"}',
        flush=True,
    )

    tree_budget = count_via_disjoint_sparsification(
        sparse,
        k=2,
        M=2,
        budgets=Section6Budgets(max_tree_nodes=1),
    )
    tree_ok = (
        tree_budget.refused
        and tree_budget.count is None
        and tree_budget.leaves == 0
        and tree_budget.boolean_vectors_scanned == 0
        and tree_budget.tree is None
        and tree_budget.leaf_counts == ()
        and 'no leaves returned' in tree_budget.reason
    )
    rows.append(
        _row(
            'disjoint-tree-budget',
            sparse.n,
            sparse.m,
            tree_ok,
            'REFUSED',
            count=None,
            node_budget=1,
            leaves=0,
            boolean_vectors_scanned=0,
            appendix_d_implemented=False,
            **completion_fields(),
            reason=tree_budget.reason,
        )
    )
    print(f'  [U] tree budget: {"REFUSED before the leaves" if tree_ok else "FAIL"}', flush=True)

    schedule = proposition_6_4_schedule(20)
    branching = choose_branching_M(2, Fraction(1))
    witness = largest_s(521_284, 19)
    schedule_ok = (
        schedule.lambda_n == 5
        and schedule.nu_n == 3
        and schedule.s == 0
        and schedule.exhaustive_search
        and schedule.sigma is None
        and schedule.rho_bound_claimed is False
        and schedule.theorem_1_1_complete is False
        and schedule.lean_proof_complete is False
        and witness == 2
        and branching.chosen
        and branching.M == 64
        and branching.C == 4
        and branching.h == 10
        and branching.gamma == 4
    )
    rows.append(
        _row(
            'proposition-6.4-schedule',
            20,
            0,
            schedule_ok,
            'CHECKED',
            lambda_n=schedule.lambda_n,
            nu_n=schedule.nu_n,
            s=schedule.s,
            sigma=schedule.sigma,
            exhaustive_search=True,
            witness_lambda=521_284,
            witness_nu=19,
            witness_s=witness,
            branching_M_at_sigma_1=branching.M,
            count_checked=False,
            rho_bound_claimed=False,
            appendix_d_implemented=False,
            **completion_fields(),
            note=(
                'n=20 has s=0, so the schedule is the exhaustive-search case. '
                's=2 is only an integer witness at lambda=521284. No rho bound is claimed.'
            ),
        )
    )
    print(
        f'  [U] schedule n=20: lambda={schedule.lambda_n}, nu={schedule.nu_n}, s={schedule.s}; '
        f'{"CHECKED" if schedule_ok else "FAIL"}',
        flush=True,
    )
    if STEP19_SCOPE and 'Theorem 1.1' not in STEP19_SCOPE:
        rows[-1]['ok'] = False
        rows[-1]['label'] = 'FAIL'
    return rows
