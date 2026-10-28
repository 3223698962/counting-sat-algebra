"""Write a finite reference report. Completion flags stay false."""

from __future__ import annotations

import json
from pathlib import Path

from .protocol import COMPLETION_FLAGS


def write_reference_report(
    path: str | Path,
    rows: list[dict],
    steps: list[int],
    *,
    failed_step: int | None = None,
) -> None:
    """Save the checks actually run, including failures, with their finite scope."""
    from .. import __version__

    report = dict(
        schema_version=1,
        version=__version__,
        goal='Reproduce the local reference exact-counting algorithm and verify its time and space bounds.',
        **COMPLETION_FLAGS,
        scope=(
            'Finite reference Sections 3 and 4 checks, sparse admitted placements, independent '
            'component counts, disjoint replacement supports, finite connected Taylor contraction, '
            'finite members of the Lemma 4.8 parameter schedule, finite Section 5.1–5.2 '
            'sparse-system checks, a finite Section 5.3–5.5 replacement sum in A_ell, and '
            'finite Section 6 equation selection with disjoint sparsification, and finite Appendix D '
            'coefficient rings, CNF networks, Cauchy deformation, tapered Hermite recovery, '
            'collision parameter schedules, packet layouts, prospective child topologies, and '
            'finite D.4.5 weighted modes. Proposition 6.4 does not claim an O_k(rho(n)) bound. '
            'Packet-tensor coefficients, Appendix D recursive contraction and Theorem 1.1 are '
            'not implemented; no general-input asymptotic claim.'
        ),
        wall_times_recorded=False,
        steps=steps,
        failed_step=failed_step,
        all_checks_passed=failed_step is None and bool(rows) and all(row['ok'] for row in rows),
        rows=rows,
    )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
