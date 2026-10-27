#!/usr/bin/env python3
"""Fail CI when a verify run looks complete but is not.

This script does not add theorem flags. It rejects a log or reference report
that is missing a finite success, or a report that claims Theorem 1.1 or the
Lean proof is complete.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ALLOWED_LABELS = frozenset({"COUNTED", "CHECKED", "REFUSED"})
REFERENCE_STEPS = list(range(9, 22))
QUICK_SUMMARY = re.compile(r"Reference \(K–W\): (\d+)/(\d+) checks passed")
THEOREM_LINE = "Theorem 1.1 is not implemented. This program does not refute ETH."
SCOPE_LINE = "No Theorem 1.1, no 2^{o(n)} claim."


def _fail(errors: list[str]) -> int:
    print("verify report check failed:", file=sys.stderr)
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    return 1


def check_quick(path: Path) -> int:
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []
    if "ALL TESTS PASSED" not in text:
        errors.append("quick log is missing ALL TESTS PASSED")
    if THEOREM_LINE not in text:
        errors.append("quick log is missing the Theorem 1.1 disclaimer")
    if SCOPE_LINE not in text:
        errors.append("quick log is missing the finite-scope line")
    match = QUICK_SUMMARY.search(text)
    passed = total = 0
    if match is None:
        errors.append("quick log is missing the Steps 9–21 summary")
    else:
        passed, total = int(match.group(1)), int(match.group(2))
        if total < 1 or passed != total:
            errors.append(f"quick reference summary is incomplete: {match.group(0)}")
    for label in sorted(ALLOWED_LABELS):
        if label not in text:
            errors.append(f"quick log never recorded {label}")
    if errors:
        return _fail(errors)
    print(f"quick log ok: {passed}/{total} reference checks; finite scope only")
    return 0


def check_reference(path: Path) -> int:
    report = json.loads(path.read_text(encoding="utf-8"))
    errors: list[str] = []
    rows = report.get("rows")
    if report.get("all_checks_passed") is not True:
        errors.append("all_checks_passed is not true")
    if report.get("failed_step") is not None:
        errors.append(f"failed_step is {report.get('failed_step')!r}")
    if report.get("theorem_1_1_complete") is not False:
        errors.append("theorem_1_1_complete must stay false")
    if report.get("lean_proof_complete") is not False:
        errors.append("lean_proof_complete must stay false")
    if report.get("steps") != REFERENCE_STEPS:
        errors.append(f"steps are {report.get('steps')!r}, expected {REFERENCE_STEPS}")
    if not isinstance(rows, list) or not rows:
        errors.append("report has no rows")
        rows = []
    elif not all(isinstance(row, dict) and row.get("ok") is True for row in rows):
        errors.append("a reference row is missing or not ok")
    labels = {row.get("label") for row in rows if isinstance(row, dict)}
    missing = ALLOWED_LABELS - labels
    extra = labels - ALLOWED_LABELS
    if missing:
        errors.append(f"report is missing labels: {sorted(missing)}")
    if extra:
        errors.append(f"report has unexpected labels: {sorted(extra)}")
    if errors:
        return _fail(errors)
    print(
        f"reference report ok: {len(rows)} finite checks; "
        "theorem_1_1_complete=false; lean_proof_complete=false"
    )
    return 0


def main(argv: list[str]) -> int:
    if len(argv) != 3 or argv[1] not in {"quick", "reference"}:
        print(
            "usage: check_verify_report.py quick LOG | reference REPORT.json",
            file=sys.stderr,
        )
        return 2
    path = Path(argv[2])
    if not path.is_file():
        print(f"missing verify output: {path}", file=sys.stderr)
        return 1
    if argv[1] == "quick":
        return check_quick(path)
    return check_reference(path)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
