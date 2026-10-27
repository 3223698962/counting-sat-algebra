"""Reference reports preserve actual checks, failures and their finite scope."""

import json

import pytest

from counting_sat_algebra import __version__, verify
from counting_sat_algebra.__main__ import main


def test_reference_output_records_executed_counts_and_scope(tmp_path, capsys):
    output = tmp_path / 'reports' / 'reference.json'
    assert main(['verify', '--step10', '--output', str(output)]) == 0
    report = json.loads(output.read_text(encoding='utf-8'))
    assert report['version'] == __version__
    assert report['steps'] == [10]
    assert report['all_checks_passed'] is True and report['failed_step'] is None
    assert report['theorem_1_1_complete'] is report['lean_proof_complete'] is False
    assert report['wall_times_recorded'] is False
    assert [row['count'] for row in report['rows']] == [1, 1, 7, 7]
    assert all(row['ok'] and row['label'] == 'COUNTED' for row in report['rows'])
    assert '4/4 checks passed' in capsys.readouterr().out


@pytest.mark.parametrize('failed_rows', [[], [dict(name='failed', ok=False, label='FAIL')]])
def test_failure_report_cannot_inherit_success_from_an_earlier_suite(tmp_path, monkeypatch, capsys, failed_rows):
    earlier = dict(name='earlier', n=0, m=0, brute=-1, ok=True, label='CHECKED')
    monkeypatch.setattr(verify, 'REFERENCE_SUITES', (
        (9, lambda: [earlier]), (14, lambda: failed_rows),
    ))
    output = tmp_path / 'failed.json'
    assert main(['verify', '--step9', '--step14', '--output', str(output)]) == 1
    report = json.loads(output.read_text(encoding='utf-8'))
    assert report['rows'] == [earlier, *failed_rows]
    assert report['steps'] == [9, 14] and report['failed_step'] == 14
    assert report['all_checks_passed'] is False
    assert 'FAIL: Step 14' in capsys.readouterr().err


@pytest.mark.parametrize('args', [[], ['--step8']])
def test_reference_output_requires_a_reference_suite(tmp_path, capsys, args):
    output = tmp_path / 'unused.json'
    with pytest.raises(SystemExit) as exc:
        main(['verify', *args, '--output', str(output)])
    assert exc.value.code == 2 and not output.exists()
    assert '--output requires' in capsys.readouterr().err
