"""DIMACS parsing and the actual count command, including refusal output."""

import json
from pathlib import Path

import pytest

from counting_sat_algebra.__main__ import main
from counting_sat_algebra.cli import count_document
from counting_sat_algebra.cnf import CNF
from counting_sat_algebra.component_count import ComponentCountResult
from counting_sat_algebra.dimacs import format_dimacs, parse_dimacs


def test_multiline_clauses_comments_multiple_clauses_and_declared_free_bits():
    cnf = parse_dimacs('\ufeffc example\np cnf 8 3\n1 -2\nc between literals\n3 0 4 0\n0\n%\n0\n')
    assert cnf.n == 8 and cnf.clauses == [[(0, False), (1, True), (2, False)], [(3, False)], []]
    assert parse_dimacs(format_dimacs(cnf)) == cnf
    assert parse_dimacs(format_dimacs(CNF(0, []))) == CNF(0, [])


@pytest.mark.parametrize('text', ['', '1 0', 'p sat 1 1\n1 0', 'p cnf -1 0',
                                  'p cnf 1 0\np cnf 1 0', 'p cnf 1 1\n2 0',
                                  'p cnf 1 1\n1', 'p cnf 1 2\n1 0',
                                  'p cnf 1 1\na 0', 'p cnf 1 0\n0'])
def test_invalid_dimacs_is_rejected(text):
    with pytest.raises(ValueError):
        parse_dimacs(text)


def test_count_command_reports_the_exact_count_and_saves_json(tmp_path, capsys):
    source = tmp_path / 'example.cnf'
    source.write_text('p cnf 8 2\n1 2 -3 0\n4 -5 6 0\n', encoding='utf-8')
    output = tmp_path / 'report.json'
    assert main(['count', str(source), '--json', '--output', str(output)]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed == json.loads(output.read_text(encoding='utf-8'))
    assert printed['status'] == 'counted' and printed['count'] == '196'
    assert printed['free_variables'] == 2 and printed['satisfiable'] is True


def test_refused_command_returns_no_count(tmp_path, capsys):
    source = tmp_path / 'example.cnf'
    source.write_text('p cnf 2 1\n1 2 0\n', encoding='utf-8')
    assert main(['count', str(source), '--vector-budget', '0', '--json']) == 2
    printed = json.loads(capsys.readouterr().out)
    assert printed['status'] == 'refused' and printed['count'] is None
    assert printed['satisfiable'] is None and printed['boolean_vectors'] == 0


def test_huge_counts_have_exact_hex_output_without_a_decimal_digit_limit():
    count = 7 ** 10923
    document = count_document(ComponentCountResult('counted', count, 'test'))
    assert document['count_encoding'] == 'hex'
    assert int(document['count'], 16) == count
    assert json.loads(json.dumps(document)) == document
