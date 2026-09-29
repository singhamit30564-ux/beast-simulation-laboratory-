import csv
import hashlib
import io
import json
import pytest
from beastlab.io import parse_dna, export_bundle, canonical
from beastlab.offtarget import scan, reference_or_fallback
from beastlab.sequtils import revcomp
from beastlab.case_study import DEMO, DEMO_GUIDE, experiment
from beastlab.animation import figure


def test_input():
    assert parse_dna('>abc\nac gt\nTG') == 'ACGTTG'
    assert parse_dna('acgu'*5, guide=True) == 'ACGT'*5
    for text in ['', 'ACGN', 'ACG123', '>a\nACGT\n>b\nACGT', 'A'*10001]:
        with pytest.raises(ValueError):
            parse_dna(text)


def test_exports():
    files = export_bundle({'sequence': 'ACGT'}, {'score': 4}, 'test')
    obj = json.loads(files['json'])
    assert obj['sha256'] == hashlib.sha256(canonical(obj['payload'])).hexdigest()
    row = next(csv.DictReader(io.StringIO(files['csv'].decode())))
    assert json.loads(row['payload_json']) == obj['payload']
    assert row['sha256'] == obj['sha256']
    assert files['pdf'].startswith(b'%PDF')


def test_scan_strands_pam_seed_coordinates():
    guide = 'ACGTACGTACGTACGTACGA'
    target = guide+'TGG'
    plus = scan(target, guide, offset=100, chromosome='2')
    assert plus[0]['mismatches'] == 0 and plus[0]['position'] == '2:100-119'
    minus = scan(revcomp(target), guide, offset=100)
    assert minus[0]['strand'] == '-' and minus[0]['position'] == 'local:103-122'
    assert not scan(guide+'ATGG', guide)  # NGG cannot be shifted away from spacer
    mutated = guide[:-1]+'C'
    assert scan(mutated+'TGG', guide)[0]['seed mismatches (last 12)'] == 1
    assert not scan('TGCA'+guide[4:]+'TGG', guide)


def test_fallback(monkeypatch):
    def fail(_):
        raise OSError('offline')
    monkeypatch.setattr('beastlab.offtarget.fetch_reference', fail)
    seq, meta, error = reference_or_fallback('HBB')
    assert len(seq) == 5386 and 'NOT human' in meta['source'] and error == 'OSError'


def test_case():
    result = experiment(DEMO, DEMO_GUIDE, 6, 75, DEMO.find('AGATAA'))
    assert result['motif_disrupted_assumption'] and result['verdict'].startswith('✅')
    assert experiment(DEMO, DEMO_GUIDE, 0, 75, DEMO.find('AGATAA'))['verdict'].startswith('⚠')
    assert experiment(DEMO, 'A'*20, 6, 75, DEMO.find('AGATAA'))['binding_site'] is None


def test_animation():
    fig = figure()
    assert len(fig.frames) == 6
    assert all(trace.type == 'scatter' for f in fig.frames for trace in f.data)
    assert len(fig.layout.updatemenus[0].buttons) == 3
    assert not figure(3, False).frames


def test_legacy_scanner_adjacent_and_overlapping_pams():
    from beastlab.crispr import find_offtargets, scan_pams
    from beastlab.data.cas_enzymes import EFFECTOR_BY_ID
    enzyme = EFFECTOR_BY_ID['SpCas9']
    guide = 'ACGTACGTACGTACGTACGA'
    near = guide[:-1]+'C'
    plus = find_offtargets(near+'TGG', guide, enzyme, 3)
    minus = find_offtargets(revcomp(near+'TGG'), guide, enzyme, 3)
    assert plus[0]['seed_mismatches'] == minus[0]['seed_mismatches'] == 1
    assert minus[0]['strand'] == '-'
    assert not find_offtargets(near+'ATGG', guide, enzyme, 3)
    assert len([s for s in scan_pams('A'*20+'TGGG', enzyme) if s['strand'] == '+']) == 2


def test_notes_stable_and_complete():
    from beastlab.field_notes import NOTES
    assert len(NOTES) == 8
    assert all(3 <= len(notes) <= 4 for notes in NOTES.values())
