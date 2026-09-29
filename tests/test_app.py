"""Render every bench both ways, plus button flows; no live network dependency."""
from pathlib import Path
import pytest
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
BENCHES = ['home', 'immunity_arena', 'genome_surgery', 'repair_lab', 'guide_design', 'cas_explorer', 'therapeutics', 'database']


def bench(name, lite=False):
    app = AppTest.from_string(f'''
import streamlit as st
from beastlab.pages.{name} import render
from beastlab.field_notes import render as notes
st.session_state['lite_mode'] = {lite!r}
notes({{'home': 'Mission control', 'guide_design': 'Guide design'}}.get({name!r}, 'Therapeutics'))
render()
''', default_timeout=30)
    return app.run()


@pytest.mark.parametrize('name', BENCHES)
@pytest.mark.parametrize('lite', [False, True])
def test_all_benches(name, lite):
    at = bench(name, lite)
    assert not at.exception


def button(at, label):
    return next(b for b in at.button if b.label == label)


def test_root():
    at = AppTest.from_file(str(ROOT / 'app.py'), default_timeout=30).run()
    assert not at.exception
    at.toggle(key='lite_mode').set_value(True).run()
    assert at.session_state['lite_mode'] and not at.exception


def test_case_run():
    at = bench('therapeutics')
    button(at, 'Run case: Cas9 cut → NHEJ → HbF').click().run()
    assert not at.exception
    assert any('Functional cure predicted' in x.value for x in at.info)


def test_guide_ml_rl_export():
    at = bench('guide_design', True)
    for label in ['Predict Efficiency', 'Show training loss', '🧠 Auto-Optimize Guide', 'Scan bundled phiX174 offline', 'Prepare exports in RAM']:
        button(at, label).click().run()
        assert not at.exception
    assert at.session_state['priority_rl'][1][0]['policy'] == 'DQN'
    assert len(at.session_state['guide_tools_exports'][1]) == 3


def test_design_and_stale_input():
    at = bench('guide_design')
    button(at, '🔍 Design and rank guides').click().run()
    assert not at.exception and at.session_state['guides']
    next(s for s in at.slider if s.label == 'Off-target mismatch budget').set_value(2).run()
    assert not at.exception
    assert not at.session_state.get('guides', None)


def test_grch38_mocked(monkeypatch):
    guide = 'ACGT'*5
    monkeypatch.setattr('beastlab.offtarget.reference_or_fallback', lambda _: (
        guide+'TGG', {'source': 'Ensembl REST (test fixture)', 'start': 100, 'end': 122, 'chromosome': '11'}, None))
    at = bench('guide_design')
    button(at, 'Scan vs GRCh38').click().run()
    assert not at.exception
    assert at.session_state['priority_scan'][1]['hits'][0]['position'] == '11:100-119'


def test_retrain_button():
    at = bench('guide_design', True)
    button(at, 'Retrain MLP · seed=42').click().run()
    assert not at.exception
    assert at.session_state['retrained_efficiency'][1]['n_guides'] == 5000


def test_torch_lazy_import():
    import subprocess
    import sys
    result = subprocess.run([sys.executable, '-c',
        "import sys; from beastlab.pages import guide_design, genome_surgery, therapeutics; import beastlab.priority_ui; assert 'torch' not in sys.modules"], capture_output=True)
    assert result.returncode == 0, result.stderr


def test_clear_session_retains_lite():
    at = AppTest.from_file(str(ROOT / 'app.py'), default_timeout=30).run()
    at.toggle(key='lite_mode').set_value(True).run()
    at.session_state['user_private_result'] = 'private'
    button(at, 'Clear session inputs & results').click().run()
    assert not at.exception and at.session_state['lite_mode']
    assert 'user_private_result' not in at.session_state
