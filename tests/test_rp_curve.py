import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest
from reentrant_fea.rp_curve import engineering_row, ccx_rp_history


def test_relative_displacement_and_area():
    row = engineering_row(1., -.8, .4, -240., 239., 10., 20.)
    assert row['strain'] == pytest.approx(.12)
    assert row['nominal_stress_MPa'] == 12.
    assert row['force_balance_y_N'] == -1.
    with pytest.raises(ValueError):
        engineering_row(1., -.8, .4, -240., 239., 10., 0.)


def abaqus_exporter():
    path = Path(__file__).resolve().parents[1]/'export_abaqus_curve.py'
    spec = importlib.util.spec_from_file_location('abaqus_curve', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_abaqus_and_ccx_same_convention():
    exporter = abaqus_exporter()
    def region(u, f):
        return SimpleNamespace(historyOutputs={
            'U2': SimpleNamespace(data=[(1., u)]),
            'RF2': SimpleNamespace(data=[(1., f)])})
    rows = exporter.curve_rows(region(-.8, -240.), region(.4, 239.), 10., 20.)
    assert rows == [engineering_row(1., -.8, .4, -240., 239., 10., 20.)]
    bad = region(0., 0.)
    bad.historyOutputs['RF2'].data = [(0., 0.)]
    with pytest.raises(ValueError, match='sampling'):
        exporter.curve_rows(region(-.8, -240.), bad, 10., 20.)


def test_dat_filters_unconverged_and_rejects_missing_pairs(tmp_path):
    from reentrant_fea.paper_contact_solver import write_paper_contact_deck
    import numpy as np
    nodes = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                      [0,0,1],[1,0,1],[1,1,1],[0,1,1]], float)
    deck = tmp_path/'slab.inp'
    write_paper_contact_deck(deck, nodes, [[0,1,2,4,5,6],[0,2,3,4,6,7]],
        {'young':69000., 'poisson':.3, 'compression_strain':.12, 'reference_area':20.},
        [(225.,0.),(300.,.2)])
    deck.with_suffix('.sta').write_text('1 1 1 3 .5 .5 .5\n1 2 1U 50 .5 .5 .5\n')
    def block(t, region, kind, value):
        return f'{kind} (x,y,z) for set {region} and time {t}\n\n26 0 {value} 0\n'
    text = ''
    for time in (.5, 1.):
        for region in ('RP_TOP_PLATEN', 'RP_BOTTOM_PLATEN'):
            text += block(time, region, 'displacements', -.06 if 'TOP' in region else 0.)
            text += block(time, region, 'forces', -240. if 'TOP' in region else 240.)
    deck.with_suffix('.dat').write_text(text)
    rows, area = ccx_rp_history(deck)
    assert len(rows) == 1 and rows[0]['strain'] == .06 and area == 20.
    deck.with_suffix('.dat').write_text(text.replace(block(.5, 'RP_BOTTOM_PLATEN', 'forces', 240.), ''))
    with pytest.raises(ValueError, match='Missing paired'):
        ccx_rp_history(deck)
