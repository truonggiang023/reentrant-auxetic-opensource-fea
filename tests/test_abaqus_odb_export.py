"""Check ODB bookkeeping without claiming an Abaqus solver execution."""
import copy
import csv
import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS

import pytest


SOURCE = Path(__file__).resolve().parents[1] / 'abaqus' / 'extract_odb.py'
SPEC = importlib.util.spec_from_file_location('abaqus_odb_export', SOURCE)
export = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(export)


def output(values, times=(0., .005, .01)):
    return NS(data=list(zip(times, values)))


def mock_odb():
    # Labels repeat in separate parts: selecting just numeric labels is wrong.
    top = NS(nodes=((NS(label=1), NS(label=2)),), instanceNames=('PLATEN',))
    bottom = NS(nodes=((NS(label=1),),), instanceNames=('BASE',))
    histories = {
        'Node PLATEN.1': NS(historyOutputs={'U2': output([0., -1.2, -2.4]),
                                          'RF2': output([0., -19., -38.])}),
        'Node PLATEN.2': NS(historyOutputs={'U2': output([0., -1.2, -2.4]),
                                          'RF2': output([0., -38., -76.])}),
        'Node BASE.1': NS(historyOutputs={'RF2': output([0., 60., 120.])}),
        'Node CORE.1': NS(historyOutputs={'U2': output([0., -30., -60.]),
                                        'RF2': output([0., -999., -999.])}),
        'Assembly ASSEMBLY': NS(historyOutputs={name: output([0., 1., 2.])
                                               for name in export.ENERGIES}),
    }
    assembly = NS(nodeSets={'TOP_LOAD': top, 'BOTTOM_LOAD': bottom},
                  instances={'PLATEN': NS(), 'BASE': NS(), 'CORE': NS()})
    step = NS(historyRegions=histories, timePeriod=.01,
              frames=[NS(frameValue=0.), NS(frameValue=.01)])
    return NS(rootAssembly=assembly, steps={'COMPRESS': step})


def test_signed_force_sum_excludes_other_instances_and_keeps_inertia_balance(tmp_path):
    rows, assessment = export.extract(mock_odb())
    assert rows[-1]['stress_mpa'] == pytest.approx(114. / 1900.)
    assert rows[-1]['strain'] == pytest.approx(.05)
    assert rows[-1]['compression_mm'] == pytest.approx(2.4)
    assert rows[-1]['top_reaction_y_n'] == -114.
    assert rows[-1]['bottom_reaction_y_n'] == 120.
    assert rows[-1]['ALLKE'] == 2.
    assert assessment['completed'] is True
    assert assessment['quasi_static_validated'] is False
    assert assessment['experimental_agreement'] is False
    target = tmp_path / 'curve.csv'
    export.write_csv(str(target), rows)
    with target.open() as stream:
        parsed = list(csv.DictReader(stream))
    assert float(parsed[-1]['stress_mpa']) == pytest.approx(.06)


def test_positive_top_reaction_remains_negative_stress():
    odb = mock_odb()
    for name in ('Node PLATEN.1', 'Node PLATEN.2'):
        odb.steps['COMPRESS'].historyRegions[name].historyOutputs['RF2'] = output([0., 1., 2.])
    rows, _ = export.extract(odb)
    assert rows[-1]['stress_mpa'] < 0


def test_mismatched_history_grids_are_rejected_instead_of_truncated():
    odb = mock_odb()
    odb.steps['COMPRESS'].historyRegions['Node BASE.1'].historyOutputs['RF2'] = output(
        [0., 60., 120.], times=(0., .004, .01))
    with pytest.raises(ValueError, match='Mismatched history time grids'):
        export.extract(odb)


def test_missing_loading_node_or_energy_history_is_rejected():
    odb = mock_odb()
    del odb.steps['COMPRESS'].historyRegions['Node PLATEN.2']
    with pytest.raises(ValueError, match='Missing nodal history'):
        export.extract(odb)
    odb = mock_odb()
    del odb.steps['COMPRESS'].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLAE']
    with pytest.raises(ValueError, match='whole-model region'):
        export.extract(odb)


def test_partial_last_frame_is_preserved_without_completion_claim():
    odb = mock_odb()
    odb.steps['COMPRESS'].frames[-1].frameValue = .009
    _, assessment = export.extract(odb)
    assert assessment['completed'] is False
    assert assessment['last_frame_time_s'] == .009


def test_point_identity_and_flat_node_set_supported_without_label_collision():
    odb = mock_odb()
    top = odb.rootAssembly.nodeSets['TOP_LOAD']
    top.nodes = [NS(label=1, instanceName='PLATEN'), NS(label=2, instanceName='PLATEN')]
    top.instanceNames = ()
    regions = odb.steps['COMPRESS'].historyRegions
    region = regions.pop('Node PLATEN.1')
    region.point = NS(node=NS(label=1), instance=NS(name='PLATEN'))
    regions['Custom top point'] = region
    rows, _ = export.extract(odb)
    assert rows[-1]['top_reaction_y_n'] == -114.


def test_nonfinite_unsorted_or_duplicate_regions_are_rejected():
    odb = mock_odb()
    regions = odb.steps['COMPRESS'].historyRegions
    regions['Node PLATEN.1'].historyOutputs['RF2'] = output([0., float('nan'), 2.])
    with pytest.raises(ValueError, match='finite samples'):
        export.extract(odb)
    odb = mock_odb()
    regions = odb.steps['COMPRESS'].historyRegions
    duplicate = copy.deepcopy(regions['Node PLATEN.1'])
    duplicate.point = NS(node=NS(label=1), instance=NS(name='PLATEN'))
    regions['duplicate'] = duplicate
    with pytest.raises(ValueError, match='Multiple history regions'):
        export.extract(odb)
