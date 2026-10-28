"""Physical input checks; these do not replace a licensed Abaqus datacheck."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile

import numpy as np
import pytest

from reentrant_fea.abaqus import mesh_audit, retained_core, write_comparison
from reentrant_fea.plastic import load_plastic


ROOT = Path(__file__).resolve().parents[1]
PLASTIC = ROOT / 'validation/materials/user_plastic_table_2026-10-26.inp'


@pytest.fixture
def two_bricks():
    # Two unit cubes sharing one complete X face, not coincident duplicate nodes.
    nodes = {}
    ids = {}
    for z in (0, 1):
        for y in (0, 1):
            for x in (0, 1, 2):
                ids[x, y, z] = len(nodes) + 1
                nodes[ids[x, y, z]] = [float(x), float(y), float(z)]
    elements = {}
    for x in (0, 1):
        elements[x + 1] = [ids[x, 0, 0], ids[x + 1, 0, 0],
                           ids[x + 1, 1, 0], ids[x, 1, 0],
                           ids[x, 0, 1], ids[x + 1, 0, 1],
                           ids[x + 1, 1, 1], ids[x, 1, 1]]
    return nodes, elements


@pytest.fixture
def config():
    return json.loads((ROOT / 'configs/dong2019_brick_damage.json').read_text())


def blocks(text):
    """Read keyword blocks for independent numerical/geometry assertions."""
    result = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            result.append((line.upper(), []))
        else:
            result[-1][1].append([value.strip() for value in line.split(',')])
    return result


def named_block(data, keyword, name):
    return next(rows for header, rows in data
                if header.startswith(keyword) and name.upper() in header)


def test_mesh_checks_geometry_and_inverted_orientation(two_bricks):
    nodes, elements = two_bricks
    audit = mesh_audit(nodes, elements)
    assert audit['volume_mm3'] == pytest.approx(2.)
    assert audit['dimensions_mm'] == [2., 1., 1.]
    assert audit['minimum_sampled_jacobian_mm3'] == pytest.approx(.125)
    assert audit['interior_facets'] == 1
    assert audit['exterior_facets'] == 10
    reflected = {ident: [-p[0], p[1], p[2]] for ident, p in nodes.items()}
    with pytest.raises(ValueError, match='Inverted'):
        mesh_audit(reflected, elements)


def test_steel_contact_faces_point_into_core(tmp_path, two_bricks, config):
    nodes, elements = two_bricks
    path = tmp_path / 'candidate.inp'
    meta = write_comparison(path, nodes, elements, config,
                            load_plastic(PLASTIC), platen_size=2.)
    data = blocks(path.read_text())
    deck_nodes = {int(row[0]): np.array(list(map(float, row[1:])))
                  for row in named_block(data, '*NODE', '')}
    deck_elements = {int(row[0]): list(map(int, row[1:]))
                     for header, rows in data if header.startswith('*ELEMENT,')
                     for row in rows}
    # Abaqus C3D8 S3=(1,2,6,5), S5=(3,4,8,7), both outward ordered.
    face_indices = {'S3': [0, 1, 5, 4], 'S5': [2, 3, 7, 6]}
    for name, face, y, direction in [('BOTTOM_INNER', 'S5', 0., 1),
                                     ('TOP_INNER', 'S3', 1., -1)]:
        for eid, label in named_block(data, '*SURFACE,', 'NAME=' + name):
            assert label == face
            points = np.array([deck_nodes[deck_elements[int(eid)][i]]
                               for i in face_indices[label]])
            np.testing.assert_allclose(points[:, 1], y)
            normal = np.cross(points[1] - points[0], points[2] - points[1])
            assert direction * normal[1] > 0
            np.testing.assert_allclose(normal[[0, 2]], 0.)
    assert meta['core_mesh']['interior_facets'] == 1
    assert any('NAME=CORE_INTERIOR' in header for header, _ in data)
    assert any('NAME=CORE_CONTACT' in header for header, _ in data)
    includes = named_block(data, '*CONTACT INCLUSIONS', '')
    assert includes == [['CORE_CONTACT', 'CORE_CONTACT'],
                        ['CORE_CONTACT', 'BOTTOM_INNER'],
                        ['CORE_CONTACT', 'TOP_INNER']]


def test_confimed_material_and_free_core_x_are_preserved(tmp_path, two_bricks, config):
    nodes, elements = two_bricks
    path = tmp_path / 'material.inp'
    meta = write_comparison(path, nodes, elements, config, load_plastic(PLASTIC))
    data = blocks(path.read_text())
    # Exact supplied numbers, no engineering-to-true conversion or fitted tail.
    np.testing.assert_array_equal(load_plastic(path), load_plastic(PLASTIC))
    assert len(load_plastic(path)) == 178
    material = next(i for i, (header, _) in enumerate(data)
                    if header == '*MATERIAL, NAME=CONFIRMED_ALSI10MG')
    assert list(map(float, data[material + 1][1][0])) == [2.67e-9]
    assert list(map(float, data[material + 2][1][0])) == [69000., .3]
    constraints = [row for header, rows in data if header.startswith('*BOUNDARY') for row in rows]
    assert ['CORE_NODES', '3', '3', '0.'] in constraints
    assert [row for row in constraints if row[0].startswith('CORE')] == [
        ['CORE_NODES', '3', '3', '0.']]
    assert not any('MASS SCALING' in header for header, _ in data)
    assert not any(header.startswith('*DAMAGE') for header, _ in data)
    assert meta['mass_scaling'] is False
    assert meta['damage_enabled'] is False
    assert meta['abaqus_datacheck_performed'] is False


def test_native_damage_uses_strain_then_triaxiality_without_ccx_extras(tmp_path, two_bricks, config):
    nodes, elements = two_bricks
    path = tmp_path / 'hypothesis.inp'
    meta = write_comparison(path, nodes, elements, config, load_plastic(PLASTIC), damage=True)
    data = blocks(path.read_text())
    rows = named_block(data, '*DAMAGE INITIATION', '')
    expected = [[eps, eta, 0.] for eta, eps in config['ductile_damage']['fracture_locus']]
    np.testing.assert_allclose(np.array(rows, float), expected, rtol=0, atol=0)
    assert named_block(data, '*DAMAGE EVOLUTION', '') == [['0.001']]
    assert any('ELEMENT DELETION=YES' in header for header, _ in data)
    assert not any(header.startswith(('*USER MATERIAL', '*DEPVAR', '*DAMAGE STABILIZATION'))
                   for header, _ in data)
    assert meta['damage_enabled'] is True
    assert meta['calibrated'] is False
    assert meta['experimental_agreement'] is False


@pytest.mark.parametrize('loading,duration,peak_factor', [('smooth', .01, 1.875),
                                                        ('linear', 19.5, 1.)])
def test_loading_is_explicit_and_units_consistent(tmp_path, two_bricks, config,
                                                  loading, duration, peak_factor):
    nodes, elements = two_bricks
    path = tmp_path / 'loading.inp'
    meta = write_comparison(path, nodes, elements, config, load_plastic(PLASTIC),
                            duration=duration, loading=loading)
    data = blocks(path.read_text())
    assert named_block(data, '*DYNAMIC, EXPLICIT', '') == [['', str(duration)]]
    boundary = named_block(data, '*BOUNDARY, AMPLITUDE=', '')
    assert float(boundary[0][-1]) == pytest.approx(-.65)
    amplitude = named_block(data, '*AMPLITUDE', '')
    np.testing.assert_allclose(np.array(amplitude, float), [[0., 0., duration, 1.]])
    assert meta['target_compression_mm'] == pytest.approx(.65)
    assert meta['constant_2_mm_per_min_duration_seconds'] == pytest.approx(19.5)
    assert meta['imposed_peak_speed_mm_per_s'] == pytest.approx(.65 / duration * peak_factor)
    if loading == 'linear':
        assert meta['imposed_mean_speed_mm_per_s'] * 60 == pytest.approx(2.)


def test_existing_inputs_or_sidecars_are_never_overwritten(tmp_path, two_bricks, config):
    nodes, elements = two_bricks
    path = tmp_path / 'keep.inp'
    write_comparison(path, nodes, elements, config, load_plastic(PLASTIC))
    before = {p: p.read_bytes() for p in [path, path.with_suffix('.json')]}
    with pytest.raises(FileExistsError):
        write_comparison(path, nodes, elements, config, load_plastic(PLASTIC), damage=True)
    assert all(p.read_bytes() == content for p, content in before.items())
    path.unlink()
    with pytest.raises(FileExistsError):
        write_comparison(path, nodes, elements, config, load_plastic(PLASTIC))
    assert not path.exists()
    assert path.with_suffix('.json').read_bytes() == before[path.with_suffix('.json')]


def test_committed_retained_mesh_identity_and_counts():
    archive = ROOT / 'validation/opencourant-thin-free-x-inputs.tar.gz'
    nodes, elements, provenance = retained_core(archive)
    assert provenance['archive_sha256'] == '5c85034a893c0bab619efb10b663e1b9c7c593166ce0fb3299f641ab4e0ff79c'
    assert provenance['starter_sha256'] == '3e2d6fc4d7ac3909046db1d597cf5e26b8d0217723a7956afc8775ab7c1a82e9'
    assert len(nodes) == 21870
    assert len(elements) == 7392
    assert max(nodes) == 21870  # The appended native steel nodes are excluded.
    np.testing.assert_allclose(np.ptp(np.array(list(nodes.values())), axis=0),
                               [76.497146699265, 48., 25.], rtol=0, atol=1e-10)


def test_cli_creates_two_standalone_inputs_and_binds_sources(tmp_path, two_bricks):
    nodes, elements = two_bricks
    lines = ['/NODE']
    lines += [' '.join(map(str, [ident] + point)) for ident, point in nodes.items()]
    lines += ['/BRICK/1']
    lines += [' '.join(map(str, [ident] + row)) for ident, row in elements.items()]
    lines += ['/END']
    raw = ('\n'.join(lines) + '\n').encode()
    archive = tmp_path / 'tiny-mesh.tar.gz'
    with tarfile.open(archive, 'w:gz') as stream:
        info = tarfile.TarInfo('tiny_0000.rad')
        info.size = len(raw)
        stream.addfile(info, io.BytesIO(raw))
    output = tmp_path / 'abaqus'
    command = [sys.executable, str(ROOT / 'export_abaqus.py'), '--mesh-archive', str(archive),
               '--output-dir', str(output), '--platen-mesh-size', '5']
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    for suffix, enabled in [('plastic', False), ('damage_hypothesis', True)]:
        path = output / ('thin_t03_' + suffix + '.inp')
        meta = json.loads(path.with_suffix('.json').read_text())
        assert meta['damage_enabled'] is enabled
        assert meta['mesh_source']['starter_sha256'] == hashlib.sha256(raw).hexdigest()
        assert meta['hardening_source_sha256'] == hashlib.sha256(PLASTIC.read_bytes()).hexdigest()
        assert meta['input_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
        assert not any(header.startswith('*INCLUDE') for header, _ in blocks(path.read_text()))
    retry = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert retry.returncode != 0
    assert 'preserved' in retry.stderr
