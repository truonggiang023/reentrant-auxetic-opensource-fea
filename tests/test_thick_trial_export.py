"""Publishing must reject partial runs and avoid following evidence symlinks."""
import json
from pathlib import Path

import pytest

from scripts import run_thick_trial as mod


def inputs(tmp_path):
    prefix = tmp_path/'native'
    config = tmp_path/'config.json'
    config.write_text('{}')
    metadata = dict(prefix='native', wall_mm=1., source_config_sha256=mod.sha(config))
    Path(str(prefix)+'.radioss.json').write_text(json.dumps(metadata))
    return prefix, config


def test_unfinished_receipt_cannot_publish_a_report(tmp_path, monkeypatch):
    prefix, config = inputs(tmp_path)
    monkeypatch.setattr(mod, 'assess', lambda *a, **kw: dict(
        completed=False, run_binding_verified=True, metadata_binding_verified=True))
    target = tmp_path/'publication'
    with pytest.raises(ValueError, match='Only completed receipt-bound'):
        mod.export_completed(prefix, config, target, 'partial', 'a'*40)
    assert not target.exists()


def test_expected_filename_symlink_cannot_export_unrelated_contents(tmp_path, monkeypatch):
    prefix, config = inputs(tmp_path)
    for suffix in mod.SOLVER_SUFFIXES:
        path = Path(str(prefix)+suffix)
        if not path.exists():
            path.write_text('solver fixture')
    unrelated = tmp_path/'unrelated.txt'
    unrelated.write_text('must not enter exported evidence')
    history = Path(str(prefix)+'T01')
    history.unlink()
    history.symlink_to(unrelated)
    monkeypatch.setattr(mod, 'assess', lambda *a, **kw: dict(
        completed=True, run_binding_verified=True, metadata_binding_verified=True))
    target = tmp_path/'publication'
    with pytest.raises(ValueError, match='regular file'):
        mod.export_completed(prefix, config, target, 'symlink', 'a'*40)
    assert not target.exists()
