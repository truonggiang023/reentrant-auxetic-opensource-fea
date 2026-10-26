"""Bind generated decks and assessments to the actual geometry/material input."""
import hashlib
import json
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def record_generation(deck, config_path, plastic_path, *, config, metadata, source_paths=()):
    deck, config_path, plastic_path = map(Path, (deck, config_path, plastic_path))
    sidecar = deck.with_suffix('.material.dat')
    record = {
        'generation_config_sha256': sha256(config_path),
        'generation_config_path': str(config_path.resolve()),
        'config_wall_mm': float(config['wall']),
        'hardening_source_sha256': sha256(plastic_path),
        'input_sha256': sha256(deck),
        'material_sidecar_sha256': sha256(sidecar) if sidecar.is_file() else None,
        'generator_sources': {str(Path(path).resolve()): sha256(path) for path in source_paths},
        'metadata': metadata,
    }
    deck.with_suffix('.generation.json').write_text(json.dumps(record, indent=2)+'\n')
    return record


def read_generation(deck):
    """Return only matching generation evidence; legacy runs may have none."""
    deck = Path(deck)
    path = deck.with_suffix('.generation.json')
    if not path.is_file():
        return None
    record = json.loads(path.read_text())
    if record.get('input_sha256') != sha256(deck):
        raise ValueError('Generation record belongs to a different input deck')
    sidecar = deck.with_suffix('.material.dat')
    expected = record.get('material_sidecar_sha256')
    if (expected is not None) != sidecar.is_file() or (expected and sha256(sidecar) != expected):
        raise ValueError('Generation record belongs to a different material sidecar')
    return record
