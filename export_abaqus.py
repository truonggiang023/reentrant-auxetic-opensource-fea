"""Export two self-contained t=.3 Abaqus/Explicit comparison candidates."""
import argparse
import hashlib
import json
from pathlib import Path

from reentrant_fea.abaqus import retained_core, write_comparison
from reentrant_fea.plastic import load_plastic


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'abaqus')
    parser.add_argument('--mesh-archive', type=Path, default=ROOT/'validation/opencourant-thin-free-x-inputs.tar.gz')
    parser.add_argument('--config', type=Path, default=ROOT/'configs/dong2019_brick_damage.json')
    parser.add_argument('--plastic', type=Path, default=ROOT/'validation/materials/user_plastic_table_2026-10-26.inp')
    parser.add_argument('--duration-seconds', type=float, default=.01)
    parser.add_argument('--loading', choices=['smooth', 'linear'], default='smooth')
    parser.add_argument('--platen-mesh-size', type=float, default=2.5)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    if abs(config['wall']-.3) > 1e-12:
        parser.error('This retained core is t=.3mm; a config alone cannot change its mesh')
    nodes, elements, provenance = retained_core(args.mesh_archive)
    plastic = load_plastic(args.plastic)
    outputs = [args.output_dir/('thin_t03_' + suffix + '.inp') for suffix in ('plastic', 'damage_hypothesis')]
    if any(p.exists() or p.with_suffix('.json').exists() for p in outputs):
        raise FileExistsError('Use a fresh --output-dir; prior inputs are preserved')
    for damage, path in zip((False, True), outputs):
        metadata = write_comparison(path, nodes, elements, config, plastic, duration=args.duration_seconds,
                                    damage=damage, loading=args.loading, platen_size=args.platen_mesh_size)
        metadata.update(mesh_source=provenance,
                        config_source_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
                        hardening_source_sha256=hashlib.sha256(args.plastic.read_bytes()).hexdigest(),
                        generator_sources_sha256={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                                  for name in ('export_abaqus.py', 'src/reentrant_fea/abaqus.py', 'src/reentrant_fea/plastic.py')})
        path.with_suffix('.json').write_text(json.dumps(metadata, indent=2)+'\n')
        print(f'{path}: {metadata["core_mesh"]["elements"]} core + {metadata["steel_elements"]} steel C3D8R elements; damage={damage}; NOT Abaqus datachecked')


if __name__ == '__main__':
    main()
