"""Prepare bounded fracture-locus hypotheses and score actual solver assessments.

This is an inverse-fit workflow, not a recovered paper material definition.
Training uses t=1 mm; t=.3 mm remains an independent holdout. Run simulations
explicitly, assess accepted results, then rank only completed training curves.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from reentrant_fea.calibration import fracture_candidate, score_assessment, select_candidate
from reentrant_fea.plastic import load_plastic, validate_material, validate_plastic_table


ROOT = Path(__file__).resolve().parent
ROLE_INTERVALS = {'training': [.04, .45], 'holdout': [.03, .65]}
NATIVE_ROLE_INTERVALS = {'training': [.04, .45], 'holdout': [.03, .64]}
NATIVE_FORCE_COLUMN = 'Contact force and energy                 3 Contact 3                                var 64'
NATIVE_TANGENTIAL_COLUMN = 'Contact force and energy                 3 Contact 3                                var 67'


def native_force_mapping(*, legacy=False, reference_points=False):
    """Declare verified TYPE25 contact impulses or observed-RP reaction histories."""
    if reference_points:
        if legacy:
            raise ValueError('RP mapping has no legacy normal-only form')
        return dict(source='rp',channel_kind='impulse',strain_measure='observed_relative_RP_DY')
    result = dict(csv_column=NATIVE_FORCE_COLUMN, compressive_sign=-1, channel_kind='impulse')
    if not legacy:
        result.update(tangential_csv_column=NATIVE_TANGENTIAL_COLUMN, component_operation='FN-FT')
    return result


def legacy_native_force_mapping(mapping):
    """Validate the declared mapping; identify a historical normal-only study.

    Legacy input/solver files remain unchanged. Reassessment always uses both
    actual bound channels, and explicitly records the interpretation upgrade.
    """
    if mapping == native_force_mapping(reference_points=True):
        return False
    if mapping == native_force_mapping():
        return False
    if mapping == native_force_mapping(legacy=True):
        return True
    raise ValueError('Native study requires verified top FN/FT impulse channels, FN-FT and sign')


def native_assess_kwargs(mapping):
    legacy_native_force_mapping(mapping)
    if mapping == native_force_mapping(reference_points=True):
        return dict(force_source='rp')
    return dict(force_column=NATIVE_FORCE_COLUMN,tangential_column=NATIVE_TANGENTIAL_COLUMN,
                force_sign=-1,force_kind='impulse')


def dump(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + '\n')


def reference_curves(path):
    curves = {'experiment_1': [], 'experiment_2': []}
    with path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            if row['curve'] in curves:
                curves[row['curve']].append([float(row['strain']), float(row['stress_MPa'])])
    if not all(curves.values()):
        raise ValueError(f'Missing experimental replicates in {path}')
    return list(curves.values())


def prepare(base, output, *, backend='ccx', duration=None, ams_time_step=None, plastic=None,
            holdout_end_strain=.65):
    """Write configurations and exact follow-up commands without launching solves."""
    if not np.isfinite(holdout_end_strain) or not .65 <= holdout_end_strain < 1:
        raise ValueError('Holdout end strain must lie in [.65, 1)')
    if backend not in ('ccx', 'opencourant'):
        raise ValueError('Select ccx or opencourant backend')
    if backend == 'ccx' and duration is not None:
        raise ValueError('Physical loading duration applies to the OpenCourant backend only')
    if backend == 'opencourant' and (duration is None or not np.isfinite(duration) or duration <= 0):
        raise ValueError('OpenCourant requires an explicit finite positive duration in seconds')
    if ams_time_step is not None and (backend != 'opencourant'
            or not np.isfinite(ams_time_step) or not 0 < ams_time_step <= duration/100):
        raise ValueError('AMS requires OpenCourant and a positive step at most duration/100')
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use an empty study output directory to preserve existing simulations')
    source = json.loads(base.read_text())
    plastic_source = Path(plastic).resolve() if plastic is not None else ROOT/'reference_setup_numeric_source.py'
    # Validate both supported input formats before creating any study files.
    validate_material(source, load_plastic(plastic_source))
    plastic_sha256 = hashlib.sha256(plastic_source.read_bytes()).hexdigest()
    records = []
    variants = [(f'locus_{low:g}_{high:g}', low, high, {})
                for low in (-.2, 0., .2) for high in (.5, 1., 1.5)]
    # Numerical sensitivity variants are not interchangeable material fits.
    if backend == 'ccx':
        variants += [('viscosity_zero', 0., 1., {'viscosity_time': 0.}),
                     ('viscosity_half', 0., 1., {'viscosity_time': .0025}),
                     ('residual_tenth', 0., 1., {'residual_stiffness': 1e-5})]
    else:
        variants += [('duration_double', 0., 1., {'duration_seconds': 2*duration})]
        if ams_time_step is not None:
            variants += [('ams_half', 0., 1., {'ams_time_step': ams_time_step/2})]
    for candidate, low, high, nuisance in variants:
        for wall, role in ((1., 'training'), (.3, 'holdout')):
            config = fracture_candidate(source, low, high)
            config['wall'] = wall
            config['compression_strain'] = .5 if role == 'training' else float(holdout_end_strain)
            if backend == 'ccx':
                config['ductile_damage'].update(nuisance)
            name = f'{candidate}_{role}'
            config_path = output/'configs'/f'{name}.json'
            dump(config_path, config)
            deck = output/'runs'/f'{name}.inp'
            assessment = output/'assessments'/f'{name}.json'
            record = {
                'candidate': candidate, 'role': role, 'wall_mm': wall,
                'backend': backend,
                'estimated_eta_low': low, 'estimated_eta_high': high,
                'numerical_sensitivity': bool(nuisance), 'numerical_overrides': nuisance,
                'config': str(config_path),
                'config_sha256': hashlib.sha256(config_path.read_bytes()).hexdigest(),
                'deck': str(deck),
                'hardening_source': str(plastic_source),
                'hardening_source_sha256': plastic_sha256,
                'assessment': str(assessment),
                'solve_argv': ['python', 'run_paper.py', '--config', str(config_path),
                               '--plastic', str(plastic_source),
                               '--output', str(deck), '--target-output-only', '--self-contact',
                               '--threads', '2', '--run'],
                'assess_argv': ['python', 'assess_run.py', str(deck), '--output', str(assessment)],
            }
            if backend == 'opencourant':
                # Native model names accept ASCII letters/digits/_/- only;
                # preserve decimal parameter values in candidate metadata.
                prefix = output/'runs'/name.replace('.', 'p')
                metadata = prefix.with_name(prefix.name+'.radioss.json')
                history = prefix.with_name(prefix.name+'T01.csv')
                actual_duration = nuisance.get('duration_seconds', duration)
                actual_ams = nuisance.get('ams_time_step', ams_time_step)
                record.update(deck=str(prefix.with_name(prefix.name+'_0000.rad')),
                              metadata=str(metadata), history=str(history),
                              assessment=str(prefix.with_name(prefix.name+'.assessment.json')),
                              duration_seconds=actual_duration, ams_time_step=actual_ams,
                              force_mapping=native_force_mapping(reference_points=bool(config.get('platen_reference_points',False))))
                record['solve_argv'] = ['python', 'run_radioss.py', '--config', str(config_path),
                                       '--plastic', record['hardening_source'], '--output', str(prefix),
                                       '--duration-seconds', str(actual_duration), '--drop-ccx-viscosity', '--run']
                if actual_ams is not None:
                    record['solve_argv'] += ['--ams-time-step', str(actual_ams)]
                record['assess_argv'] = ['python', 'assess_courant.py', str(metadata), str(history),
                                        '--force-column', NATIVE_FORCE_COLUMN, '--force-sign', '-1',
                                        '--tangential-column', NATIVE_TANGENTIAL_COLUMN,
                                        '--force-kind', 'impulse', '--output', record['assessment']]
                if config.get('platen_reference_points'):
                    record['assess_argv'] = ['python','assess_courant.py',str(metadata),str(history),
                                            '--force-source','rp','--output',record['assessment']]
            records.append(record)
    manifest = {
        'status': 'unrun hypotheses; no calibrated material', 'calibrated': False,
        'base_config': str(base.resolve()), 'working_directory': str(ROOT), 'backend': backend,
        'hardening_source': str(plastic_source), 'hardening_source_sha256': plastic_sha256,
        'fit_intervals': deepcopy_intervals(),
        'holdout_end_strain': float(holdout_end_strain),
        'native_force_source': ('observed RP REACY and relative DY' if source.get('platen_reference_points')
                                else 'TYPE25 total contact FN-FT') if backend=='opencourant' else None,
        'deformation_targets': {'training': .5, 'holdout': float(holdout_end_strain)},
        'late_response_scope': 'End strain controls run coverage; late densification is assessed separately from the historical fit interval. No extrapolation or altered training/holdout role.',
        'geometry_qualification': {
            'qualified': False,
            'reason': 'Inherited centerline geometry with wall changed; verify dimensions, relative density and mesh against each paper specimen before material fitting',
        },
        'fixed_paper_inputs': {'fracture_point': [.333, .047], 'failure_displacement_mm': .001},
        'estimated_inputs': ['eta for fracture strain .08', 'eta for fracture strain .018'],
        'notes': [
            'Bounds are investigator assumptions, not measured confidence intervals.',
            'Fit completed thick-wall runs first. Transfer chosen parameters unchanged to thin-wall holdout.',
            'Backend-specific numerical sensitivity variants test robustness and are excluded from material-candidate selection.',
            'A low curve error cannot compensate for invalid deletion/contact or establish a unique fracture locus.',
            'Changing wall alone changes outer dimensions and relative density; thick-wall hypotheses need geometry qualification before candidate selection.',
            'Thick-wall reference Fig12(a) ends at strain .5; training fit interval .04–.45. Thin run end is declared separately in deformation_targets; fit intervals do not define full densification coverage.',
            'The paper reports tensile fracture strain near .047 and uses .047 in ductile initiation; do not infer that its plotted total strain is a measured equivalent plastic-strain history.',
        ], 'runs': records,
    }
    if backend == 'opencourant':
        manifest['runtime_setup'] = 'source /workspace/opencourant/env.sh before running solve_argv'
        manifest['fit_intervals'] = {role: values.copy() for role, values in NATIVE_ROLE_INTERVALS.items()}
        manifest['backend_differences'] = [
            'Native INIEVO erosion and all-face TYPE25 contact replace the CCX residual-element material.',
            'CCX pseudo-time viscosity is intentionally omitted; it is not converted to physical seconds.',
            'CCX compression cutoff is not transferred. Native endpoint extrapolation is a constant boundary convention, not additional measured data.',
            'Native characteristic length, solid formulation and initiation-increment integration differ from CCX.',
            'Accelerated loading and AMS require actual-core energy, speed and AMS response sensitivity; a passing global screen is insufficient.',
            'TYPE25 IMS=2 top-secondary force is FN-FT: FNY var64 and FTY var67 are cumulative impulses. Compression is -Delta(FNY-FTY)/Delta(time), not the normal component alone.',
            'Native holdout force fit ends at .64; interval-average force samples precede the configured deformation end and are never extrapolated to it.',
        ]
        manifest['native_qualification'] = {'qualified': False,
            'reason': 'Actual full-lattice core energy, loading-speed and AMS response sensitivity remain unresolved'}
    dump(output/'manifest.json', manifest)
    return manifest


def score_native_record(row, config, saved_assessment, experiments, interval):
    """Reassess actual native files; cached JSON flags/curves cannot qualify a run."""
    from assess_courant import assess
    metadata_path, history = Path(row['metadata']), Path(row['history'])
    metadata = json.loads(metadata_path.read_text())
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    source = Path(row.get('hardening_source', ROOT/'reference_setup_numeric_source.py'))
    validate_plastic_table(load_plastic(source))
    source_hash = digest(source)
    legacy_mapping = legacy_native_force_mapping(row.get('force_mapping'))
    if (saved_assessment.get('metadata_sha256') != digest(metadata_path)
            or saved_assessment.get('history_sha256') != digest(history)):
        raise ValueError('Native assessment belongs to different metadata or history')
    if (metadata.get('source_config_sha256') != row['config_sha256']
            or metadata.get('wall_mm') != row['wall_mm']
            or metadata.get('hardening_source_sha256') != source_hash
            or row.get('hardening_source_sha256', source_hash) != source_hash
            or Path(metadata['starter']).resolve() != Path(row['deck']).resolve()):
        raise ValueError('Native generation configuration, wall, hardening or deck does not match candidate')
    law = config['ductile_damage']
    if (not np.allclose(metadata.get('fracture_locus', []), law['fracture_locus'], rtol=1e-10, atol=1e-12)
            or metadata.get('duration_seconds') != row['duration_seconds']
            or metadata.get('ams', {}).get('requested_time_step_seconds') != row['ams_time_step']
            or metadata.get('ccx_viscosity_omitted') != law.get('viscosity_time', 0)):
        raise ValueError('Native material assumptions, duration or AMS do not match candidate')
    # assess() verifies generated deck hashes and the .run.json binding to the
    # actual starter/engine/T01/CSV files, and then reconstructs the force curve.
    actual = assess(metadata_path, history, **native_assess_kwargs(row.get('force_mapping')))
    result = score_assessment(actual, experiments, interval)
    result.update(run_binding_verified=actual['run_binding_verified'],
                  metadata_binding_verified=actual.get('metadata_binding_verified', False),
                  metadata_sha256=actual['metadata_sha256'], history_sha256=actual['history_sha256'],
                  generation_config_matches=True, force_mapping=actual['force_mapping'],
                  legacy_manifest_force_mapping=legacy_mapping,
                  force_mapping_recomputed_from_bound_channels=True,
                  global_energy_screen=actual['global_energy_screen'],
                  core_energy_screen=actual['core_energy_screen'],
                  duration_seconds=metadata['duration_seconds'], ams=metadata.get('ams'),
                  deformation_target_nominal_strain=metadata['target_nominal_strain'],
                  last_observed_nominal_strain=actual.get('last_observed_nominal_strain'),
                  coverage_and_completion_passed=result['eligible_for_selection'])
    # A full explicit run and a global energy screen do not demonstrate a
    # quasi-static material response. No JSON "qualified" flag bypasses this
    # unresolved actual-core/speed/AMS response study.
    blockers = ['actual full-lattice core energy and work balance not qualified',
                'loading-speed response sensitivity not qualified',
                'mesh/contact/erosion response sensitivity not qualified']
    if not actual.get('metadata_binding_verified'):
        blockers.append('execution-time metadata hash was not recorded')
    if metadata.get('ams', {}).get('enabled'):
        blockers.append('actual full-lattice AMS response sensitivity not qualified')
    if not actual['core_energy_screen']['sampled_energy_screen_passed']:
        blockers += actual['core_energy_screen']['reasons']
    result['qualification_blockers'] = blockers
    result['eligible_for_selection'] = False
    return result


def diagnostic_training_ranking(evaluated, geometry):
    """Rank only complete bound training curves, without qualifying a material."""
    ranked = []
    for row in evaluated:
        result = row['score']
        bound = (result.get('run_binding_verified') and result.get('metadata_binding_verified')
                 if row.get('backend') == 'opencourant' else result.get('execution_binding_verified'))
        if (row['role'] != 'training' or not np.isclose(row['wall_mm'], 1.) or row['numerical_sensitivity']
                or not result.get('solver_completed') or not result.get('full_interval_covered')
                or not result.get('generation_config_matches') or not bound):
            continue
        blockers = list(result.get('qualification_blockers', []))
        if not geometry.get('qualified'):
            blockers.append('specimen geometry has not been qualified against paper')
        blockers.append('diagnostic curve ranking does not establish a unique calibration or independent holdout agreement')
        ranked.append(dict(candidate=row['candidate'], wall_mm=row['wall_mm'],
            overlap_rmse_MPa=result['overlap_rmse_MPa'], interval_scored=result['interval_scored'],
            coverage_fraction=result['coverage_fraction'], eligible_for_selection=False,
            qualification_blockers=blockers))
    ranked.sort(key=lambda row: (row['overlap_rmse_MPa'], row['candidate']))
    return [dict(rank=i+1, **row) for i, row in enumerate(ranked)]


def score(manifest, reference_dir):
    references = {role: reference_curves(reference_dir/f'dong2019_fig12_{kind}.csv')
                  for role, kind in (('training', 'thick'), ('holdout', 'thin'))}
    evaluated, unrun = [], []
    for row in manifest['runs']:
        path = Path(row['assessment'])
        if not path.is_file():
            unrun.append({'candidate': row['candidate'], 'role': row['role'], 'reason': 'assessment absent'})
            continue
        config_path = Path(row['config'])
        if hashlib.sha256(config_path.read_bytes()).hexdigest() != row['config_sha256']:
            raise ValueError(f'Configuration changed since study creation: {config_path}')
        assessment = json.loads(path.read_text())
        config = json.loads(config_path.read_text())
        backend = row.get('backend', manifest.get('backend', 'ccx'))
        if backend == 'opencourant':
            item = {key: row[key] for key in ('candidate', 'role', 'wall_mm', 'numerical_sensitivity')}
            item['backend'] = backend
            interval = manifest.get('fit_intervals', NATIVE_ROLE_INTERVALS)[row['role']]
            item['score'] = score_native_record(row, config, assessment, references[row['role']], interval)
            evaluated.append(item)
            continue
        # Do not rank a stale assessment or attach a different material to the
        # candidate merely because its filename happens to match the manifest.
        deck = Path(row['deck'])
        material = deck.with_suffix('.material.dat')
        if (not deck.is_file() or not material.is_file()
                or hashlib.sha256(deck.read_bytes()).hexdigest() != assessment.get('input_sha256')
                or hashlib.sha256(material.read_bytes()).hexdigest() != assessment.get('material_sidecar_sha256')):
            raise ValueError(f'Missing or mismatched deck/material evidence for {row["candidate"]}')
        header, *table = material.read_text().splitlines()
        magic, nh, nf = header.split()
        values = np.array([[float(v) for v in line.split()] for line in table if line.strip()])
        if (magic != 'REENTRANT_TABLES_V1' or values.shape != (int(nh)+int(nf), 2)
                or int(nf) != len(config['ductile_damage']['fracture_locus'])
                or not np.allclose(values[-int(nf):], config['ductile_damage']['fracture_locus'], rtol=1e-10, atol=1e-12)):
            raise ValueError(f'Fracture locus does not match candidate {row["candidate"]}')
        plastic_source = Path(row.get('hardening_source', ROOT/'reference_setup_numeric_source.py'))
        source_sha256 = hashlib.sha256(plastic_source.read_bytes()).hexdigest()
        if row.get('hardening_source_sha256', source_sha256) != source_sha256:
            raise ValueError(f'Hardening source changed for candidate {row["candidate"]}')
        hardening = validate_plastic_table(load_plastic(plastic_source))
        if hardening.shape != (int(nh), 2) or not np.allclose(values[:int(nh)], hardening, rtol=1e-10, atol=1e-12):
            raise ValueError(f'Hardening table does not match candidate {row["candidate"]}')
        lines = deck.read_text().splitlines()
        material_count = 0
        for index, line in enumerate(lines):
            if not line.startswith('*USER MATERIAL'):
                continue
            material_count += 1
            props = []
            for data in lines[index+1:]:
                if data.startswith('*'):
                    break
                props.extend(float(value.strip()) for value in data.split(',') if value.strip())
            law = config['ductile_damage']
            expected = [config['young'], config['poisson'], law['failure_displacement'],
                        law.get('residual_stiffness', 1e-4), law.get('compression_cutoff', -1/3),
                        law.get('viscosity_time', 0.)]
            if len(props) != 9 or not np.allclose(np.array(props)[[0, 1, 2, 3, 7, 8]], expected):
                raise ValueError(f'Material descriptors do not match candidate {row["candidate"]}')
        if not material_count:
            raise ValueError(f'No native material in candidate deck {row["candidate"]}')
        item = {key: row[key] for key in ('candidate', 'role', 'wall_mm', 'numerical_sensitivity')}
        item['backend'] = 'ccx'
        try:
            interval = manifest.get('fit_intervals', ROLE_INTERVALS)[row['role']]
            receipt = path.with_name(path.name+'.execution.json')
            execution_bound = False
            if receipt.is_file():
                # The runner's receipt binds raw CCX outputs; a cached JSON
                # completed flag and matching input deck alone are insufficient.
                from run_calibration import verify_completed
                assessment = verify_completed(row)
                execution_bound = True
            item['score'] = score_assessment(assessment, references[row['role']], interval)
            item['score']['execution_binding_verified'] = execution_bound
            binding = (assessment.get('generation_config_sha256') == row['config_sha256']
                       and assessment.get('config_wall_mm') == row['wall_mm']
                       and assessment.get('hardening_source_sha256') == source_sha256)
            item['score']['generation_config_matches'] = binding
            binding_blockers = []
            if not execution_bound:
                binding_blockers.append('CCX execution receipt is absent; cached assessment is diagnostic only')
            if not binding:
                binding_blockers.append('assessment lacks matching generation-config, wall or hardening-source provenance')
            if binding_blockers:
                item['score']['eligible_for_selection'] = False
                item['score']['binding_blocker'] = '; '.join(binding_blockers)
        except ValueError as error:
            item['score'] = {'eligible_for_selection': False, 'reason': str(error)}
        evaluated.append(item)
    training = [row for row in evaluated if row['role'] == 'training' and not row['numerical_sensitivity']]
    geometry = manifest.get('geometry_qualification', {'qualified': False, 'reason': 'geometry not qualified'})
    if geometry.get('qualified') is True:
        selection = select_candidate(training)
    else:
        selection = {'selected_candidate': None, 'calibrated': False,
                     'status': 'blocked: specimen geometry has not been qualified against paper'}
    held_out = [row for row in evaluated if row['role'] == 'holdout'
                and row['candidate'] == selection['selected_candidate']]
    return {'calibrated': False, 'selection': selection, 'thin_wall_holdout': held_out,
            'diagnostic_training_ranking': diagnostic_training_ranking(evaluated, geometry),
            'fit_intervals': manifest.get('fit_intervals', NATIVE_ROLE_INTERVALS
                                           if manifest.get('backend') == 'opencourant' else deepcopy_intervals()),
            'legacy_interval_note': ('Legacy single fit_interval ignored: paper thick/thin curves have different strain extents'
                                     if 'fit_interval' in manifest else None),
            'geometry_qualification': geometry,
            'evaluated': evaluated, 'unrun': unrun,
            'identifiability': 'not established; evaluate response sensitivities and numerical convergence',
            'experimental_agreement': 'not established by candidate ranking alone',
            'reference_uncertainty': 'digitized paper curves; see dong2019_fig12_provenance.json',
            'assessment_binding': 'CCX deck/material/provenance checked; native metadata/deck/execution/history reassessed from actual files'}


def deepcopy_intervals():
    return {role: interval.copy() for role, interval in ROLE_INTERVALS.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    generate = commands.add_parser('prepare')
    generate.add_argument('--config', type=Path, default=ROOT/'configs/dong2019_brick_damage.json')
    generate.add_argument('--output', type=Path, required=True)
    generate.add_argument('--plastic', type=Path, default=ROOT/'reference_setup_numeric_source.py',
                          help='True stress/plastic strain: Abaqus *Plastic block or literal PLASTIC_DATA Python source; neither is executed')
    generate.add_argument('--backend', choices=('ccx', 'opencourant'), default='ccx')
    generate.add_argument('--duration', '--duration-seconds', dest='duration', type=float,
                          help='Explicit loading duration in seconds for native OpenCourant studies')
    generate.add_argument('--ams-time-step', type=float,
                          help='Native AMS step; changes inertial response and requires sensitivity checks')
    generate.add_argument('--holdout-end-strain', type=float, default=.65,
                          help='End strain for thin holdout; use .78 to cover the visible late Fig12(b) rise. Historical fit interval is unchanged.')
    rank = commands.add_parser('score')
    rank.add_argument('manifest', type=Path)
    rank.add_argument('--reference-dir', type=Path, default=ROOT/'validation/paper_reference')
    rank.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'prepare':
        result = prepare(args.config, args.output, backend=args.backend,
                         duration=args.duration, ams_time_step=args.ams_time_step, plastic=args.plastic,
                         holdout_end_strain=args.holdout_end_strain)
        print(f"Prepared {len(result['runs'])} configurations; simulations have not run.")
    else:
        result = score(json.loads(args.manifest.read_text()), args.reference_dir)
        dump(args.output, result)
        print(json.dumps(result['selection'], indent=2))


if __name__ == '__main__':
    main()
