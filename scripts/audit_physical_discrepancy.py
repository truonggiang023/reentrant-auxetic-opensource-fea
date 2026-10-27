#!/usr/bin/env python3
"""Reproduce read-only geometry, early stiffness and coupon-raster diagnostics.

This audit cannot calibrate a fracture law or certify a quasi-static run. Input
paper images remain outside the repository. Coupon pixels give a diagnostic
slope estimate only; their strain/stress definitions and toe compliance are not
established by the figure alone.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image

from reentrant_fea.geometry import core_polygon, parameters
from reentrant_fea.radioss_results import read_history_csv
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from plot_courant_frames import convert, read_frame, strain_at_time


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def curve_diagnostic(values, interval=(.005, .02)):
    data = np.asarray(values, dtype=float)
    selected = data[(data[:, 0] >= interval[0]) & (data[:, 0] <= interval[1])]
    if (len(selected) < 3 or not np.isfinite(data).all()
            or np.any(np.diff(data[:, 0]) <= 0)):
        raise ValueError('Insufficient finite ordered curve samples in the diagnostic interval')
    slope, intercept = np.polyfit(selected[:, 0], selected[:, 1], 1)
    points = [.005, .01, .02, .03, .05, .12]
    return dict(interval_nominal_strain=list(interval), samples=len(selected),
                least_squares_slope_MPa=float(slope), intercept_MPa=float(intercept),
                diagnostic_stress_MPa={str(p): float(np.interp(p, data[:, 0], data[:, 1]))
                                      for p in points if data[0, 0] <= p <= data[-1, 0]})


def coupon_diagnostic(pdf):
    with tempfile.TemporaryDirectory(prefix='dong-coupon-audit-') as tmp:
        prefix = Path(tmp)/'page2'
        subprocess.run(['pdfimages', '-f', '2', '-l', '2', '-j', str(pdf), str(prefix)], check=True)
        image_path = Path(str(prefix)+'-001.jpg')
        rgb = np.asarray(Image.open(image_path).convert('RGB')).astype(int)
        if rgb.shape != (718, 1422, 3):
            raise ValueError('Expected original Fig. 2 RGB raster 1422x718')
        image_hash = digest(image_path)
    r, g, b = np.moveaxis(rgb, 2, 0)
    masks = {'black': (r < 100) & (g < 100) & (b < 100),
             'colored': ((b-r > 30) & (b-g > 30)) | ((r-g > 30) & (b-g > 30))}
    result = {}
    for name, mask in masks.items():
        values = []
        for x in range(686, 729):
            ys = np.flatnonzero(mask[250:535, x])+250
            if len(ys) >= 2:
                y = np.median(ys)
                values.append([(x-678)*.05/(1393-678), (547-y)*350/(547-23)])
        a = np.asarray(values)
        if len(a) < 20:
            raise ValueError('Insufficient initial-branch coupon pixels')
        slope, intercept = np.polyfit(a[:, 0], a[:, 1], 1)
        result[name] = dict(samples=len(a), initial_figure_strain_range=a[[0, -1], 0].tolist(),
                            free_intercept_slope_MPa=float(slope), intercept_MPa=float(intercept),
                            through_origin_slope_MPa=float(np.sum(a[:, 0]*a[:, 1])/np.sum(a[:, 0]**2)))
    return dict(source='Dong et al. (2019), Fig. 2(b), PDF page 2',
                embedded_image_sha256=image_hash,
                original_image_size_pixels=[1422, 718],
                axes_pixels=dict(x0=678, x1=1393, y0=547, y1=23),
                axis_range=dict(strain=[0, .05], stress_MPa=[0, 350]),
                method='Per-column median dark/colored trace pixels x=686..728, y=250..534; unweighted stress versus strain line fit',
                estimates=result,
                uncertainty='Diagnostic raster estimate, not a measured modulus or confidence interval. At least +/-2 pixels per coordinate (strain .000140, stress 1.34 MPa), plus trace width, branch choice, possible toe/grip compliance, and unverified engineering/true definitions.',
                action='Exploratory printed-curve diagnostic only. Keep user-confirmed E=69000 MPa and the already true-stress/plastic-strain hardening input unchanged; do not infer an input error or repeat engineering-to-true conversion from this raster slope.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paper', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--native-prefix', type=Path,
                        default=ROOT/'results/opencourant/thin_free_x_ams1e6')
    parser.add_argument('--converter', type=Path,
                        default=Path('/workspace/opencourant/runtime/OpenCourant/exec/anim_to_vtk_linux64_gf'))
    args = parser.parse_args()
    config_path = ROOT/'configs/dong2019_brick_damage.json'
    config = json.loads(config_path.read_text())
    confirmation_path = ROOT/'validation/materials/user_material_confirmation.json'
    confirmation = json.loads(confirmation_path.read_text())
    fixed = confirmation['fixed_inputs']
    if (config['young'] != fixed['young_MPa'] or config['poisson'] != fixed['poisson']
            or config['density'] != fixed['density_Mg_per_mm3']):
        raise ValueError('Audit configuration differs from user-confirmed fixed material inputs')
    poly = core_polygon(config)
    xmin, ymin, xmax, ymax = poly.bounds
    h, shift, pitch = parameters(config)
    reference = ROOT/'validation/paper_reference/dong2019_fig12_thin.csv'
    with reference.open() as stream:
        experiment = [[float(r['strain']), float(r['stress_MPa'])] for r in csv.DictReader(stream)
                      if r['curve'] == 'experiment_2']
    ccx_path = ROOT/'validation/native-final.json'
    ccx = json.loads(ccx_path.read_text())
    native_assessment = args.native_prefix.with_name(args.native_prefix.name+'.assessment.json')
    native_metadata = args.native_prefix.with_name(args.native_prefix.name+'.radioss.json')
    history_path = args.native_prefix.with_name(args.native_prefix.name+'T01.csv')
    native = json.loads(native_assessment.read_text())
    metadata = json.loads(native_metadata.read_text())
    curves = dict(paper_experiment_2=curve_diagnostic(experiment))
    for name, assessment in [('historical_static_ccx', ccx), ('native_free_x', native)]:
        curves[name] = curve_diagnostic([[r['strain'], r['nominal_stress_MPa']]
                                         for r in assessment['reaction_history']])
    history = read_history_csv(history_path)
    samples = np.asarray(metadata['displacement_fraction_samples'])
    strain = np.interp(history['time'], samples[:, 0], samples[:, 1])*metadata['target_nominal_strain']
    early = (strain >= .005) & (strain <= .02)
    eroded_column = next(name for name in history if name.split() == ['Part', '1', 'ERODED'])
    native_states = []
    for index in (10, 17):
        path = args.native_prefix.with_name(args.native_prefix.name+f'A{index:03d}')
        frame = read_frame(convert(path, args.converter))
        core = frame['part_id'] == 1
        native_states.append(dict(animation=path.name, sha256=digest(path),
                                  actual_time_seconds=frame['time'],
                                  actual_nominal_strain=strain_at_time(frame['time'], metadata),
                                  max_core_element_damage=float(frame['damage'][core].max()),
                                  eroded_core_elements=int(np.sum(core & ~frame['active']))))
    geometry = dict(rows=config['rows'], columns=config['cols'], wall_mm=config['wall'],
                    depth_mm=config['depth'], actual_width_mm=xmax-xmin, actual_height_mm=ymax-ymin,
                    centerline_cell_height_mm=h, bevel_mm=config['bevel'], shift_mm=shift, pitch_mm=pitch,
                    area_mm2=poly.area, solid_fraction_actual_envelope=poly.area/((xmax-xmin)*(ymax-ymin)),
                    paper_reported_relative_density=.109, paper_nominal_width_mm=76,
                    paper_nominal_height_mm=48, nominal_stress_area_mm2=1900,
                    conclusion='Thin geometry envelope/density deviations are about percent scale; no factor-of-three dimensional or area error found.')
    report = dict(status='read-only physical discrepancy audit; no calibration or experimental agreement claim',
                  reproduction=dict(script='scripts/audit_physical_discrepancy.py',
                                    script_sha256=digest(Path(__file__)),
                                    command='python scripts/audit_physical_discrepancy.py /path/to/user_supplied_paper.pdf --output validation/opencourant_physical_discrepancy_audit.json',
                                    required_retained_native_prefix=str(args.native_prefix),
                                    converter_sha256=digest(args.converter)),
                  paper_sha256=digest(args.paper),
                  paper_page_evidence={
                      'page_3_section_2_1': 'Six transverse and five horizontal cells; nominal 48x76 mm, depth25 mm; wall .3 mm; H8 mm, L1=12 mm, L3=5.66 mm, 2L2=12 mm, theta45 degrees; relative density .109.',
                      'page_7_section_3_4': 'Crushing stress is recorded crosshead force divided by cross-section 76x25 mm2.',
                      'page_8_section_4_1': 'Elastic steel platens; C3D8R; horizontal/bevel edges19/10 elements; total10068; out-of-plane displacement constrained; platen friction .35; velocity2 mm/min.',
                      'page_2_figure_2_b_and_page_3_section_2_2': 'Three coupon tensile traces; exploratory printed slopes do not override the user-confirmed E and already-converted plastic input. Raw test histories and plotted strain definitions are not provided here.'},
                  source_hashes={str(p.relative_to(ROOT)): digest(p) for p in [config_path, reference, ccx_path, confirmation_path]},
                  native_evidence_hashes={p.name: digest(p) for p in [native_assessment, native_metadata, history_path]},
                  geometry_units=geometry, early_curve_diagnostics=curves,
                  native_early_eroded_count_max=float(history[eroded_column][early].max()),
                  native_actual_animation_checks=native_states,
                  confirmed_material_inputs=dict(confirmation_file=str(confirmation_path.relative_to(ROOT)),
                                                  confirmation_sha256=digest(confirmation_path),
                                                  fixed_inputs=fixed,
                                                  plastic_interpretation=confirmation['plastic_table']['interpretation'],
                                                  engineering_to_true_reconversion=False,
                                                  historical_source_note='The old source comment called elasticity an editable default. Explicit user confirmation on 2026-10-27 now fixes these inputs for this task; that comment is not evidence of an input error.'),
                  coupon_raster_diagnostic=coupon_diagnostic(args.paper),
                  native_loading=dict(duration_seconds=metadata['duration_seconds'],
                                      physical_duration_at_2mm_per_min_seconds=metadata['physical_duration_seconds'],
                                      peak_speed_mm_per_s=metadata['peak_imposed_speed_mm_per_s'],
                                      ams_requested_time_step_seconds=metadata['ams']['requested_time_step_seconds']),
                  interpretation=[
                      'Early stress/stiffness mismatch precedes material damage in inspected actual frames. Changing fracture-locus coordinates cannot repair the undamaged elastic slope.',
                      'Removing the optional core X pin did not remove the early native stiffness excess: inspect both independent actual run curves before assigning cause.',
                      'Beam bending controls passed in stable cases, but they do not qualify the actual skewed, one-depth-layer lattice with contact and AMS.',
                      'Matched elastic small-strain lattice controls across CCX/native, timestep/AMS, duration and contact variants should precede inverse material fitting.',
                      'E=69000 MPa, nu=.3 and density=2.67e-9 Mg/mm3 are user-confirmed fixed inputs. The supplied hardening table is already true stress versus equivalent true plastic strain and must not be re-converted.',
                      'The exploratory coupon-raster slope is not a diagnosed material input error. Do not change E or fit material interpretation to the thin-wall validation curve.',
                      'No new definite dimensional/unit conversion bug was established by this read-only audit.'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(curves, indent=2))


if __name__ == '__main__':
    main()
