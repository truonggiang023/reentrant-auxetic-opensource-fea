# Re-entrant honeycomb FEM with CalculiX and OpenCourant

Research implementation for the Dong et al. (2019) compression problem.
**Experimental agreement and the complete Bi-V → Z → collapse sequence are not
established.** CalculiX and OpenCourant perform the FEM solves; Python handles geometry, meshing,
automation, tests and postprocessing.

GitHub Actions runs [OpenCourant native checks](.github/workflows/opencourant-smoke.yml)
on pushes to `main` and `feature/opencourant-erosion-calibration`, on pull requests,
and through `workflow_dispatch`. It verifies the pinned official runtime, runs
Python/native component regressions, and compares direct-face and surface-to-RP
coupon reactions through strain 0.12 (relative force RMS tolerance 0.1%). Download
the run artifact for solver logs, input/output hashes, CSV curves, the comparison
plot and JUnit results. The fast coupon does not qualify quasi-static behavior or
complete the paper calibration. CalculiX integrations skip without their binaries;
Abaqus scripts receive contract checks, with licensed execution still required.

The primary workflow now includes a compiled finite-strain J2 ductile-softening
material for CalculiX 2.20, separate steel platens with sliding contact, and
optional disjoint wall-contact patches. See [model and evidence](DONG2019_CALCULIX.md)
and the [native material contract](native/README.md) for qualifications and limits.

## Linux cloud setup

The prepared machine uses `/workspace/fea-venv`, `/workspace/bin/ccx` (packaged
CalculiX 2.20) and `/workspace/bin/ccx-damage` (the qualified external-material
build). No service needs to remain running.

```sh
python -m venv /workspace/fea-venv
/workspace/fea-venv/bin/pip install -e '.[dev,post]'
bash scripts/build_ccx_damage.sh
export PATH=/workspace/bin:/workspace/fea-venv/bin:$PATH
OMP_NUM_THREADS=2 CCX_EXECUTABLE=/workspace/bin/ccx \
  CCX_DAMAGE_EXECUTABLE=/workspace/bin/ccx-damage pytest -q -rx
```

The build uses signed Debian package/source metadata and checks the vendor files
before applying scoped hourglass and output-time retry compatibility patches. It installs locally in
`/workspace`, preserves TLS/package verification, and records artifact hashes.
It currently targets Debian amd64 with GCC 14 and the downloaded Fortran 12
compiler. The bundled Windows 2.13 executable is not qualified for this material.

## Generate and run

Để kiểm tra chéo trước khi tiếp tục calibration, có sẵn hai **input
Abaqus/Explicit cho thin wall t = 0,3 mm**: [bản đàn dẻo](abaqus/thin_t03_plastic.inp)
và [bản damage giả thuyết](abaqus/thin_t03_damage_hypothesis.inp).
[Hướng dẫn chạy và trích ODB](abaqus/README.md) có lệnh datacheck, chạy job,
xuất CSV và vẽ so với paper. Chúng dùng đủ 178 điểm Plastic bạn cung cấp,
giữ nguyên core mesh hiện tại và chưa được chạy/datacheck bằng Abaqus.

An additional **experimental OpenCourant backend** now uses native solid
erosion and TYPE25 all-face contact. Its installation, actual component tests,
commands and unresolved lattice qualification are documented in
[the OpenCourant workflow](validation/opencourant_workflow.md). Native deletion
has been demonstrated; complete fracture-contact coverage, calibration and
experimental agreement have not been established.

The [paper-reference data](validation/paper_reference/README.md) and
[calibration workflow](validation/calibration_workflow.md) keep digitized
experiments, published FE and current simulations distinct. Generate the
historical thin-wall comparison with `python plot_paper_comparison.py`.

The completed OpenCourant free-X baseline now has a corrected
[total-force comparison](validation/opencourant-thin-total-force-comparison.png).
TYPE25 secondary force uses FN−FT; the former normal-only plots are historical.
At strain .12 the corrected stress is 1.390 MPa versus approximately .281/.259
MPa in the experiments. Halving AMS changes the total-force curve by 14.97%
relative L2. See the [assessment](validation/opencourant-thin-total-force-qualification.json).

A matched inner-Y-only platen-contact control reached force strain .12044;
at .12 it gives 1.252 MPa. This is an intentionally partial diagnostic, not a
calibrated trajectory. The exporter now supports `platen_contact_surface=inner_y`
via [this configuration](configs/dong2019_inner_y_contact.json). The legacy
`all` setting remains the default so old studies remain reproducible.
Core all-face self-contact is retained. Calibration and experimental agreement
remain unresolved; the reported good early Abaqus response is useful external
context, not an executed OpenCourant qualification.

```sh
# Builtin J2 plasticity, wedge mesh, platen contact:
python run_paper.py --config configs/dong2019_contact.json --run

# Native coupled softening, all-quad C3D8R mesh, wall-contact candidate:
python run_paper.py --config configs/dong2019_brick_damage.json \
  --self-contact --target-output-only --run

# Plot accepted solver frames; missing targets are shown explicitly:
MPLCONFIGDIR=/workspace/mplconfig python plot_run.py results/dong2019_contact.inp
```

The damage configuration is **explicitly uncalibrated**: the paper publishes one
fracture-locus point and a range for other fracture strains, without the missing
triaxiality coordinates. Its assumptions and viscosity are written in the JSON.
Use `--imperfection-mm .003` only as a labelled imperfection sensitivity study;
that value is not a measured specimen defect.

The [user-supplied 178-point tensile flow table](validation/materials/README.md)
is preserved verbatim and audited against every value in the existing source.
Both generators accept it directly with
`--plastic validation/materials/user_plastic_table_2026-10-26.inp`.
It agrees with six-significant-digit formatting of the existing table; it adds
no fracture-locus measurements. Previously executed inputs remain unchanged.

The writer produces `.inp` and `.material.dat` together. The runner sets
`REENTRANT_MATERIAL_DATA` for that solver process and chooses `ccx-damage` when
a damage law is present. Do not replace the nine descriptor properties with a
large inline `*USER MATERIAL` table: CalculiX 2.20 has 21-entry local buffers.

The runner saves logs, preserves old result files, and verifies a completed status
file and material-load signature. A zero process exit alone is insufficient.
`frd.committed_fields` filters out unconverged “best solution” frames after a
failure. Target-output mode includes exact strains 0.12, 0.35 and 0.65; the example
also requests a regular output grid for reaction histories. Failed or unrun
frames must not be presented as deformation evidence.

## Tự chạy calibration với OpenCourant

Workflow này tạo và chạy một study tham số, tính sai lệch với paper và xuất
plot. **Hiện chưa có bộ tham số đã được xác nhận experimental agreement.**
Thick wall `t=1 mm` dùng để fitting; thin wall `t=0.3 mm` là holdout để kiểm tra
khả năng dự đoán. Không fitting cả hai rồi gọi thin wall là validation độc lập.

Các lệnh dưới chạy từ checkout trong cloud environment đã setup. Trên máy
khác, thay đường dẫn Python/runtime tương ứng. Nếu thiếu native runtime, cài
bằng `bash scripts/setup_opencourant.sh`; script kiểm tra checksum bản chính thức.

```sh
cd /workspace/reentrant-auxetic-opensource-fea
export PATH=/workspace/bin:/workspace/fea-venv/bin:$PATH
source /workspace/opencourant/env.sh
export MPLCONFIGDIR=/workspace/mplconfig
python -m pytest -q -rx
```

Giữ `E=69000 MPa`, `nu=0.3`, `density=2.67e-9 Mg/mm3`. Bảng được cung cấp đã
là **true stress–plastic strain**, không chuyển đổi lần nữa. Tham số phá hủy
ngoài điểm `(eta=0.333, strain=0.047)` vẫn là giả thuyết cần hiệu chỉnh.

**1. Tạo study mới.** Đổi tên thư mục khi cần chạy lại; không ghi đè kết quả cũ.

```sh
STUDY_DIR="$PWD/results/calibration/native-v1"
python calibration_study.py prepare \
  --backend opencourant --duration .01 --ams-time-step 1e-6 \
  --plastic validation/materials/user_plastic_table_2026-10-26.inp \
  --output "$STUDY_DIR"
```

Lệnh này chỉ tạo `manifest.json` và 22 configurations: 9 giả thuyết fracture
locus, 2 biến thể kiểm tra duration/AMS, mỗi loại có training và holdout.
`.01 s` và AMS `1e-6 s` là giá trị khởi đầu để khảo sát, chưa phải lựa chọn đã
hội tụ; tốc độ thí nghiệm là `2 mm/min`. Hai lần chạy thin wall hiện có mất
khoảng 15–30 phút/lần trên 4 threads; cả study có thể mất nhiều giờ.

**2. Xem lệnh, rồi chạy một training baseline.** `locus_0_1` là ví dụ, không
phải bộ tham số tốt nhất đã được xác nhận.

```sh
python run_calibration.py "$STUDY_DIR/manifest.json" \
  --candidate locus_0_1 --role training --dry-run

python run_calibration.py "$STUDY_DIR/manifest.json" \
  --candidate locus_0_1 --role training --threads 4 --timeout-seconds 3600
```

Runner chạy solver, chuyển T01 sang CSV và assess tự động, lưu log trong
`execution-logs/` và trạng thái trong `execution-journal.json`. Một process
exit code 0 không đủ: runner kiểm tra output thực và execution hashes.
Ctrl-C dừng cả CLI và solver con, giữ dữ liệu dở dang, ghi
`interrupted_preserved` vào journal và trả exit code 130.

**3. Chạy các đối chứng số trước khi fitting vật liệu.**

```sh
python run_calibration.py "$STUDY_DIR/manifest.json" \
  --candidate duration_double --candidate ams_half --role training \
  --threads 4 --timeout-seconds 3600
```

So sánh toàn bộ đường lực, năng lượng của core, erosion và dạng biến dạng;
không chỉ so plateau trung bình. Báo cáo hiện tại cho thấy giảm AMS một nửa
làm thay đổi đường stress khoảng 16.8% theo relative L2, nên chưa đạt hội tụ.
Nếu các đối chứng còn nhạy, cần giải quyết chúng trước khi kết luận calibration.

**4. Chạy các training candidates còn lại và tính sai lệch.**

```sh
python run_calibration.py "$STUDY_DIR/manifest.json" \
  --role training --threads 4 --timeout-seconds 3600 --resume --keep-going

python calibration_study.py score "$STUDY_DIR/manifest.json" \
  --output "$STUDY_DIR/scores.json"
```

`--resume` chỉ bỏ qua run hoàn tất có bằng chứng hợp lệ. Run timeout/failed vẫn
được giữ nguyên và không được ghi đè; tạo study mới nếu muốn thử lại. Lệnh bỏ
`--candidate` sẽ chọn tất cả candidates thuộc role đã chỉ định.
`scores.json` chứa RMSE, bias, coverage, `evaluated`, `unrun` và
`diagnostic_training_ranking`. Ranking chỉ so các training curves hoàn tất và
đủ coverage, loại numerical sensitivities; nó không chứng minh model đã đúng.
Hiện `selection.selected_candidate` vẫn bị chặn bởi geometry/numerical
qualification chưa đạt. Không đổi cờ `qualified` để bỏ qua bằng chứng này.

**5. Chạy thin-wall holdout và vẽ plot.** Giữ nguyên tham số của candidate
đã khảo sát ở training. Ví dụ dưới dùng `locus_0_1`; nếu chọn candidate khác,
dùng đường dẫn assessment/metadata ghi trong `manifest.json`.

```sh
python run_calibration.py "$STUDY_DIR/manifest.json" \
  --candidate locus_0_1 --role holdout --threads 4 --timeout-seconds 3600

python plot_paper_comparison.py \
  --assessment "$STUDY_DIR/runs/locus_0_1_holdout.assessment.json" \
  --output "$STUDY_DIR/plots/locus_0_1_thin" \
  --case-label 'Calibration candidate; experimental agreement unverified'

python plot_courant_frames.py "$STUDY_DIR/runs/locus_0_1_holdout.radioss.json" \
  --output "$STUDY_DIR/plots/locus_0_1_frames" \
  --case-label 'Calibration candidate; experimental agreement unverified'

python calibration_study.py score "$STUDY_DIR/manifest.json" \
  --output "$STUDY_DIR/scores-with-holdout.json"
```

Plot tạo PNG/PDF/JSON; thiếu dữ liệu strain sẽ được ghi rõ, không ngoại suy.
Candidate names có dấu chấm được mã hóa thành `p` trong tên file native;
hãy lấy path từ manifest thay vì đoán. `thin_wall_holdout` chỉ được gắn vào
selection khi selection hợp lệ; kết quả holdout chưa qualified vẫn nằm trong
`evaluated`. Xem thêm [quy trình và điều kiện chấp nhận](validation/calibration_workflow.md),
[native runtime/contact](validation/opencourant_workflow.md) và
[dữ liệu paper](validation/paper_reference/README.md).

Các bước prepare → solve → resume → score → plot đã được chạy thử bằng solver
thật trên một mô hình một-cell riêng; xem
[bằng chứng workflow](validation/calibration_runner_smoke.json). Đây là kiểm
tra công cụ, không phải calibration của lattice 6 hàng × 5 cột. Cả batch 22 runs
chưa được chạy đầy đủ. Nếu cần kiểm tra sai lệch ban đầu trước fitting, dùng
[đối chứng đàn hồi cùng lưới](validation/matched_elastic_workflow.md).

## CalculiX model scope

Units are mm, N and MPa. Compression is along Y and extrusion along Z. The
paper-oriented deck fixes all core U3 DOFs and removes one in-plane rigid
translation with a centreline roller. Only the outer steel platen faces receive
the prescribed/support displacement; their contact faces remain elastically
free. Static time is pseudo-time, not experimental seconds.

The C3D8R mesher uses Gmsh quad recombination and conforming refinement. Wedges
and fully integrated C3D8 bricks are also available. Geometry and mesh convergence
remain to be demonstrated; an element-type match is not a reproduction.

Damage is coupled to stress and state inside the compiled material. It uses
characteristic-length softening with retained residual stiffness. **It does not
automatically delete elements or create new crack-contact surfaces.** Retained
failed faces can affect collapse/contact. The wall algorithm uses initially
facing disjoint patches and cannot guarantee all new contacts after severe
rotations. The test suite exposes known frictional mortar convergence failures
as XFAILs; penalty contact is the selected formulation.

The legacy `reentrant-fea`/`--advanced` writers and standalone Python material-point
prototypes remain for comparison. The compiled material workflow is
`run_paper.py` with an explicit `ductile_damage` configuration, rather than the
legacy advanced writer's unsupported damage switch.

## Paper-following calibration revision

See [the revised evidence report](validation/paper_calibration_revision.md) and
[full Fig.12(b) comparison](validation/paper_calibration_revision.png).
The matched no-AMS elastic control still exceeds the static C3D8I reference
stiffness by 72.27%; inverse fracture fitting remains blocked. New study
preparation accepts `--holdout-end-strain .78` to assess the late visible
thin-wall response, preserving thick-wall training and independent holdout roles.

The supplied paper paragraph does not report the solver/step formulation.
[The loading correction](validation/calibration_loading_revision.md) records
constant 2 mm/min, linear displacement, self-contact, core U3 and platen friction.
The optional wall-aligned mesh reduces the measured **elastic component**
work-stiffness discrepancy to 0.86%; full collapse calibration is still separate.
Generate a physical-speed input with `run_radioss.py --physical-speed-mm-min 2`;
this is distinct from an accelerated explicit duration. The Abaqus physical-speed
variant and its execution limitations are described in
[the Abaqus instructions](abaqus/README_PAPER_VELOCITY.md).
