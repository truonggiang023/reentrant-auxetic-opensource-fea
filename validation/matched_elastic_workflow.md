# Đối chứng đàn hồi trước calibration

Mục đích là tách ảnh hưởng phần tử/lưới khỏi vật liệu phá hủy và contact.
Giữ nguyên 21.870 nodes và 7.392 bricks của baseline thin wall. Tất cả U3=0,
đáy Uy=0, đỉnh Uy=-0.01 mm, một node đáy khóa Ux để bỏ rigid translation;
các Ux còn lại tự do. E=69000 MPa, nu=.3. Không có platen, contact, plasticity
hoặc damage. Đây không phải mô hình nén hoàn chỉnh của paper.

Từ checkout, dùng input archive đã commit; không phụ thuộc ignored results
của máy tác giả:

```sh
export PATH=/workspace/bin:/workspace/fea-venv/bin:$PATH
source /workspace/opencourant/env.sh
ELASTIC_DIR="$PWD/results/calibration/elastic-check-v1"
mkdir -p "$ELASTIC_DIR/source"
tar -xzf validation/opencourant-thin-free-x-inputs.tar.gz \
  -C "$ELASTIC_DIR/source"

python scripts/probe_matched_elastic_lattice.py \
  --source "$ELASTIC_DIR/source/thin_free_x_ams1e6_0000.rad" \
  --output "$ELASTIC_DIR/ccx" --ccx /workspace/bin/ccx --threads 2

python scripts/probe_courant_elastic_lattice.py \
  --mesh "$ELASTIC_DIR/ccx/mesh_and_bcs.json" \
  --output "$ELASTIC_DIR/native24" --solid 24 --duration .004 \
  --ams-time-step 1e-6 --threads 2 --timeout-seconds 1800

python scripts/probe_courant_elastic_lattice.py \
  --mesh "$ELASTIC_DIR/ccx/mesh_and_bcs.json" \
  --output "$ELASTIC_DIR/native18" --solid 18 --duration .004 \
  --ams-time-step 1e-6 --threads 2 --timeout-seconds 1800
```

Chạy tuần tự. Mỗi native control hiện mất khoảng 13 phút trên 2 threads;
dùng output directory mới để giữ bằng chứng cũ. Nếu chỉ cần xuất shared mesh,
thêm `--prepare-only` vào lệnh CCX. Các script không sửa cấu hình vật liệu hoặc
exporter chính.

Các đối chứng static CCX đã chạy hoàn tất với zero errors/warnings:

| Phần tử | k từ phản lực, N/mm | Effective modulus, MPa |
|---|---:|---:|
| C3D8R | 448.4803 | 11.3300 |
| C3D8I | 584.0973 | 14.7561 |
| C3D8 | 2067.6930 | 52.2365 |

[Báo cáo CCX](matched_elastic_ccx_audit.json) giữ riêng phản lực và internal
energy: ELSE của C3D8R thấp hơn F·delta/2 khoảng 10.30%; không mặc định rằng
stabilization energy đã được tính vào ELSE.

Cả hai native controls cũng hoàn tất với execution-time metadata/input/output
binding. Tuy nhiên, cộng rồi vi phân hàng trăm nodal impulses cho lực bị nhiễu
và không vượt qua kiểm tra cân bằng ở cuối ramp. Không dùng trung bình lực
này để tuyên bố đã khớp static stiffness. Script báo riêng work-equivalent
stiffness `2*external_work/displacement^2` và regression của cumulative impulse
với integral(displacement), velocity và offset. Các đại lượng này là diagnostics,
vẫn cần duration/AMS refinement; global/core KE nhỏ không tự chứng minh hội tụ.

Xem [báo cáo native](matched_elastic_native_audit.json) và input archive được
ghi hash trong báo cáo. Archive có generation script gốc của mỗi run; metadata
và solver output đã quan sát không được sửa lại để gán provenance mới.

Kết quả xác nhận sự nhạy với discretization/formulation trên cùng geometry,
không xác định được fracture locus và không chứng minh experimental agreement.
Giải quyết sensitivity này trước khi chọn material candidate từ một curve fit.
