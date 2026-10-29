# Chạy mẫu thin wall t = 0,3 mm bằng Abaqus

Hai file `.inp` là input độc lập cho **Abaqus/Explicit**; không cần UMAT,
plugin hoặc CalculiX/OpenCourant để chạy. Bắt đầu với `thin_t03_plastic.inp`
để kiểm tra đáp ứng đàn dẻo. Sau đó chạy `thin_t03_damage_hypothesis.inp`
để kiểm tra thêm damage và xóa phần tử.

**Chưa chạy Abaqus hoặc Abaqus datacheck trong môi trường này vì không có
Abaqus được cấp phép.** Đã kiểm tra cấu trúc lưới, hướng phần tử, vật liệu,
các tập tải/contact và logic trích ODB bằng Python. Đây là input nghiên cứu
để kiểm tra chéo; calibration và experimental agreement chưa được xác lập.

## 1. Chạy input trực tiếp

Copy hai file `.inp` và `extract_odb.py` vào thư mục làm việc của máy có Abaqus.
Trong Abaqus Command, chạy bản đàn dẻo trước:

```sh
abaqus job=thin_plastic_check input=thin_t03_plastic.inp datacheck interactive
abaqus job=thin_plastic input=thin_t03_plastic.inp cpus=4 double=both interactive
```

Chỉ chạy analysis sau khi datacheck không báo lỗi. Dùng tên job riêng như trên
để giữ lại kết quả kiểm tra. Tiếp theo, nếu muốn kiểm tra giả thuyết damage:

```sh
abaqus job=thin_damage_check input=thin_t03_damage_hypothesis.inp datacheck interactive
abaqus job=thin_damage input=thin_t03_damage_hypothesis.inp cpus=4 double=both interactive
```

Trong Abaqus/CAE có thể dùng **Job → Create → Source: Input file**, chọn file
`.inp`, rồi Data Check/Submit. Nên chạy input gốc trực tiếp: import thành Model
rồi xuất lại có thể thay đổi hoặc bỏ các keyword contact mặt trong. Nếu phiên
bản Abaqus báo keyword không hỗ trợ, giữ file `.dat`/`.msg` để xác định lỗi.

Không có mass scaling được yêu cầu. Explicit có thể cần nhiều increment;
không tăng mass scaling để làm nhanh rồi coi kết quả là quasi-static.

## 2. Thông số đang có trong file

| Đại lượng | Giá trị |
| --- | --- |
| Hệ đơn vị | mm, s, Mg (tấn), N, MPa |
| Wall | 0,3 mm; 6 hàng × 5 cột |
| Core thực tế | 76,497146699 × 48 × 25 mm |
| Core mesh | 21.870 node; 7.392 C3D8R, giữ nguyên ID/tọa độ/connectivity của core native đã lưu |
| Platen | Thép đàn hồi, dày 10 mm; 2.560 C3D8R tổng cộng, lưới có kích thước mục tiêu 2,5 mm |
| Toàn mô hình | 25.500 node; 9.952 C3D8R |
| AlSi10Mg | E = 69.000 MPa, ν = 0,3, ρ = 2,67e-9 Mg/mm³ |
| Plastic | Đúng 178 điểm true stress–plastic strain bạn cung cấp; không quy đổi lại |
| Steel | E = 210.000 MPa, ν = 0,3, ρ = 7,85e-9 Mg/mm³ |
| Ma sát core/platen | 0,35 |
| Ma sát core/core | 0,35: giả thuyết kế thừa case native; paper chỉ nêu rõ core/platen |
| Tải | Smooth step quintic, ép 31,2 mm trong **0,01 s**, nominal strain = 0,65 |
| Chuẩn hóa stress | Diện tích danh nghĩa 76 × 25 = **1.900 mm²**, không lấy diện tích vật liệu đặc |
| Hourglass / bulk viscosity | ENHANCED / 0,06 và 1,2; không tương đương HEPH của native |

Core bị khóa **U3 = 0** ở toàn bộ node; các chuyển vị X/Y do contact quyết
định, không có pin X trong core. Mặt ngoài platen dưới khóa U1/U2/U3; mặt
ngoài platen trên khóa U1/U3 và áp U2. Các node khác của platen vẫn biến dạng
đàn hồi được.

General contact bao gồm core tự tiếp xúc và core với **hai mặt Y hướng vào
mẫu** của platen. Surface core là hợp của mặt ngoài và mặt trong để đưa các
mặt mới lộ ra do xóa phần tử vào miền contact. Các mặt bên của platen không
tham gia contact. Field output có `CFORCE`, `CSTRESS`, `CDISP`; bản damage
có thêm `SDEG`, `STATUS`.

Đây là **cùng lưới core với mô hình hiện tại**, không phải bản phục dựng lưới
gốc 10.068 phần tử của paper. Hai platen được chia lưới lại; case native cũ
dùng một brick/platen và có cả mặt bên trong miền contact. Abaqus cũng dùng
hourglass/contact và cơ chế damage riêng, không tái tạo AMS của OpenCourant.
Hai file JSON bên cạnh input ghi audit hình học và SHA256 nguồn/input.

## 3. Bản damage còn giả định gì?

`thin_t03_damage_hypothesis.inp` dùng native ductile damage, linear displacement
softening, `uf = 0,001 mm`, bật element deletion. Các hàng input theo thứ tự
**fracture strain, triaxiality, plastic strain rate**:

```text
0.080, 0.000, 0.
0.047, 0.333, 0.
0.018, 1.000, 0.
```

Paper chỉ xác định đầy đủ điểm **(η = 0,333; εf = 0,047)** và uf. Các tọa độ
η = 0 và 1 là giả thuyết, chưa phải dữ liệu calibration từ paper. Không có
compression cutoff, residual stiffness hoặc viscosity tùy biến của CalculiX
được chuyển sang Abaqus. Không dùng bản damage để tuyên bố đã tái lập paper.

## 4. Xuất đường stress–strain và so với paper

Sau khi analysis tạo ODB, chạy trong Abaqus Command:

```sh
abaqus python extract_odb.py thin_plastic.odb --output thin_plastic_curve.csv
abaqus python extract_odb.py thin_damage.odb --output thin_damage_curve.csv
```

Script lấy tổng **RF2 ở TOP_LOAD**, không cộng các node core hoặc platen
dưới. RF2 của Abaqus là lực tức thời tính bằng N, không phải impulse:

```text
compression_mm = -mean(U2 tại TOP_LOAD)
nominal_strain = compression_mm / 48
nominal_stress_MPa = -sum(RF2 tại TOP_LOAD) / 1900
```

CSV giữ cả phản lực platen dưới và các năng lượng ALLIE, ALLKE, ALLAE, ALLVD,
ALLFD, ALLWK, ETOTAL. File `.assessment.json` được xuất cùng CSV. Phản lực
mặt tải còn chứa quán tính platen trong quá trình tăng tốc; không đồng nhất
nó với lực contact riêng lên core. Không lấy trị tuyệt đối để che sai dấu.

Tại thư mục gốc repository, dùng Python thông thường đã cài Matplotlib để vẽ
so với hai thí nghiệm và FE đã số hóa từ Fig. 12(b):

```sh
python plot_paper_comparison.py \
  --assessment abaqus/thin_plastic_curve.assessment.json \
  --case-label "Abaqus plastic baseline; actuator reaction" \
  --output results/abaqus/thin_plastic_vs_paper
```

Sửa đường dẫn `--assessment` theo nơi bạn lưu kết quả. Lệnh tạo PNG/PDF/JSON.
Đường paper là số hóa từ hình in, không phải dữ liệu thí nghiệm thô. Chưa có
đường kết quả Abaqus nào được tạo sẵn hoặc suy đoán trong gói này.

## 5. Kiểm tra tốc độ trước khi calibration

**0,01 s là tải tăng tốc**, cùng thời lượng case native đã lưu; không phải
tốc độ 2 mm/phút trong phần FE của paper. Với 48 mm và strain 0,65, tốc độ
2 mm/phút tương ứng **936 s với displacement tuyến tính**. Phần thử nghiệm
của paper lại ghi 1 mm/phút; hai tốc độ này cần được phân biệt.

Sau khi bản đầu chạy được, nên so với thời lượng gấp đôi (0,02 s), kiểm tra
năng lượng, cân bằng lực và độ nhạy lưới; riêng ALLKE/ALLIE nhỏ chưa đủ để
kết luận hội tụ. Không điều chỉnh E hoặc bảng Plastic để bù lỗi contact/lưới.

Có thể tái tạo hai input với thời lượng khác từ repository bằng Python thường
(không cần Abaqus để export; cần dependency NumPy của repository):

```sh
python export_abaqus.py --output-dir results/abaqus_duration_double \
  --duration-seconds 0.02

# Tốc độ FE 2 mm/phút: có thể rất tốn thời gian với Explicit, không phải run nhanh.
python export_abaqus.py --output-dir results/abaqus_physical_speed \
  --duration-seconds 936 --loading linear
```

Exporter đọc core từ archive đã commit, không chạy lại Gmsh. Dùng thư mục
output mới cho mỗi biến thể để giữ nguyên dữ liệu đã chạy.

## Script CAE từ file gốc, strain 0,12

Bộ riêng [run_strain012.py](run_strain012.py) được tạo từ
`reference_setup_numeric_source.py`, giữ coupling mặt ngoài platen → RP và
xuất RF2/U2. Xem [hướng dẫn chạy](README_VI.md) và
[exporter ODB](export_abaqus_curve.py). Bộ này giữ các lựa chọn vật liệu,
damage, mesh và mass scaling của script gốc; không thay thế hai deck thin_t03
ở trên và không phải một mô hình đồng nhất với native OpenCourant.
