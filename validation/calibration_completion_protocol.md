# Quy trình calibration đề xuất từ các kiểm chứng đã chạy

Mục tiêu là một mô hình hội tụ về số, fit trên mẫu dày rồi dự đoán được mẫu
0,3 mm với cùng bộ vật liệu. Chưa có cơ sở tuyên bố calibration hoàn chỉnh.
Các ngưỡng đề xuất bên dưới là tiêu chí nghiên cứu của workflow này, không
phải các giá trị paper công bố.

## Vì sao cần đổi trình tự

- Coupon direct/RP lệch phản lực RMS 0,00375%; lấy curve qua RP không tự sửa
  được sai lệch honeycomb.
- Mesh cũ lệch stiffness đàn hồi 72,27% so với tham chiếu static cùng mesh.
  Mesh bám thành giảm chênh lệch work stiffness xuống 0,86%. Kết quả này
  chỉ kiểm chứng thành phần đàn hồi, chưa kiểm chứng nonlinear hoặc fracture.
- Giảm tốc độ bằng cách tăng thời gian 0,01 lên 0,04 s làm curve thay đổi
  4,6873%, nhưng stress gần strain 0,12 tăng 1,185 lên 1,281 MPa. Core KE/IE
  giảm xuống 2%, còn global KE/IE 18,3% nên vẫn chưa đạt quasi-static.
- Vật liệu hiện tại độc lập strain rate. Chưa có dữ liệu để fit một luật
  rate-dependent; không dùng tham số tốc độ để bù sai số mesh.

Xem [kết quả tải](loading_rate012.json) và [đồ thị](loading_rate012.png).

## 1. Khóa vật lý, chuẩn hóa hình học từng mẫu

Giữ E=69000 MPa, nu=0,3, density=2,67e-9 Mg/mm³ và bảng true stress/plastic
strain đã cung cấp. Ghi hash của bảng dùng trong mỗi run; không chuyển đổi
lần hai hoặc thay E để ép curve trùng hình.

Kiểm tra riêng mẫu wall=1 mm và 0,3 mm: kích thước ngoài, chiều sâu 25 mm,
cell arrangement, bề dày thành, khối lượng/relative density và định nghĩa
diện tích danh nghĩa. Đổi riêng `wall` hiện làm envelope thay đổi: khoảng
76,497×48 mm cho thin và 78,203×48,7 mm cho thick. Diện tích envelope×depth
là 1912,43/1955,09 mm², còn curve hiện dùng A0=1900 mm². Phải xác định
quy ước đo của specimen/paper; không tự đổi A0 để fit. Mức chênh này nhỏ
hơn nhiều sai lệch stress đang quan sát.

Giữ elastic steel plates, self-contact, U3 của toàn bộ core bị khóa và
core/platen friction=0,35 theo đoạn paper. Core self-friction chưa được
paper cho biết; giá trị 0,35 là giả thiết cần kiểm tra độ nhạy. Không thêm
pin X trong core. Không coi C3D8R là bằng chứng paper dùng Explicit.

## 2. Kiểm chứng mesh, contact và tải bằng J2 không damage

Bước giải tiếp theo có giá trị nhất là structured-wall J2/RP đến strain 0,12,
không phải sweep fracture trên mesh cũ. Hai config và bốn đầu vào đã chuẩn bị:

- [Mesh cơ sở](../configs/dong2019_structured_j2_qualification012.json).
- [Mesh mịn](../configs/dong2019_structured_j2_qualification012_refined.json).
- [Manifest và lệnh chạy](next_qualification_study.json).

| Case | Thời gian, s | AMS target, s | Refinement level | Mục đích |
| --- | ---: | ---: | ---: | --- |
| qA_T004 | 0,04 | 1e-6 | 1 | Mốc kiểm chứng mới trên mesh bám thành |
| qB_T008 | 0,08 | 1e-6 | 1 | Chỉ giảm tốc độ |
| qC_AMShalf | 0,08 | 5e-7 | 1 | Chỉ giảm AMS target |
| qD_meshrefined | 0,08 | 5e-7 | 2 | Chỉ refine mesh trong mặt phẳng |

`structured_refinement_levels=1/2` cho 2/4 phần tử qua các thành thẳng;
level 2 đồng thời refine dọc thành và junction. Chỉ giảm `mesh_size` không
đủ để kiểm tra convergence qua bề dày thành. Giữ default HEPH pressure
treatment: Icpre1 từng làm mềm sai coupon độc lập.

Các deck trên đã qua native starter, chưa chạy engine. Chúng là accelerated
smoothstep diagnostics, không phải tốc độ vật lý 2 mm/min. Nếu duration
hoặc AMS chưa hội tụ thì tiếp tục tăng thời gian hoặc giảm AMS; không mặc
định 0,08 s là đủ. J2-only được dùng để kiểm chứng số, không dùng nó để thay
thế fracture khi đánh giá fit thực nghiệm.

Tiêu chí trước khi fit: đủ observed RP motion và strain coverage; sai lệch
curve RMS chuẩn hóa không quá 5% giữa các control trên cùng khoảng strain;
stiffness đàn hồi sai khác khoảng 2% hoặc nhỏ hơn trên benchmark độc lập;
core/global KE/IE và hourglass/IE không quá 5% sau startup cho J2 control;
kiểm tra work/energy, penetration và contact work. Xem toàn bộ histories,
không cắt dao động hoặc lọc đỉnh để vượt ngưỡng. Các fracture transient về
sau cần được báo riêng và kiểm chứng, không tự động bỏ qua.

Kiểm tra contact/platen tiếp theo gồm penalty/gap sensitivity, refinement
elastic steel plates (hiện mỗi platen chỉ có một brick), rồi chạy các
control phù hợp đến strain collapse. Một chứng nhận đến 0,12 không được
dùng để chứng nhận toàn curve đến 0,78.

## 3. Thống nhất tốc độ và xác minh lực/energy

Với mỗi specimen, T_physical=60*H_actual*epsilon_target/2. Thin H=48 mm
cho T=172,8 s tại 0,12 và 1123,2 s tại 0,78. Nếu dùng acceleration factor
alpha để giải nhanh mô hình rate-independent, T=T_physical/alpha cho từng
chiều cao/hành trình. Không dùng một T chung cho thick strain 0,5 và thin
strain 0,78. Ghi rõ linear hay smoothstep; smoothstep có peak speed bằng
1,875 lần mean speed. Alpha lớn chỉ được chấp nhận sau convergence.

Với physical alpha=1, dùng linear ở 2 mm/min. Accelerated smoothstep là
giả thiết numerical khác, cần kiểm tra tốc độ/shape sensitivity ở nghiệm
đã hội tụ. Không đánh đồng một run tốc độ cao với thí nghiệm vật lý.

Study generator đã có speed mode để áp dụng quy tắc này, ví dụ lệnh chỉ
chuẩn bị hypotheses ở tốc độ vật lý, chưa chạy engine:

```sh
python calibration_study.py prepare --backend opencourant \
  --config configs/dong2019_structured_velocity012.json \
  --physical-speed-mm-min 2 --holdout-end-strain .78 \
  --output results/future_physical_hypotheses
```

Với hình học hiện tại, generator cho T_thick=730,5 s tại strain 0,5 và
T_thin=1123,2 s tại 0,78. `--acceleration-factor` có thể dùng sau khi đã
kiểm chứng acceleration thích hợp; manifest ghi mean/peak speed và duration
riêng từng role. Fixed-duration mode cũ được giữ để tái lập các control.
Geometry và native selection gates vẫn chưa được mở.

Lấy strain từ relative RP DY và stress từ signed reaction/A0, giữ cả
top/bottom RP và contact histories. Native REACY là impulse nên đạo hàm
theo thời gian đúng một lần. Giải quyết ý nghĩa master-RP inertia bằng
coupon platen không contact rồi có contact, đối chiếu impulse với Mv và
force với Ma. Chỉ hiệu chỉnh quán tính khi phép cân bằng đã được kiểm chứng;
không bỏ global energy screen hoặc lấy trị tuyệt đối để che sai dấu.

## 4. Kiểm chứng damage regularization trước inverse fit

Giữ anchor fracture strain 0,047 tại triaxiality 0,333 và failure
displacement 0,001 mm. Các triaxiality ứng với 0,08–0,018 chưa được paper
công bố. Endpoint clamping và compression damage là giả thiết, chưa phải
material đã xác định.

Native INIEVO ILEN=1 dùng characteristic length từ initial critical step;
uf bằng nhau chưa bảo đảm cùng fracture work với Abaqus hoặc CCX. Chạy
tensile/shear coupon có aspect ratio, hướng và kích thước giống wall cells,
kiểm tra initiation, softening displacement, dissipated energy/area,
deletion và contact-after-deletion khi refine mesh. Unit cube test hiện có
chưa đủ cho việc chuyển uf sang mesh honeycomb.

Ở vị trí first fracture, lưu triaxiality, PEEQ và damage để biết nhánh locus
nào được kích hoạt. Không dùng damage để sửa sai initial stiffness.

## 5. Fit thick, khóa tham số, validate thin 0,3 mm

Fit trên hai experimental curves của mẫu thick đến ít nhất strain 0,5.
Chỉ tối ưu phần locus còn thiếu; bắt đầu ít tham số, sau đó dùng completed
perturbation runs và sensitivity/SVD để xem chúng có nhận dạng được không.
Nếu nhiều bộ cùng phù hợp, báo miền tham số/bất định thay vì một nghiệm duy
nhất. Giữ friction core/platen và vật liệu nền theo dữ liệu đã biết.

Tính loss trên cùng lưới strain, dùng cả hai thí nghiệm; báo RMS/bias,
initial slope, plateau, absorption và densification riêng. Dùng replicate
spread cùng digitization resolution làm mức bất định; hai thí nghiệm không
đủ để tự suy ra confidence interval chặt. Không lấy đường FE màu xanh của
paper làm target thí nghiệm. Không extrapolate vùng dữ liệu thiếu.

Khóa toàn bộ tham số sau khi fit thick, chuyển sang thin 0,3 mm không retune.
Kiểm tra đến khoảng strain 0,78 nếu reference thực sự phủ được vùng đó;
đối chiếu curve lẫn deformation tại observed strain 0,12/0,35/0,65 và vùng
densification. Plateau error khoảng 10% và densification shift khoảng 0,03
có thể dùng làm diagnostic ban đầu, cần chốt cùng độ phân tán/digitization
trước lựa chọn cuối. Một điểm strain 0,12 không đủ chứng minh calibration.

## Điều kiện hoàn tất và phần mềm còn cần làm

Chỉ đánh dấu hoàn tất khi: hình học đúng; elastic/nonlinear/contact/damage
controls đạt trong toàn strain range cần dùng; tải/AMS/mesh đã hội tụ;
fit thick có identifiability/uncertainty hợp lý; thin holdout dự đoán đạt
curve và deformation với cùng tham số.

Hiện native scorer trong `calibration_study.py` luôn chặn selection để bảo
vệ khỏi các kiểm chứng chưa đạt. Sau khi có đủ evidence, cần bổ sung verifier
cho qualification receipts, kiểm tra hashes của geometry/mesh/contact/
material/loading, paired controls và strain range. Không mở gate bằng cờ
`qualified=true` thủ công và không bỏ blocker chỉ để workflow báo thành công.

Đích thực tế là calibration có bằng chứng và giới hạn sử dụng rõ ràng.
Độ thiếu của paper/dữ liệu hiện có không cho phép hứa tái tạo duy nhất mọi
chi tiết của mô hình ABAQUS gốc.
