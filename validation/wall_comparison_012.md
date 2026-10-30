# Đối chiếu mẫu 1 mm và 0,3 mm: kết quả đã chạy

Các curve dưới đây là J2 **không damage**, dùng cùng bảng vật liệu gốc,
RP trên mặt ngoài platen, contact `inner_y`, core-only AMS 1e-6 s và
smoothstep. Đây là numerical controls đến strain xấp xỉ 0,12, chưa phải
calibration đầy đủ cho Fig.12(a)/(b).

Xem [đồ thị phóng to](wall_comparison_012_zoom.png),
[toàn vùng paper với coverage mô phỏng được giữ nguyên](wall_comparison_012_full_paper.png),
[JSON có hashes](wall_comparison_012.json) và [CSV signed](wall_comparison_012.csv).
[Gói evidence](wall_comparison_012_inputs.tar.gz) giữ input/output/receipts của
bốn RP controls và hai thick elastic controls mới, kèm manifest SHA256 đã
kiểm tra lại từng file. Thin elastic outputs vẫn thuộc các study trước;
không tuyên bố archive này tự chứa toàn bộ lịch sử kiểm chứng.

## Plasticity/contact: các run hoàn tất

| Mẫu | Thời gian tải, s | Stress ở interval cuối gần 0,12, MPa | Max core KE/IE | Max global KE/IE |
| --- | ---: | ---: | ---: | ---: |
| 1 mm, nhanh | 0,01014583 | 19,99368 | 3,177% | 11,130% |
| 1 mm, chậm | 0,04058333 | 19,56748 | 0,192% | 0,616% |
| 0,3 mm, nhanh | 0,01 | 1,18505 | 35,217% | 222,267% |
| 0,3 mm, chậm | 0,04 | 1,28147 | 2,003% | 18,296% |

Printed energy-error screen cũng được giữ riêng: max absolute error sau
startup lần lượt là 5,4%/3,4% cho thick nhanh/chậm và 27,3%/13,3% cho thin
nhanh/chậm. Chỉ thick chậm đạt ngưỡng 5% của screen này. Các giá trị lớn
đầu quá trình cần được xem cùng absolute work/energy; không lọc bỏ chúng
để thay đổi kết luận chấp nhận. Không thấy zero/negative-volume hoặc
nonfinite warning trong bốn engine outputs này.

JSON bên ngoài archive đã được bổ sung các diagnostics này từ cùng raw
outputs. Archive là snapshot của bản report ban đầu tại commit d9b0030;
raw inputs/outputs và receipts trong archive không thay đổi.

Observed force intervals kết thúc hơi trước strain 0,12; không extrapolate
để tạo giá trị đúng 0,12. Hai đường thực nghiệm tại 0,12 là 5,625/4,340 MPa
cho mẫu 1 mm, và 0,28123/0,25912 MPa cho mẫu 0,3 mm. Đây là dữ liệu số hóa
từ hình, không phải raw experimental data.

Trong mỗi cặp, cùng geometry/material/AMS/profile và chỉ thay thời gian.
Hai mẫu có mean crosshead speed bằng nhau ở từng mức nhanh/chậm:
576/144 mm/s, khác xa tốc độ vật lý 2 mm/min. Đã dùng actual H=48,7/48 mm
để tính duration. Smoothstep peak speed bằng 1,875 lần mean.

Sai khác curve L2 trên strain 0,005–0,119 khi đổi thời gian gấp bốn là
**2,427% cho mẫu 1 mm**, so với **4,687% cho mẫu 0,3 mm**. Mẫu 1 mm chậm
đạt cả hai energy screens, nhưng chưa có đủ mesh/AMS/contact/work-balance
evidence để gọi toàn mô hình là quasi-static validated. Tải nhanh không
giải thích được toàn bộ sai lệch stress quan sát được.

Paper mẫu 1 mm có first peak khoảng 13,49/13,44 MPa tại strain 0,0574/0,0594,
sau đó drop mạnh. Curve J2 vẫn tăng đến hơn 21 MPa quanh vùng 0,08. Việc
không có fracture trong control này phải được phân biệt với lỗi numerical:
bài đối chứng cùng thiết lập có native damage đang được chạy. Chưa kết luận
damage là nguyên nhân duy nhất khi chưa có kết quả đối chứng đó.

## Đàn hồi: sai lệch không giống nhau ở hai bề dày

[Các cặp core không contact/damage](wall_elastic_comparison.json) dùng cùng
mesh, vật liệu và BC giữa native HEPH và static C3D8I; hashes và các bộ BC
được kiểm tra lại từ input/output thực tế.

| Mẫu/mesh | Native work stiffness, N/mm | Static C3D8I, N/mm | Chênh lệch |
| --- | ---: | ---: | ---: |
| 1 mm, original | 41242,718 | 40526,530 | +1,767% |
| 1 mm, structured | 37956,442 | 36636,850 | +3,602% |
| 0,3 mm, original | 1006,209 | 584,097 | +72,267% |
| 0,3 mm, structured | 590,170 | 585,135 | +0,860% |

Không có bằng chứng cho một hệ số stiffness sai chung áp dụng cả hai mẫu.
Mesh/formulation artifact của thành mỏng lớn hơn rõ rệt. Tuy vậy, static
thick stiffness đổi khoảng -9,60% giữa hai mesh family cùng thể tích, nên
thick mesh convergence vẫn chưa được xác nhận.

Không so trực tiếp stiffness tại displacement 0,01 mm của benchmark đàn hồi
với slope macro curve 0,005–0,02: vùng macro đó có thể đã plastic và có
experimental seating/toe. Không suy ra E mới từ một phép chia hai slope.

## Những lỗi chung đã kiểm tra và những giả thiết còn mở

[Audit hình học](wall_geometry_equivalence.json) xác nhận core Python bằng
core từ thuật toán offset/ring/ligament trong file Abaqus gốc ở cả hai wall:
symmetric-difference area bằng 0. Cả hai dùng actual height để tính stroke
và A0=1900 mm² để tính stress. Chưa thấy lỗi chuyển CAD hay đơn vị đủ giải
thích sai lệch lớn.

Điều này chỉ chứng minh bằng file người dùng, chưa chứng minh đúng specimen
paper. Literal centerline height 47,7 mm cho envelope 48 mm ở wall0,3 và
48,7 mm ở wall1. Diện tích envelope×depth cao hơn A0 khoảng 0,65%/2,90%,
không đủ giải thích mức stress đang sai. Cần đối chiếu các kích thước và E
với paper trước inverse fit.

Các khác biệt cross-solver còn cần control: platen native một brick,
thickness10 mm/width envelope+2 mm; file Abaqus gốc thickness25 mm/width76 mm
với mesh2 mm; khối lượng và contact coverage vì thế khác nhau. Không coi
những khác biệt này là nguyên nhân đã chứng minh khi chưa tách thí nghiệm.

Các bài bổ sung đang ưu tiên là thin structured J2 đến 0,12 và thick damage
hypothesis tại cùng điều kiện với thick J2 chậm. Workflow
`thick-wall-trials.yml` chạy thêm hai bài 1 mm đến 0,5 trên GitHub, J2 và
damage hypothesis, ở cùng accelerated mean speed576 mm/s. Các tọa độ locus
còn thiếu vẫn là giả thiết; không có tham số nào được fit trong các run này.
Kết quả chỉ được xuất bản như completed khi receipt và observed RP motion
đạt; failure/partial histories được giữ trong artifact.

Calibration vẫn chưa hoàn tất. Xem [protocol](calibration_completion_protocol.md)
để biết các kiểm chứng cần đạt trước khi chọn tham số.
