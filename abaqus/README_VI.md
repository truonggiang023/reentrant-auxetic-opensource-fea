# Chạy Abaqus tại strain 0,12

File run_strain012.py là bản độc lập từ reference_setup_numeric_source.py của bạn.
Giữ nguyên hình học, bảng plasticity, damage, contact, mass scaling và platen của
file cũ. Thay đổi: strain cuối 0,12; tên model/job riêng; lưu CAE và INP;
history RF2/U2 mỗi 0,002 s. Coupling KINEMATIC từ mặt ngoài mỗi platen về RP.

1. Giải nén vào một thư mục làm việc mới.
2. Trong Abaqus Command, chạy:

   abaqus cae noGUI=run_strain012.py

   Hoặc Abaqus/CAE → File → Run Script → chọn run_strain012.py.
   Script tạo ReEntrant_RP_Strain012.cae, Compression_RP_Strain012.inp,
   rp_curve_metadata.json. Script không tự submit.
3. Chạy job trong Job Manager hoặc:

   abaqus job=Compression_RP_Strain012 input=Compression_RP_Strain012.inp interactive

4. Khi job hoàn thành:

   abaqus python export_abaqus_curve.py Compression_RP_Strain012.odb --metadata rp_curve_metadata.json

CSV dùng strain=-(U2_top-U2_bottom)/H0 và stress=-RF2_top/A0, nén dương,
đơn vị mm/N/MPa. A0=1900 mm², H0 lấy từ model. Đường này là đáp ứng kết cấu,
không phải bảng true stress/plastic strain của vật liệu.

Chưa chạy kiểm chứng trong Abaqus vì môi trường cloud không có Abaqus/license.
Đây là Explicit có mass scaling theo file gốc: kiểm tra ALLKE/ALLIE và cân bằng
RF2 trước khi kết luận quasi-static. File cũ có damage bật và platen dày 25 mm;
lần chạy CalculiX J2 trước dùng platen 10 mm, nên chưa phải cặp mô hình đồng nhất.
