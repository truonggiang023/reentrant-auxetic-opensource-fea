# Điều kiện tải 2 mm/phút từ đoạn paper

`run_paper_velocity012.py` được tạo riêng từ script cũ của bạn. Paper không
nêu Standard hay Explicit trong đoạn đã cung cấp. Script này chọn Explicit
như một lựa chọn triển khai; không khẳng định paper đã dùng Explicit.

- Thành 0,3 mm, strain cuối 0,12, coupling mặt ngoài platen → RP.
- Thời lượng bằng stroke/(2/60). Với H0=48 mm: stroke 5,76 mm, T=172,8 s.
- Tabular amplitude tuyến tính, smooth=0: vận tốc không đổi 2 mm/phút.
  SmoothStep của script cũ có vận tốc thay đổi và đạt đỉnh 1,875 lần trung bình.
- Khóa U3 trên mọi node core đã mesh, không chỉ hai mặt extrusion.
- Giữ contact general/self-contact và ma sát 0,35 theo script cũ. Việc dùng
  cùng ma sát 0,35 cho core self-contact là giả định: đoạn paper chỉ quy định
  honeycomb–indenter/pedestal.
- Không yêu cầu mass scaling. Phân tích Explicit ở thời gian vật lý có thể
  rất tốn chi phí. Script tạo CAE/INP và không tự submit job.
- Giữ vật liệu và damage một điểm của file cũ. Các vị trí triaxiality còn
  thiếu trong paper vẫn chưa xác định; đây chưa phải vật liệu đã calibration.

Chạy ở thư mục mới trong Abaqus Command:

```bat
abaqus cae noGUI=run_paper_velocity012.py
abaqus job=PaperSpeed_check input=Compression_PaperSpeed012.inp datacheck interactive
```

Chỉ sau khi datacheck thành công và đã chọn tài nguyên cho bài Explicit dài:

```bat
abaqus job=Compression_PaperSpeed012 input=Compression_PaperSpeed012.inp interactive
abaqus python export_abaqus_curve.py Compression_PaperSpeed012.odb --metadata paper_speed_curve_metadata.json
```

Đã kiểm tra cú pháp, giữ nguyên bảng Plastic và các lựa chọn tải bằng Python.
Chưa chạy Abaqus/datacheck vì cloud không có Abaqus/license. Platen/lưới của
script cũ và mô hình native vẫn có khác biệt ngoài các điều kiện tải nêu trên;
không coi hai deck là cặp mô hình đồng nhất đã được kiểm chứng.
