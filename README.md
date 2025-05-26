# 📂 dynamic_roboticarm

- The files `asdasd` and `tinhdaoham` directly compute the torques applied to each joint of the robotic arm.
- From these results, the dynamic matrices **M** (inertia), **C** (Coriolis/centrifugal), and **G** (gravity) can be derived.
- In cases where the **C** matrix is too complex to simplify manually, you can use the `tinhC` code (requires the **M** matrix as input).

- File `asdasd` và file `tinhdaoham` sẽ trực tiếp tính ra torque cấp vào cho mỗi khớp trong cánh tay.
- Từ kết quả đó, ta có thể rút gọn ra các ma trận **M**, **C**, **G** trong phương trình động lực học.
- Trong trường hợp ma trận **C** quá khó để rút gọn bằng tay, có thể sử dụng code `tinhC` (yêu cầu có sẵn ma trận **M**).

