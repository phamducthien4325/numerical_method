def regula_falsi(f, a, b, tol=1e-5, max_iter=100):
    if f(a) * f(b) > 0:
        print("Hàm không có dấu trái ngược tại a và b!")
        return None
    
    for n in range(max_iter):
        # Áp dụng công thức điểm sai
        c = (a * f(b) - b * f(a)) / (f(b) - f(a))
        
        # Kiểm tra điều kiện dừng
        if abs(f(c)) < tol:
            return c
        
        # Cập nhật giá trị a và b
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    
    return c

# Hàm f(x) = x^2 - 2
def f(x):
    return x**2 - 2

# Chọn giá trị khởi tạo a và b
a = 1.0
b = 2.0

# Tìm nghiệm của f(x) = 0 (căn bậc hai của 2)
result = regula_falsi(f, a, b)
print("Nghiệm gần đúng là:", result)
