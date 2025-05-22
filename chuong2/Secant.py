def secant_method(f, p0, p1, tol=1e-5, max_iter=100):
    i = 2
    q0 = f(p0)
    q1 = f(p1)
    while i <= max_iter:
        p = p1 - q1 * (p1 - p0) / (q1 - q0)
        if abs(p - p1) < tol:
            return p
        i += 1
        p0 = p1
        q0 = q1
        p1 = p
        q1 = f(p)
    print("Khong tim thay nghiem sau {} vong lap".format(max_iter))
    return None

# Hàm f(x) = x^2 - 2
def f(x):
    return x**2 - 2

# Chọn hai giá trị khởi tạo x0 và x1
x0 = 1.0
x1 = 2.0

# Tìm nghiệm của f(x) = 0 (căn bậc hai của 2)
result = secant_method(f, x0, x1)
print("Nghiệm gần đúng là:", result)
