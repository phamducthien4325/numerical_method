def fixed_point_iteration(g, p0, tol=1e-5, max_iter=100):
    i = 1
    while i <= max_iter:
        p = g(p0)
        if p - p0 < tol:
            return p
        else:
            i += 1
            p0 = p
    print("Không tìm thấy nghiệm sau {} vòng lặp".format(max_iter))
    return None

# Hàm g(x) = (x + 2) / 2
def g(x):
    return (x + 2) / 2

# Chọn giá trị khởi tạo x0
x0 = 1.0

# Tìm nghiệm của f(x) = 0 (căn bậc hai của 2)
result = fixed_point_iteration(g, x0)
print("Nghiệm gần đúng là:", result)
