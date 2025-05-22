def newton_method(f, f_prime, p0, tol=1e-5, max_iter=100):
    i = 1
    while i < max_iter:
        p = p0 - f(p0) / f_prime(p0)
        if abs(p - p0) < tol:
            return p
        i += 1
        p0 = p
    print("Không tìm thấy nghiệm sau {} vòng lặp".format(max_iter))
    return None

# Hàm f(x) = x^2 - 2 và đạo hàm f'(x) = 2x
def f(x):
    return x**2 - 2

def f_prime(x):
    return 2 * x

# Chọn giá trị khởi tạo x0
x0 = 1.0

# Tìm nghiệm của f(x) = 0 (căn bậc hai của 2)
result = newton_method(f, f_prime, x0)
print("Nghiệm gần đúng là:", result)
