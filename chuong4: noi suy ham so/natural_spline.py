import numpy as np

def cubic_spline_coefficients(x, y):
    """
    Tính hệ số spline bậc ba cho các điểm (x, y)
    Trả về các mảng b, c, d sao cho
    với mỗi đoạn [x[i], x[i+1]], spline là:
    S_i(x) = y[i] + b[i]*(x - x[i]) + c[i]*(x - x[i])^2 + d[i]*(x - x[i])^3
    """
    n = len(x)
    h = np.diff(x)

    # Tính alpha theo công thức spline
    alpha = [0] + [3 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1]) for i in range(1, n-1)]

    # Khởi tạo các mảng hỗ trợ
    l = np.ones(n)
    mu = np.zeros(n)
    z = np.zeros(n)

    # Giải hệ phương trình tam giác để tìm c
    for i in range(1, n-1):
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]

    b = np.zeros(n-1)
    c = np.zeros(n)
    d = np.zeros(n-1)

    # Tính b, c, d từ dưới lên
    for j in range(n-2, -1, -1):
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = (y[j+1] - y[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j]) / (3*h[j])

    return b, c, d

def print_cubic_spline(x, y, b, c, d):
    """
    In từng đoạn spline dưới dạng hàm bậc ba:
    S_i(x) = y[i] + b[i]*(x - x[i]) + c[i]*(x - x[i])^2 + d[i]*(x - x[i])^3
    """
    n = len(x)
    for i in range(n-1):
        print(f"Đoạn [{x[i]}, {x[i+1]}]: S_{i}(x) = {y[i]:.4f} + {b[i]:.4f}*(x - {x[i]:.4f}) + "
              f"{c[i]:.4f}*(x - {x[i]:.4f})^2 + {d[i]:.4f}*(x - {x[i]:.4f})^3")
        print()

# Hàm đánh giá giá trị của spline tại các điểm mới
# x_new là mảng các điểm cần đánh giá
def spline_evaluate(x, y, b, c, d, x_new):
    n = len(x)
    y_new = np.zeros_like(x_new)
    for idx, xi in enumerate(x_new):
        # tìm đoạn spline chứa xi
        i = np.searchsorted(x, xi) - 1
        if i < 0:
            i = 0
        elif i >= n-1:
            i = n-2
        dx = xi - x[i]
        y_new[idx] = y[i] + b[i]*dx + c[i]*dx**2 + d[i]*dx**3
    return y_new

if __name__ == "__main__":
    x = np.array([0, 2, 4, 6])
    y = np.array([1, 9, 41, 41])

    b, c, d = cubic_spline_coefficients(x, y)
    print_cubic_spline(x, y, b, c, d)

    test_points = [0.5, 1.5, 2.5, 3.5]
    y_test = spline_evaluate(x, y, b, c, d, test_points)
    print("Giá trị của spline tại các điểm mới:")
    for xi, yi in zip(test_points, y_test):
        print(f"Spline({xi}) = {yi:.4f}")