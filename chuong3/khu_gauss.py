# 10,000
import numpy as np

def gauss(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = len(A)
    Ab = np.hstack([A, b])

    for i in range(n):
        if Ab[i, i] == 0:
            for j in range(i + 1, n):
                if Ab[j, i] != 0:
                    Ab[[i, j]] = Ab[[j, i]]
                    break

        Ab[i] = Ab[i] / Ab[i, i]
        for j in range(i + 1, n):
            Ab[j] = Ab[j] - Ab[i] * Ab[j, i]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (Ab[i, -1] - np.sum(Ab[i, i + 1:n] * x[i + 1:n]))
    return x

if __name__ == '__main__':
    # n = int(input("Bạn muốn giải hệ phương trình bao nhiêu ẩn? "))
    # A = np.random.rand(n, n)
    # b = np.random.rand(n).reshape(-1, 1)
    A = np.array([[2, 1], [1, -3]], dtype=float)
    b = np.array([[-1], [5]], dtype=float)
    print("Ma trận mở rộng Ab:")
    print(np.hstack([A, b]))
    x = gauss(A, b)
    print("Nghiệm của hệ phương trình:", x)
    print("Kiểm tra:")
    print("Ax - b =", np.dot(A, x.reshape(-1, 1)) - b)