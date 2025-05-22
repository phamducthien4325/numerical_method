#  20,000

import numpy as np

def gauss_seidel(A: np.ndarray,
                 b: np.ndarray,
                 x0: np.ndarray = None,
                 tol: float = 1e-3,
                 max_iter: int = 1000000) -> np.ndarray:
    n = len(A)

    if x0 is None:
        x0 = np.zeros(n, dtype=float64)

    for m in range(max_iter):
        x = x0.copy()
        for j in range(n):
            x[j] = (b[j] - np.sum(A[j, :j] * x[:j]) - np.sum(A[j, j + 1:] * x0[j + 1:])) / A[j, j]
        if np.max(np.abs(x - x0)) < tol:
            print(f"Solution obtained after {m} iteration steps.")
            return x
        x0 = x.copy()
    print(f"No solution satisfying the tolerance condition obtained after {max_iter} iteration steps.")
    return None

def generate_diagonally_dominant_system(n, min_val=-10, max_val=10):
    A = np.random.randint(min_val, max_val, (n, n)).astype(np.float64)
    for i in range(n):
        A[i, i] = sum(abs(A[i])) + np.random.randint(1, 5)  
    b = np.random.randint(min_val, max_val, (n, 1)).astype(np.float64)
    return A, b

if __name__ == '__main__':
    n = int(input("Bạn muốn giải hệ phương trình bao nhiêu ẩn? "))
    # A = np.random.rand(n, n)
    # b = np.random.rand(n).reshape(-1, 1)
    A, b = generate_diagonally_dominant_system(n)
    print("Ma trận mở rộng Ab:")
    print(np.hstack([A, b]))
    x = gauss_seidel(A, b, np.zeros(n, np.float64))  
    print("Nghiệm của hệ phương trình:", x)
    print("Kiểm tra:")
    print("Ax - b =", np.dot(A, x.reshape(-1, 1)) - b)