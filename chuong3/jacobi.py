import numpy as np

def jacobi(A: np.ndarray,
           b: np.ndarray,
           x0: np.ndarray = None,
           tol: float = 1e-3,
           max_iter: int = 1000000) -> np.ndarray:
    n = len(A)
    b = b.flatten()
    if x0 is None:
        x0 = np.zeros(n, dtype=float64)
    x = x0.copy()
    for m in range(max_iter):
        for j in range(n):
            x[j] = (b[j] - np.sum(A[j, :j] * x0[:j]) - np.sum(A[j, j + 1:] * x0[j + 1:])) / A[j, j]
        if np.max(np.abs(x - x0)) < tol:
            print(f"Solution obtained after {m} iteration steps.")
            return x
        x0 = x.copy()
    print(f"No solution satisfying the tolerance condition obtained after {max_iter} iteration steps.")
    return x



def generate_diagonally_dominant_system(n, min_val=-10, max_val=10):
    A = np.random.randint(min_val, max_val, (n, n)).astype(float)
    b = np.random.randint(min_val, max_val, (n, 1)).astype(float)
    
    for i in range(n):
        A[i, i] = abs(A[i]).sum() + np.random.randint(1, 5)  # Đảm bảo chéo trội
    
    return A, b


if __name__ == '__main__':
    # n = int(input("Bạn muốn giải hệ phương trình bao nhiêu ẩn? "))
    # # A = np.random.rand(n, n)
    # # b = np.random.rand(n).reshape(-1, 1)
    # A, b = generate_diagonally_dominant_system(n)
    A = np.array([
        [3, -1, 1],
        [3, 6, 2],
        [3, 3, 7]
    ], dtype=float)
    b = np.array([[1], [0], [4]], dtype=float)

    print("Ma trận mở rộng Ab:")
    print(np.hstack([A, b]))
    x = jacobi(A, b, np.zeros(len(A), np.float64), max_iter=2)  
    print("Nghiệm của hệ phương trình:", x)
    # print("Kiểm tra:")
    # print("Ax - b =", np.dot(A, x.reshape(-1, 1)) - b)