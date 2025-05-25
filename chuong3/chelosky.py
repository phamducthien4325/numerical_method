# code Chatgpt, không code lại vì chelosky chỉ tính được ma trận là là đối xứng xác định dương

import numpy as np

def cholesky_decomposition(A):
    """Thực hiện phân rã Cholesky: A = L * L.T"""
    n = A.shape[0]
    L = np.zeros_like(A, dtype=np.float64)
    
    for i in range(n):
        # Tính phần tử đường chéo
        L[i, i] = np.sqrt(A[i, i] - np.sum(L[i, :i] ** 2))
        
        for j in range(i + 1, n):
            # Tính phần tử ngoài đường chéo
            L[j, i] = (A[j, i] - np.sum(L[j, :i] * L[i, :i])) / L[i, i]
    
    return L

def solve_cholesky(A, b):
    L = cholesky_decomposition(A)
    y = np.zeros_like(b, dtype=np.float64)
    for i in range(len(b)):
        y[i] = (b[i] - np.sum(L[i, :i] * y[:i])) / L[i, i]
    x = np.zeros_like(y, dtype=np.float64)
    for i in range(len(y) - 1, -1, -1):
        x[i] = (y[i] - np.sum(L.T[i, i + 1:] * x[i + 1:])) / L[i, i]
    return x

if __name__ == '__main__':
    A = np.array([
        [4, 2, 4, 2],
        [2, 3, 4, 2],
        [4, 3, 6, 3],
        [2, 3, 3, 9]
    ], dtype=float)
    b = np.array([[20], [36], [60], [122]], dtype=float)
    x = solve_cholesky(A, b)
    print("Nghiệm của hệ phương trình:", x.reshape(-1))
    print("Kiểm tra:")
    print("Ax - b =", np.dot(A, x.reshape(-1, 1)) - b)