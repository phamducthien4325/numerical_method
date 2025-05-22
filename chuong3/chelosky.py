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

# Ví dụ sử dụng
A = np.array([[4, 2, 4],
              [2, 10, 8],
              [4, 8, 17]], dtype=np.float64)

L = cholesky_decomposition(A)
print("Ma trận L:")
print(L)
print("Kiểm tra L * L.T:")
print(np.dot(L, L.T))