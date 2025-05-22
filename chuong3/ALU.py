import numpy as np
# Con goi la phuong phap doolittle
def lu_decomposition(A):
    n = len(A)
    L = np.eye(n) 
    U = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            U[i, j] = A[i, j] - sum(L[i, k] * U[k, j] for k in range(i))
        for j in range(i+1, n):
            L[j, i] = (A[j, i] - sum(L[j, k] * U[k, i] for k in range(i))) / U[i, i]
    
    return L, U

def forward_substitution(L, b):
    n = len(L)
    y = np.zeros(n)
    
    for i in range(n):
        y[i] = (b[i] - sum(L[i, k] * y[k] for k in range(i))) / L[i, i]
    
    return y

def backward_substitution(U, y):
    n = len(U)
    x = np.zeros(n)
    
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - sum(U[i, k] * x[k] for k in range(i+1, n))) / U[i, i]
    
    return x

def solve_lu(A, b):
    """ Giải hệ phương trình Ax = b bằng LU """
    L, U = lu_decomposition(A)
    y = forward_substitution(L, b)
    x = backward_substitution(U, y)
    return L, U, x

# if __name__ == '__main__':
#     n = int(input("Bạn muốn giải hệ phương trình bao nhiêu ẩn? "))
#     A = np.random.rand(n, n)
#     b = np.random.rand(n).reshape(-1, 1)
#     print("Ma trận mở rộng Ab:")
#     print(np.hstack([A, b]))
#     x = solve_lu(A, b)
#     print("Nghiệm của hệ phương trình:", x)
#     print("Kiểm tra:")
#     print("Ax - b =", np.dot(A, x) - b)

if __name__ == '__main__':
    A = np.array([[0.00111, 20, 4],
         [0, 0.00112, 32],
         [41, 1, 2]])
    b = np.array([37,128.0112, 9]).reshape(-1, 1)
    L, U, x = solve_lu(A, b)
    print("L:", L)
    print("U:", U)
    print("Nghiệm của hệ phương trình:", x)
    print("Kiểm tra:")
    print("Ax - b =", np.dot(A, x.reshape(-1, 1)) - b)
    
