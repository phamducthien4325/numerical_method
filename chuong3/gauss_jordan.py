import numpy as np

def gauss_jordan(A, b):
    n = len(A)
    Ab = np.hstack([A, b])

    for i in range(n):
        if Ab[i, i] == 0:
            for k in range(i + 1, n):
                if Ab[k, i] != 0:
                    Ab[[i, k]] = Ab[[k, i]]  
                    break
                
        Ab[i] = Ab[i] / Ab[i, i]
        for j in range(n):
            if i != j:
                Ab[j] = Ab[j] - Ab[j, i] * Ab[i]

        print(f"Ma trận sau bước {i + 1}:")
        print(Ab)

    return Ab[:, -1]

if __name__ == '__main__':
    # n = int(input("Bạn muốn giải hệ phương trình bao nhiêu ẩn? "))
    # A = np.random.rand(n, n)
    # b = np.random.rand(n).reshape(-1, 1)
    A = np.array([[4, -1, 1], [2, 5, 2], [1, 2, 4]], dtype=float)
    b = np.array([[8], [3], [11]], dtype=float)
    print("Ma trận mở rộng Ab:")
    print(np.hstack([A, b]))
    x = gauss_jordan(A, b)
    print("Nghiệm của hệ phương trình:", x)
    print("Kiểm tra:")
    print("Ax - b =", np.dot(A, x.reshape(-1, 1)) - b)
