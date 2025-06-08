import numpy as np

def power_method(A, x, TOL, N):
    n = A.shape[0]
    k = 1
    p = np.argmax(np.abs(x))
    x = x / x[p]
    
    while k <= N:
        y = A @ x
        mu = y[p]
        p = np.argmax(np.abs(y))
        if y[p] == 0:
            return "A has the eigenvalue 0, select a new vector x and restart", x
        ERR = np.linalg.norm(x - y / y[p], ord=np.inf)
        x = y / y[p]
        if ERR < TOL:
            return mu, x
        k += 1
    return "The maximum number of iterations was exceeded", None

if __name__ == "__main__":
    A = np.array([[4, 1],
                  [2, 3]], dtype=float)
    x = np.array([1, 1], dtype=float)
    TOL = 1e-6
    N = 100

    result = power_method(A, x, TOL, N)
    print("Result:", result)