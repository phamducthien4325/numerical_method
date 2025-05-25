import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

import numpy as np
from chuong3 import gauss

def find_k(fi: list[float], k0: float, kn: float, h: float) -> list[float]:
    """
    :param fi: list of function values at xi
    :param k0: value of the first derivative at the first point
    :param k4: value of the last derivative at the last point
    :param h: step size (assumed constant)
    :return: list of second derivatives at xi
    """
    n = len(fi) - 1
    A = np.zeros((n - 1, n - 1), dtype=float)
    b = np.zeros((n - 1, 1), dtype=float)
    k = [0] * (n + 1)
    for i in range(0, n - 1):
        if i == 0:
            A[0][0] = 4
            A[0][1] = 1
        elif i == n - 2:
            A[i][i - 1] = 1
            A[i][i] = 4
        else:
            A[i][i - 1] = 1
            A[i][i] = 4
            A[i][i + 1] = 1
    for i in range(0, n - 1):
        b[i] = 3 / h * (fi[i + 2] - fi[i])
    b[0] -= k0
    b[n - 2] -= kn
    k[1:n] = gauss(A, b).flatten().tolist()
    k[0] = k0
    k[n] = kn
    return k

def cubic_spline(fi: list[float], k: list[float], h: float):
    """
    :param fi: list of function values at xi
    :param k: list of second derivatives at xi
    :param h: step size (assumed constant)
    :return: list of cubic spline coefficients
    """
    n = len(fi) - 1
    a = [[0] * 4 for _ in range(n)]
    for j in range(n):
        a[j][0] = fi[j]
        a[j][1] = k[j]
        a[j][2] = 3 / h**2 * (fi[j + 1] - fi[j]) - (2 * k[j] + k[j + 1]) / h
        a[j][3] = 2 / h**3 * (fi[j] - fi[j + 1]) + (k[j] + k[j + 1]) / h**2
    return a

if __name__ == '__main__':
    xi = [0, 2, 4, 6]
    fi = [2, -2, 2, 78] 
    k0 = 0
    kn = 0
    n = len(xi) - 1
    k = find_k(fi, k0, kn, xi[1] - xi[0])
    a = cubic_spline(fi, k, xi[1] - xi[0])
    print('Cac he so cua spline cubic la:')
    for i in range(len(a)):
        print(f'a[{i}] = {a[i]}')