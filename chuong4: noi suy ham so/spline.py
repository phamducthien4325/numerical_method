import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

import numpy as np
from chuong3 import gauss

def find_k(fi: list[float], k0: float, k4: float, h: float) -> list[float]:
    """
    :param fi: list of function values at xi
    :param k0: value of the first derivative at the first point
    :param k4: value of the last derivative at the last point
    :param h: step size (assumed constant)
    :return: list of second derivatives at xi
    """
    A = np.array([
        [4, 1, 0],
        [1, 4, 1],
        [0, 1, 4],
    ], dtype=float)
    b = 3 / h * np.array([[fi[2] - fi[0]],
                          [fi[3] - fi[1]],
                          [fi[4] - fi[2]]], dtype=float) - np.array([[k0], [0], [k4]], dtype=float)
    k = [k0, 0, 0, 0, k4]
    k[1:4] = gauss(A, b).flatten().tolist()
    return k

def cubic_spline(fi: list[float], k: list[float], h: float):
    """
    :param fi: list of function values at xi
    :param k: list of second derivatives at xi
    :param h: step size (assumed constant)
    :return: list of cubic spline coefficients
    """
    a = [[0] * 4 for _ in range(4)]
    for j in range(4):
        a[j][0] = fi[j]
        a[j][1] = k[j]
        a[j][2] = 3 / h**2 * (fi[j + 1] - fi[j]) - (2 * k[j] + k[j + 1]) / h
        a[j][3] = 2 / h**3 * (fi[j] - fi[j + 1]) + (k[j] + k[j + 1]) / h**2
    return a

if __name__ == '__main__':
    k0 = 0
    kn = 0
    xi = [-2, -1, 0, 1, 2]
    fi = [0, 0, 1, 0, 0]
    k = find_k(fi, k0, kn, xi[1] - xi[0])
    a = cubic_spline(fi, k, xi[1] - xi[0])
    print('Cac he so cua spline cubic la:')
    for i in range(len(a)):
        print(f'a[{i}] = {a[i]}')