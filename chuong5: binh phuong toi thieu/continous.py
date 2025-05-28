import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from chuong3 import gauss
import numpy as np
from scipy.integrate import quad
import math

def error(fx, a: float, b: float, px: list) -> float:
    """
    :param fx: ham so
    :param a: diem bat dau
    :param b: diem ket thuc
    :param px: danh sach cac he so cua da thuc
    :return: sai so
    """
    def p(px, x):
        return sum(px[j] * x ** j for j in range(len(px)))
    
    lamb = lambda x: (fx(x) - p(px, x)) ** 2
    integral, _ = quad(lamb, a, b)
    return integral

def continous(n: int, f, a: float, b: float) -> list:
    A = [[0] * (n + 1) for _ in range(n + 1)]
    B = [0] * (n + 1)
    for j in range(n + 1):
        for k in range(n + 1):
            A[j][k] = (b ** (j + k + 1) - a ** (j + k + 1)) / (j + k + 1)
        lamb = lambda x: f(x) * x ** j
        B[j], error = quad(lamb, a, b)
    B = np.array(B).reshape(-1, 1)
    A = np.array(A)
    return gauss(A, B)

def fx(x):
    return math.log(x + 2)
    
if __name__ == '__main__':
    # n = int(input('Nhap bac n cua p(x): '))
    # a = float(input('Nhap a: '))
    # b = float(input('Nhap b: '))
    n = 2
    a = -1
    b = 1
    px = continous(n, fx, a, b)
    print('Da thuc p(x) = ', end='')
    for i in range(n + 1):
        if i == 0:
            print(f'{px[i]:.4f} + ', end='')
        elif i == n:
            print(f'{px[i]:.4f} * x^{i}', end='')
        else:
            print(f'{px[i]:.4f} * x^{i} + ', end='')
    print(f"; error = {error(fx, a, b, px)}") 