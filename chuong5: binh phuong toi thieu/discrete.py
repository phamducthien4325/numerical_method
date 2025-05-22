import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from chuong3 import gauss
import numpy as np

def discrete(n: int, m: int, x: list, fx: list) -> list:
    if n == 0:
        n = 1
        fx = np.log(fx)
    elif n == -1:
        n = 1
        fx = np.log(fx)
        x = np.log(x)
    A = [[0] * (n + 1) for _ in range(n + 1)]
    b = [0] * (n + 1)
    for i in range(n + 1):
        for j in range(n + 1):
            A[i][j] = sum([x[k] ** (i + j) for k in range(m)])
        b[i] = sum([fx[k] * x[k] ** i for k in range(m)])
    # Giai he phuong trinh Ax = b
    A = np.array(A)
    b = np.array(b).reshape(-1, 1)
    return gauss(A, b)

def error(n: int, m: int, x: list, fx: list, px: list) -> float:
    if n == 0:
        return np.sum((px[0] * np.exp(px[1] * np.array(x, dtype=float)) - fx) ** 2)
    elif n == -1:
        return np.sum((px[0] * np.array(x, dtype=float) ** px[1] - fx) ** 2)
    err = 0
    for i in range(m):
        err += abs(fx[i] - sum([px[j] * x[i] ** j for j in range(n + 1)])) ** 2
    return err

if __name__ == '__main__':

    x = [2, 4, 6]
    fx = [1.7, 4.1, 7]
    # if n = 0: form b*e^(ax)
    # if n = -1: form b*x^a
    # if n > 0 polynomial
    n = 1
    m = 3
    # n = int(input('Nhap bac n cua p(x): '))
    # m = int(input('Nhap so luong mau: '))
    # x = []
    # fx = []
    # x = list(map(float, input(f'Nhap {m} mau x: ').split()))
    # while len(x) != m:
    #     print('So luong mau x khong hop le, vui long nhap lai!')
    #     x = list(map(float, input(f'Nhap {m} mau x: ').split()))
    # fx = list(map(float, input(f'Nhap {m} mau f(x): ').split()))
    # while len(fx) != m:
    #     print('So luong mau f(x) khong hop le, vui long nhap lai!')
    #     fx = list(map(float, input(f'Nhap {m} mau f(x): ').split()))
    px = discrete(n, m, x, fx)
    if n == 0:
        px[0] = np.exp(px[0])
    elif n == -1:
        px[0] = np.exp(px[0])
    print('Da thuc p(x) = ', end='')
    if n > 0:
        for i in range(n + 1):
            if i == 0:
                print(f'{px[i]:.4f} + ', end='')
            elif i == n:
                print(f'{px[i]:.4f} * x^{i}', end='')
            else:
                print(f'{px[i]:.4f} * x^{i} + ', end='')
        print()
    elif n == 0:
        print(f'{px[0]:.4f} * e^({px[1]:.4f} * x)', end='')
        print()
    elif n == -1:
        print(f'{px[0]:.4f} * x^{px[1]:.4f}', end='')
        print()
    print('Sai so: ', end='')
    print(f'{error(n, m, x, fx, px):.4f}')