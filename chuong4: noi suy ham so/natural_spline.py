import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

import numpy as np
from chuong3 import gauss

def natural_cubic_spline(x, a):
    n = len(x) - 1  # number of intervals
    h = [x[i + 1] - x[i] for i in range(n)]

    # Step 2: compute alpha
    alpha = [0] * n
    for i in range(1, n):
        alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1])

    # Step 3: initialize l, mu, z
    l = [1] + [0] * n
    mu = [0] * (n + 1)
    z = [0] * (n + 1)

    # Step 4: forward sweep
    for i in range(1, n):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    # Step 5: set last values
    l[n] = 1
    z[n] = 0
    c = [0] * (n + 1)
    b = [0] * n
    d = [0] * n

    # Step 6: backward sweep
    for j in range(n - 1, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    # Step 7: output
    return [(a[j], b[j], c[j], d[j]) for j in range(n)]


if __name__ == '__main__':
    xi = [0, 2, 4, 6]
    fi = [2, -2, 2, 78] 
    n = len(xi) - 1
    a = cubic_spline(fi, k, xi[1] - xi[0])
    print('Cac he so cua spline cubic la:')
    for i in range(len(a)):
        print(f'a[{i}] = {a[i]}')