#  __________________________
# |                 f(xi)   | 
# |x(i + 1) = xi - -------- | chon x0 sao cho f(x0) * f''(x0) < 0
# |                 f'(xi)  | 
# |_________________________| 

from typing import Optional
import math
import sympy as sp

# Con goi la phuong phap tiep tuyen
def newton(f,
           f_prime,
           p0: float,
           tol: float = 1e-5,
           max_iter: int = 10000) -> Optional[float]:
    i = 1
    while i <= max_iter:
        if f_prime(p0) == 0:
            print("f'(p0) = 0, try another p0")
            return None
        p = p0 - f(p0) / f_prime(p0)
        if abs(p - p0) < tol:
            print("Sau n = {}".format(i))
            return p
        p0 = p
        i += 1
    print("Khong tim thay nghiem sau n la {} va p la {}".format(i - 1, p0))
    return None

def g1(x: float) -> float:
    return x - 0.8 - 0.2 * math.sin(x)

def g1_prime(x: float) -> float:
    return 1 - 0.2 * math.cos(x)

def find_x0(f, f_2prime, a: float, b: float) -> float:
    if f(a) * f_2prime(a) > 0:
        return a
    else:
        return b

if __name__ == "__main__":
    x = sp.symbols('x')
    f = x ** 3 - 2 * x ** 2 - 5
    a = 1
    b = 4
    f_prime = sp.diff(f, x)
    f_2prime = sp.diff(f_prime, x)
    x0 = find_x0(sp.lambdify(x, f), sp.lambdify(x, f_2prime), a, b)
    print(newton(sp.lambdify(x, f), sp.lambdify(x, f_prime), x0, tol=1e-5, max_iter=1000000))