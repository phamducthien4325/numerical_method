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
            print("Sau n = {}".format(i), end="; ")
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
    elif f(b) * f_2prime(b) > 0:
        return b
    else:
        print("a, b khong thoa man x0")
        return None

if __name__ == "__main__":
    x = sp.symbols('x')
    f = x + 1 - 2 * sp.sin(sp.pi * x)
    a = 0.5
    b = 1
    f_prime = sp.diff(f, x)
    f_2prime = sp.diff(f_prime, x)
    print("Phuong phap newton:", end="")
    x0 = find_x0(sp.lambdify(x, f), sp.lambdify(x, f_2prime), a, b)
    print(f"tai x0 = {x0}; ", end="")
    print(newton(sp.lambdify(x, f), sp.lambdify(x, f_prime), x0, tol=1e-7, max_iter=1000000))