#  _________________________________________
# |                 f(xi)*(xi - x(i - 1))  | 
# |x(i + 1) = xi - ---------------------   | if f(xi)*f(x(i-1)) < 0 x(i-1) = x(i-1)
# |                 f(xi) - f(x(i - 1))    | else x(i-1) = x(i-2)
# |________________________________________| ap dung voi i >= 2

from typing import Optional
import math
import sympy as sp

def diemsai(f,
            p0: float,
            p1: float,
            tol: float = 1e-5,
            max_iter: int = 10000) -> Optional[float]:
    i = 1
    q0 = f(p0)
    q1 = f(p1)
    while i < max_iter:
        p = p1 - q1 * (p1 - p0) / (q1 - q0)
        if abs(p - p1) < tol:
            print("Sau n = {}".format(i), end="; ")
            return p
        i += 1
        q = f(p)
        if q * q1 < 0:
            p0 = p1
            q0 = q1
        p1 = p
        q1 = q
    print("Khong tim thay nghiem sau n la {} va p la {}".format(i - 1, p1))
    return None

if __name__ == "__main__":
    x = sp.symbols('x')
    f = x + 1 - 2 * sp.sin(sp.pi * x)
    a = 0.5
    b = 1
    print("Phuong phap diem sai:", end="")
    print(diemsai(sp.lambdify(x, f), a, b, tol=1e-7, max_iter=1000000))