from typing import Optional
import math

def daycung(f,
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
            print("Sau n = {}".format(i))
            return p
        p0 = p1
        q0 = q1
        p1 = p
        q1 = f(p)
        i += 1
    print("Khong tim thay nghiem sau n la {} va p la {}".format(max_iter, p1))
    return None

def g1(x: float) -> float:
    return - x ** 3 - math.cos(x)

if __name__ == "__main__":
    print(daycung(g1, -1, 0, tol=1e-4, max_iter=3))