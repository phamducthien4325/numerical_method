from typing import Optional
import math

# Con goi la phuong phap lap don
def fixed_point_iteration(g,
                          p0: float,
                          tol: float = 1e-5,
                          max_iter: int = 10000) -> Optional[float]:
    i = 1
    while i < max_iter:
        p = g(p0)
        if abs(p - p0) < tol:
            print("Sau n = {}".format(i))
            return p
        p0 = p
        i += 1
    print("Khong tim thay nghiem sau n la {} va p la {}".format(i, p0))
    return None

def g1(x: float) -> float:
    return 5 ** (-x)

if __name__ == "__main__":
    print(fixed_point_iteration(g1, 4, tol=1e-5, max_iter=10000000))