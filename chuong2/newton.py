from typing import Optional
import math

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

if __name__ == "__main__":
    print(newton(g1, g1_prime, math.pi / 2, tol=1e-5, max_iter=1000000))