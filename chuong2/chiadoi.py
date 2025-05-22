from typing import Optional
import math

def bisection(f,
              a: float, 
              b: float,
              tol: float = 1e-5,
              max_iter: int = 10000) -> Optional[float]:
    i = 1
    while i <= max_iter:
        p = a + (b - a) / 2
        if f(p) == 0 or (b - a) / 2 < tol:
            return p
        elif f(p) * f(a) < 0:
            b = p
        else:
            a = p
        i += 1
    print(f"Number of iterations: {i - 1}")
    print(f"Final approximation: {p}")
    return None

def func(x: float) -> float:
    return x ** 3 - x - 1

if __name__ == "__main__":
    res = bisection(func, 1, 2, tol=1e-4, max_iter=10000000)
    if res != None:
        print(res)
    else:
        print("none")