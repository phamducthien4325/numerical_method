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
            print(f"n={i}", end = ';')
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
    return  + 1- 2 * math.sin(math.pi * x)

if __name__ == "__main__":
    print("Phuong phap chia doi:", end="")
    res = bisection(func, 0, 1.5, tol=1e-7, max_iter=10000000)
    if res != None:
        print(res)
    else:
        print("none")