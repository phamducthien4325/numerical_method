from scipy.optimize import minimize_scalar
import math

#  ____________________________________________
# |          f(x0 + h) - f(x0)     h           |
# |f'(x0) = ------------------- - --- * f''(ξ) |
# |                  h             2           | # với ξ ∈ (x0, x0 + h)
# |____________________________________________|
# neu h < 0 thi goi la backward difference
def forward_difference(f, x0: float, h: float) -> float:
    return (f(x0 + h) - f(x0)) / h
    
def forward_difference_2(fx0: float, fxh: float, x0: float, h: float) -> float:
    return (fxh - fx0) / h

def fx(x):
    return math.sin(x)

def fprime(x):
    return math.cos(x)

def abs_f2prime(x):
    return math.fabs(-math.sin(x))

def error(f, x0: float, fx0: float) -> float:
    return abs(f(x0) - fx0)

def error_bound(f2prime, a: float, b: float) -> float:
    if a > b:
        a, b = b, a
    res_min = minimize_scalar(f2prime, bounds=(a, b), method='bounded')
    res_max = minimize_scalar(lambda x: -f2prime(x), bounds=(a, b), method='bounded')
    return res_min.fun, res_max.fun

if __name__ == '__main__':
    x0 = 0.0
    fx0 = 0
    h = 0.2
    fxh = 0.7414
    f_prime = forward_difference_2(fx0, fxh, x0, h)
    print(f"f'({x0}) = {f_prime:.4f}")
    print(f"Error: {error(fprime, 0.5, f_prime):.4f}")
    a, b = error_bound(abs_f2prime, x0, x0 + h)
    a = abs(h / 2 * a)
    b = abs(h / 2 * b)
    print(f"Error bound: [{a:.4f}, {b:.4f}]")