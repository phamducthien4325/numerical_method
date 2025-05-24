#  Cong thuc Composite Simpson Rule, khi n = 2 thi cong thuc tro thanh Trapezoidal Rule
#  _______________________________________________________________________________________________________________________
# |                      h    [           ⎛(n/2)-1       ⎞       ⎛(n/2)          ⎞        ]    (b - a) * h^4             | ξ ∈ (a, b)
# |∫[a to b] f(x) dx  = --- * [f(a) + 2 * ⎜   Σ    f(x2j)⎟ + 4 * ⎜  Σ  f(x(2j-1))⎟  + f(b)] - --------------- * f''''(ξ) | h = (b - a) / n
# |                      3    [           ⎝  j=1         ⎠       ⎝ j=1           ⎠        ]         180                  | xj = a + j * h
# |______________________________________________________________________________________________________________________| j = 0, 1, 2, ..., n
import sympy as sp
import math
from scipy.optimize import minimize_scalar

def simpson(f, a: float, b: float, n: int) -> float:
    """
    :param f: function to be integrated
    :param a: lower limit of integration
    :param b: upper limit of integration
    :param n: number of subintervals
    :return: approximate value of the integral
    """
    h = (b - a) / n
    return h / 3 * (f(a) + 2 * sum(f(a + 2 * j * h) for j in range(1, n // 2)) +
                   4 * sum(f(a + (2 * j - 1) * h) for j in range(1, n // 2 + 1)) + f(b))

def find_min_max_scipy(f, a: float, b: float) -> float:
    """
    :param f2prime: second derivative of f
    :param a: left bound
    :param b: right bound
    :return: error bound of f'(x0)
    """
    if a > b:
        a, b = b, a
    res_min = minimize_scalar(f, bounds=(a, b), method='bounded')
    res_max = minimize_scalar(lambda x: -f(x), bounds=(a, b), method='bounded')
    if res_min.fun < res_max.fun:
        res_min, res_max = res_max, res_min
    return res_min.fun, res_max.fun

def find_abs_min_max(f, x, a: float, b: float):
    try:
        f_prime = sp.diff(f, x)
        critical_points = sp.solve(f_prime, x)
        zero_points = sp.solve(f, x)
        def is_valid(pt):
            return pt.is_real and a <= float(pt.evalf()) <= b
        domain_points = [pt for pt in critical_points + zero_points if is_valid(pt)]
        domain_points.extend([a, b])
        domain_points = list(set(domain_points))
        values = [(pt, sp.Abs(f).subs(x, pt).evalf()) for pt in domain_points]
        min_point, min_val = min(values, key=lambda t: t[1])
        max_point, max_val = max(values, key=lambda t: t[1])
    except Exception as e:
        f = sp.lambdify(x, f)
        min_val,  max_val = find_min_max_scipy(f, a, b)
        if min_val < 0 and max_val < 0:
            min_val, max_val = max_val, min_val
        elif min_val * max_val < 0:
            min_val, max_val = 0, max(abs(min_val), abs(max_val))
    return min_val, max_val


def error_bound(f, x, a: float, b: float, n: int) -> float:
    """
    :param f: function to be integrated
    :param a: lower limit of integration
    :param b: upper limit of integration
    :param n: number of subintervals
    :return: error bound of the integral
    """
    if a > b:
        a, b = b, a
    h = (b - a) / n
    f_4prime = sp.diff(f, x, 4)
    Min, Max = find_abs_min_max(f_4prime, x, a, b)
    Min = abs(Min)
    Max = abs(Max)
    if Min > Max:
        Min, Max = Max, Min
    return Min, Max

if __name__ == "__main__":
    x = sp.symbols('x')
    f = 1/(x * sp.ln(x))
    a = math.e
    b = math.e + 1
    n = 2
    h = (b - a) / n
    integ = simpson(sp.lambdify(x, f), a, b, n)
    print(f"Simpson Rule: {integ}")
    print(f"Error: {abs(sp.integrate(f, (x, a, b)).evalf() - integ)}")
    Min, Max = error_bound(f, x, a, b, n)
    Min = abs((b - a) * h ** 4 / 180 * Min)
    Max = abs((b - a) * h ** 4 / 180 * Max)
    print(f"Error Bound: [{Min:.10f}, {Max:.10f}]")
