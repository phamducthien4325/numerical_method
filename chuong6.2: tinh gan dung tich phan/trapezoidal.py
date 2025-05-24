#  Cong thuc Composite Trapezoidal Rule, khi n = 1 thi cong thuc tro thanh Trapezoidal Rule
#  ____________________________________________________________________________________
# |                      h                 n-1               (b - a) * h^2            | ξ ∈ (a, b)
# |∫[a to b] f(x) dx  = --- * [f(a) + 2 *   Σ f(xj)+ f(b)] - --------------- * f''(ξ) | h = (b - a) / n
# |                      2                 j=1                     12                 | xj = a + j * h
# |___________________________________________________________________________________| j = 0, 1, 2, ..., n
import sympy as sp

def trapezoidal(f, a: float, b: float, n: int) -> float:
    """
    :param f: function to be integrated
    :param a: lower limit of integration
    :param b: upper limit of integration
    :param n: number of subintervals
    :return: approximate value of the integral
    """
    h = (b - a) / n
    return h / 2 * (f(a) + 2 * sum(f(a + j * h) for j in range(1, n)) + f(b))

import sympy as sp

def find_abs_min_max(f, x, a: float, b: float):
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
    f_2prime = sp.diff(f, x, 2)
    Min, Max = find_abs_min_max(f_2prime, x, a, b)
    Min = abs(Min)
    Max = abs(Max)
    if Min > Max:
        Min, Max = Max, Min
    return Min, Max

if __name__ == "__main__":
    x = sp.symbols('x')
    f = x ** 2 * sp.exp(-x)
    a = 0
    b = 1
    n = 1
    h = (b - a) / n
    integ = trapezoidal(sp.lambdify(x, f), a, b, n)
    print(f"Trapezoidal Rule: {integ}")
    print(f"Error: {abs(sp.integrate(f, (x, a, b)).evalf() - integ)}")
    Min, Max = error_bound(f, x, a, b, n)
    Min = abs((b - a) * h ** 2 / 12 * Min)
    Max = abs((b - a) * h ** 2 / 12 * Max)
    print(f"Error Bound: [{Min:.7f}, {Max:.7f}]")
