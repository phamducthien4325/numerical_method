import sympy as sp
import math
from scipy.optimize import minimize_scalar

def newton_inter(x: list[float],
                 fi: list[float],
                 x0: float) -> float:
    n = len(x)
    fx = fi.copy()
    result = fx[0]
    mul = 1
    for i in range(1, n):
        mul *= (x0 - x[i - 1])
        for j in range(n - i):
            fx[j] = (fx[j + 1] - fx[j]) / (x[j + i] - x[j])
        result += mul * fx[0]
    return result

def newton_forward(x: list[float],
                   fi: list[float],
                   x0: float) -> float:
    n = len(x)
    fx = fi.copy()
    h = x[1] - x[0]
    u = (x0 - x[0]) / h
    result = fx[0]
    mul = 1
    for i in range(1, n):
        mul *= (u - (i - 1)) / i
        for j in range(n - i):
            fx[j] = (fx[j + 1] - fx[j])
        result += mul * fx[0]
    return result

def newton_backward(x: list[float],
                   fi: list[float],
                   x0: float) -> float:
    n = len(x)
    fx = fi.copy()
    h = x[1] - x[0]
    u = (x0 - x[n - 1]) / h
    result = fx[n - 1]
    mul = 1
    for i in range(1, n):
        mul *= (u + (i - 1)) / i
        for j in range(n - 1, i - 1, -1):
            fx[j] = fx[j] - fx[j - 1]
        result += mul * fx[n - 1]
    return result

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


def error_bound(f, x, x0: float, xi: list[float]) -> tuple[float, float]:
    """
    :param f: function to be integrated
    :param x: sympy symbol for the variable
    :param xi: list of x values for interpolation
    :return: tuple of minimum and maximum absolute error bounds
    """
    a, b = min(xi), max(xi)
    n = len(xi) - 1
    f_prime = sp.diff(f, x, n + 1)
    Min, Max = find_abs_min_max(f_prime, x, a, b)
    prod = 1
    for i in range(n + 1):
        prod *= (x0 - xi[i]) / (i + 1)
    Min = abs(Min * prod)
    Max = abs(Max * prod)
    if Min > Max:
        Min, Max = Max, Min
    return Min, Max


if __name__ == '__main__':
    # n = int(input('Nhap bac n cua p(x): '))
    # x = []
    # fx = []
    # for i in range(n + 1):
    #     x.append(float(input(f'x[{i}]: ')))
    #     fx.append(float(input(f'f(x[{i}]): ')))
    # x0 = float(input('Nhap x0: '))
    x = sp.Symbol('x')
    # f = sp.ln(sp.exp(x) + 2)
    xi = [-1.0, -0.5, 0.0, 0.5]
    # fi = [f.subs(x, i).evalf() for i in xi]
    fi = [0.86199480, 0.95802009, 1.0986123, 1.2943767]
    n = 1  # Degree of polynomial
    x0 = 0.43
    for n in range(1, len(xi)):
        newton_val = newton_backward(xi[len(xi) - n - 1:], fi[len(xi) - n - 1:], x0)
        print(f'{n} degree backward f({x0}) = {newton_val}')
        # newton_val = newton_inter(xi[:n + 1], fi[:n + 1], x0)
        # print(f'{n} degree inter f({x0}) = {newton_val}')
        # newton_val = newton_forward(xi[:n + 1], fi[:n + 1], x0)
        # print(f'{n} degree forward f({x0}) = {newton_val}')
        # print(f'Absolute error: {abs(f.subs(x, x0).evalf() - newton_val)}')
        # Min, Max = error_bound(f, x, x0, xi[:n + 1])
        # print(f'Error bound: [{Min}, {Max}]')