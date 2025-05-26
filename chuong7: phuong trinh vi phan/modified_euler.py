import sympy as sp
from scipy.optimize import minimize_scalar
import math

def modified_euler_method(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> list[tuple[float, float]]:
    """
    Euler's method for solving ordinary differential equations.
    
    :param x0: initial x value
    :param y0: initial y value
    :param h: step size
    :param xn: final x value
    :param y_prime: the derivative of y with respect to x (dy/dx)
    :param x: sympy symbol for x
    :param y: sympy symbol for y
    :return: list of tuples (x, y) representing the solution
    """
    xi = x0
    yi = y0
    result = [(xi, yi)]
    
    while xi < xn:
        k1 = y_prime.subs({x: xi, y: yi}).evalf()
        y_star = yi + h * k1
        k2 = y_prime.subs({x: xi + h, y: y_star}).evalf()
        yi += h / 2 * (k1 + k2)
        xi += h
        result.append((xi, yi))
    
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


# def error_bound(f, x, a: float, b: float) -> tuple[float, float]:
#     """
#     :param f: function to be integrated
#     :param x: sympy symbol for the variable
#     :param x0: point at which to evaluate the error bound
#     :param xi: list of x values for interpolation
#     :return: tuple of minimum and maximum absolute error bounds
#     """
#     f_2prime = sp.diff(f, x, 2)
#     Min, Max = find_abs_min_max(f_2prime, x, a, b)
#     Min = abs(Min)
#     Max = abs(Max)
#     if Min > Max:
#         Min, Max = Max, Min
#     return Min, Max

if __name__ == '__main__':
    x = sp.symbols('x')
    y = sp.symbols('y')
    y_prime = x * sp.exp(3 * x) - 2 * y
    x0 = 0
    xn = 1
    y0 = 0
    h = 0.5
    f = 1 / 5 * x * sp.exp(3 * x) - 1 / 25 * sp.exp(3 * x) + 1 / 25 * sp.exp(-2 * x)
    solution = modified_euler_method(x0, y0, h, xn, y_prime, x, y)
    print("Solution using Euler's method:")
    for x_val, y_val in solution:
        print(f"x: {x_val}, y: {y_val}", end='; ')
        print(f"Error: {abs(y_val - f.subs(x, x_val).evalf())}")
    # Min, Max = error_bound(f, x, x0, xn)
    # Min = abs(h**2 / 2 * Min)
    # Max = abs(h**2 / 2 * Max)
    # print(f"Error Bound: [{Min}, {Max}]")
