from scipy.optimize import minimize_scalar
import math
import sympy as sp

# 5 End point formula
#  _______________________________________________________________________________________________________________
# |            1                                                                                h^4              |
# |f'(x0) = -------- * [-25f(x0) + 48f(x0 + h) - 36f(x0 + 2h) + 16f(x0 + 3h) - 3f(x0 + 4h)] +  ----- * f'''''(ξ) | ξ ∈ (x0, x0 + 4h)
# |          12 * h                                                                              5               |
# |______________________________________________________________________________________________________________|

# 5 Mid point formula(most accurate)
#  ______________________________________________________________________________________________
# |            1                                                               h^4              |
# |f'(x0) = -------- * [f(x0 - 2h) - 8f(x0 - h) + 8f(x0 + h) - f(x0 + 2h)] +  ----- * f'''''(ξ) | ξ ∈ (x0 - 2h, x0 + 2h)
# |          12 * h                                                            30               |
# |_____________________________________________________________________________________________|

def five_end_point(f, x0: float, h: float) -> float:
    """
    :param f: function to be derived
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (1 / (12 * h)) * (-25 * f(x0) + 48 * f(x0 + h) - 36 * f(x0 + 2 * h) + 16 * f(x0 + 3 * h) - 3 * f(x0 + 4 * h))

def five_mid_point(f, x0: float, h: float) -> float:
    """
    :param f: function to be derived
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (1 / (12 * h)) * (f(x0 - 2 * h) - 8 * f(x0 - h) + 8 * f(x0 + h) - f(x0 + 2 * h))

def five_end_point_2(fx0: float, fxh: float, fx2h: float, fx3h: float, fx4h: float, x0: float, h: float) -> float:
    """
    :param fx0: f(x0)
    :param fxh: f(x0 + h)
    :param fx2h: f(x0 + 2h)
    :param fx3h: f(x0 + 3h)
    :param fx4h: f(x0 + 4h)
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (1 / (12 * h)) * (-25 * fx0 + 48 * fxh - 36 * fx2h + 16 * fx3h - 3 * fx4h)
    

def five_mid_point_2(fxm2h: float, fxmh: float, fxh: float, fx2h: float, x0: float, h: float) -> float:
    """
    :param fxm2h: f(x0 - 2h)
    :param fxmh: f(x0 - h)
    :param fxh: f(x0 + h)
    :param fx2h: f(x0 + 2h)
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (1 / (12 * h)) * (fxm2h - 8 * fxmh + 8 * fxh - fx2h)

def error(n: float, n_pred: float) -> float:
    """
    :param n: true value
    :param n_pred: predicted value
    :return: error of f'(x0)
    """
    return abs(n - n_pred)

def find_min_max(f, x, a, b):
    f_prime = sp.diff(f, x)
    critical_points = sp.solve(f_prime, x)
    domain_points = [pt for pt in critical_points if pt.is_real and a <= pt <= b]
    domain_points.extend([a, b])
    values = [(pt, f.subs(x, pt).evalf()) for pt in domain_points]
    min_point, min_val = min(values, key=lambda t: t[1])
    max_point, max_val = max(values, key=lambda t: t[1])
    
    return (min_point, min_val), (max_point, max_val)

def error_bound(f5prime, x, a: float, b: float) -> float:
    """
    :param f2prime: second derivative of f
    :param a: left bound
    :param b: right bound
    :return: error bound of f'(x0)
    """
    if a > b:
        a, b = b, a
    mi, ma = find_min_max(f5prime, x, a, b)
    y_min = (mi[1])
    y_max = (ma[1])
    if y_min > y_max:
        y_min, y_max = y_max, y_min
    return y_min, y_max

def f(x):
    """
    :param x: point to be derived
    :return: f(x)
    """
    return math.tan(x)

if __name__ == '__main__':
    x = [-3.0, -2.8, -2.6, -2.4, -2.2, -2.0]
    fx = [16.08554, 12.64465, 9.863738, 7.623176, 5.825013, 4.389056]
    y = sp.symbols('y')
    fy = sp.exp(-y) - 1 + y
    fprime = sp.diff(fy, y)
    f5prime = sp.diff(fy, y, 5)
    h = x[1] - x[0]
    for i in range(len(x)):
        x0 = x[i]
        fx0 = fx[i]
        if i == 0 or i == 1:
            fxh = fx[i + 1]
            fx2h = fx[i + 2]
            fx3h = fx[i + 3]
            fx4h = fx[i + 4]
            print("5 End point formula")
            f_prime = five_end_point_2(fx0, fxh, fx2h, fx3h, fx4h, x0, h)
        elif i == len(x) - 1 or i == len(x) - 2:
            fxh = fx[i - 1]
            fx2h = fx[i - 2]
            fx3h = fx[i - 3]
            fx4h = fx[i - 4]
            print("5 End point formula")
            f_prime = five_end_point_2(fx0, fxh, fx2h, fx3h, fx4h, x0, -h)
        else:
            fxh = fx[i + 1]
            fxmh = fx[i - 1]
            fx2h = fx[i + 2]
            fxm2h = fx[i - 2]
            print("5 Mid point formula")
            f_prime = five_mid_point_2(fxm2h, fxmh, fxh, fx2h, x0, h)
        print(f"f'({x0}) = {f_prime}")
        print(f"Error: {error(fprime.subs({y: x0}), f_prime)}")
        if i == 0 or i == 1:
            a, b = error_bound(f5prime, y, x0, x0 + 4 * h)
            a = abs((h ** 4) / 5 * a)
            b = abs((h ** 4) / 5 * b)
        elif i == len(x) - 1 or i == len(x) - 2:
            a, b = error_bound(f5prime, y, x0 - 4 * h, x0)
            a = abs((h ** 4) / 5 * a)
            b = abs((h ** 4) / 5 * b)
        else:
            a, b = error_bound(f5prime, y, x0 - h, x0 + h)
            a = abs((h ** 4) / 30 * a)
            b = abs((h ** 4)/ 30 * b)
        print(f"Error bound: [{a:.7f}, {b:.7f}]")