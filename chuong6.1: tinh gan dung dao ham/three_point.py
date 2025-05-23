from scipy.optimize import minimize_scalar
import math

# End point formula
#  ___________________________________________________________________________
# |            1                                              h^2            |
# |f'(x0) = ------- * [-3f(x0) + 4f(x0 + h) - f(x0 + 2h)] +  ----- * f'''(ξ) | ξ ∈ (x0, x0 + 2h)
# |          2 * h                                             3             |
# |__________________________________________________________________________|

# Mid point formula(most accurate)
#  _______________________________________________________________
# |            1                                  h^2            |
# |f'(x0) = ------- * [-f(x0 - h) + f(x0 + h)] - ----- * f'''(ξ) | ξ ∈ (x0 - h, x0 + h)
# |          2 * h                                 6             |
# |______________________________________________________________|

def end_point(f, x0: float, h: float) -> float:
    """
    :param f: function to be derived
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (1 / (2 * h)) * (-3 * f(x0) + 4 * f(x0 + h) - f(x0 + 2 * h))

def mid_point(f, x0: float, h: float) -> float:
    """
    :param f: function to be derived
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (1 / (2 * h)) * (-f(x0 - h) + f(x0 + h))

def end_point_2(fx0: float, fxh: float, fx2h: float, x0: float, h: float) -> float:
    """
    :param fx0: f(x0)
    :param fxh: f(x0 + h)
    :param fx2h: f(x0 + 2h)
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (1 / (2 * h)) * (-3 * fx0 + 4 * fxh - fx2h)

def mid_point_2(fxmh: float, fxh: float, x0: float, h: float) -> float:
    """
    :param fx0: f(x0)
    :param fxh: f(x0 + h)
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (1 / (2 * h)) * (-fxmh + fxh)

def fx(x):
    """
    :param x: point to be derived
    :return: f(x)
    """
    return x ** 2 * math.log(x) + 1

def fprime(x):
    """
    :param x: point to be derived
    :return: f'(x)
    """
    return 2 * math.e ** (2 * x)

def abs_f3prime(x):
    """
    :param x: point to be derived
    :return: |f''(x)|
    """
    return 8 * math.e ** (2 * x)

def error(f, x0: float, fx0: float) -> float:
    """
    :param f: function to be derived
    :param x0: point to be derived
    :param fx0: f(x0)
    :return: error of f'(x0)
    """
    return abs(f(x0) - fx0)

def error_bound(f3prime, a: float, b: float) -> float:
    """
    :param f2prime: second derivative of f
    :param a: left bound
    :param b: right bound
    :return: error bound of f'(x0)
    """
    if a > b:
        a, b = b, a
    res_min = minimize_scalar(f3prime, bounds=(a, b), method='bounded')
    res_max = minimize_scalar(lambda x: -f3prime(x), bounds=(a, b), method='bounded')
    return res_min.fun, res_max.fun

if __name__ == '__main__':
    x = [1.1, 1.2, 1.3, 1.4]
    fx = [9.025013, 11.02318, 13.46374, 16.44465]
    h = x[1] - x[0]
    for i in range(len(x)):
        x0 = x[i]
        fx0 = fx[i]
        if i == 0:
            fxh = fx[i + 1]
            fx2h = fx[i + 2]
            print("End point formula")
            f_prime = end_point_2(fx0, fxh, fx2h, x0, h)
        elif i == len(x) - 1:
            h = -h
            fxh = fx[i - 1]
            fx2h = fx[i - 2]
            print("End point formula")
            f_prime = end_point_2(fx0, fxh, fx2h, x0, h)
        else:
            fxh = fx[i + 1]
            fxmh = fx[i - 1]
            print("Mid point formula")
            f_prime = mid_point_2(fxmh, fxh, x0, h)
        print(f"f'({x0}) = {f_prime}")
        print(f"Error: {error(fprime, x[i], f_prime)}")
        if i == 0 or i == len(x) - 1:
            a, b = error_bound(abs_f3prime, x0, x0 + 2 * h)
            a = abs(h ** 2 / 3 * a)
            b = abs(h ** 2 / 3 * b)
        else:
            a, b = error_bound(abs_f3prime, x0 - h, x0 + h)
            a = abs(h ** 2 / 6 * a)
            b = abs(h ** 2 / 6 * b)
        print(f"Error bound: [{a:.4f}, {b:.4f}]")