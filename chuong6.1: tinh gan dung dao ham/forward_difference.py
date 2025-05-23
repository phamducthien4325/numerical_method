from scipy.optimize import minimize_scalar
import math

#  ____________________________________________
# |          f(x0 + h) - f(x0)     h           |
# |f'(x0) = ------------------- - --- * f''(ξ) |
# |                  h             2           | # với ξ ∈ (x0, x0 + h)
# |____________________________________________|
# neu h < 0 thi goi la backward difference
def forward_difference(f, x0: float, h: float) -> float:
    """
    :param f: function to be derived
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (f(x0 + h) - f(x0)) / h
    
def forward_difference_2(fx0: float, fxh: float, x0: float, h: float) -> float:
    """
    :param fx0: f(x0)
    :param fxh: f(x0 + h)
    :param x0: point to be derived
    :param h: step size
    :return: derivative of f at x0
    """
    return (fxh - fx0) / h

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
    return 2 * x * math.log(x) + x

def abs_f2prime(x):
    """
    :param x: point to be derived
    :return: |f''(x)|
    """
    return abs(2 * math.log(x) + 3)

def error(f, x0: float, fx0: float) -> float:
    """
    :param f: function to be derived
    :param x0: point to be derived
    :param fx0: f(x0)
    :return: error of f'(x0)
    """
    return abs(f(x0) - fx0)

def error_bound(f2prime, a: float, b: float) -> float:
    """
    :param f2prime: second derivative of f
    :param a: left bound
    :param b: right bound
    :return: error bound of f'(x0)
    """
    if a > b:
        a, b = b, a
    res_min = minimize_scalar(f2prime, bounds=(a, b), method='bounded')
    res_max = minimize_scalar(lambda x: -f2prime(x), bounds=(a, b), method='bounded')
    return res_min.fun, res_max.fun

if __name__ == '__main__':
    x = [1.0, 1.2, 1.4]
    fx = [1.0000, 1.2625, 1.6595]
    h = x[1] - x[0]
    for i in range(len(x) - 1):
        print("Forward difference:")
        f_prime = forward_difference_2(fx[i], fx[i + 1], x[i], h)
        print(f"f'({x[i]}) = {f_prime}")
        print(f"Error: {error(fprime, x[i], f_prime):.4f}")
        a, b = error_bound(abs_f2prime, x[i], x[i] + h)
        a = abs(h / 2 * a)
        b = abs(h / 2 * b)
        print(f"Error bound: [{a:.4f}, {b:.4f}]")
    for i in range(len(x) - 1, 0, -1):
        print("Backward difference:")
        f_prime = forward_difference_2(fx[i], fx[i - 1], x[i], -h)
        print(f"f'({x[i]}) = {f_prime}")
        print(f"Error: {error(fprime, x[i], f_prime):.4f}")
        a, b = error_bound(abs_f2prime, x[i] - h, x[i])
        a = abs(h / 2 * a)
        b = abs(h / 2 * b)
        print(f"Error bound: [{a:.4f}, {b:.4f}]")