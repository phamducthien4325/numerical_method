import sympy as sp
from scipy.optimize import minimize_scalar
import math

def modified_euler_method(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> list[tuple[float, float]]:
    xi = x0
    yi = y0
    result = [(xi, yi)]
    
    step = (xn - x0) / h

    for _ in range(int(step)):
        k1 = y_prime.subs({x: xi, y: yi}).evalf()
        y_star = yi + h * k1
        k2 = y_prime.subs({x: xi + h, y: y_star}).evalf()
        yi += h / 2 * (k1 + k2)
        xi += h
        result.append((xi, yi))
    
    return result

if __name__ == '__main__':
    x = sp.symbols('x')
    y = sp.symbols('y')
    ##############################
    y_prime = (-x * y + 4 * x / y)
    x0 = 0
    xn = 1
    y0 = 1
    h = 0.1
    f = (4 - 3 * sp.exp(-x ** 2)) ** (1 / 2)
    ############################
    solution = modified_euler_method(x0, y0, h, xn, y_prime, x, y)
    print("Solution using Modified Euler's method:")
    for x_val, y_val in solution:
        print(f"x: {x_val}, y: {y_val}", end='; ')
        print(f"Error: {abs(y_val - f.subs(x, x_val).evalf())}")
    # Min, Max = error_bound(f, x, x0, xn)
    # Min = abs(h**2 / 2 * Min)
    # Max = abs(h**2 / 2 * Max)
    # print(f"Error Bound: [{Min}, {Max}]")
