import sympy as sp
import math

def rk2(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> list[tuple[float, float]]:
    xi = x0
    yi = y0
    result = [(xi, yi)]

    step = round((xn - x0) / h)

    for _ in range(int(step)):
        k1 = y_prime.subs({x: xi, y: yi}).evalf()
        k2 = y_prime.subs({x: xi + h, y: yi + k1 * h}).evalf()
        yi = yi + h / 2 * (k1 + k2)
        xi += h
        result.append((xi, yi))
    
    return result

def rk3(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> list[tuple[float, float]]:
    xi = x0
    yi = y0
    result = [(xi, yi)]

    step = round((xn - x0) / h)

    for _ in range(int(step)):
        k1 = y_prime.subs({x: xi, y: yi}).evalf()
        k2 = y_prime.subs({x: xi + h / 2, y: yi + h / 2 * k1}).evalf()
        k3 = y_prime.subs({x: xi + h, y: yi - h * k1 + 2 * h * k2}).evalf()
        yi = yi + h / 6 * (k1 + 4 * k2 + k3)
        xi += h
        result.append((xi, yi))

    return result

def rk4(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> list[tuple[float, float]]:
    xi = x0
    yi = y0
    result = [(xi, yi)]

    step = round((xn - x0) / h)
    
    for _ in range(int(step)):
        k1 = y_prime.subs({x: xi, y: yi}).evalf()
        k2 = y_prime.subs({x: xi + h / 2, y: yi + 1 / 2 * k1 * h}).evalf()
        k3 = y_prime.subs({x: xi + h / 2, y: yi + 1 / 2 * k2 * h}).evalf()
        k4 = y_prime.subs({x: xi + h, y: yi + k3 * h}).evalf()
        yi = yi + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        xi += h
        result.append((xi, yi))
    
    return result

if __name__ == '__main__':
    x = sp.symbols('x')
    y = sp.symbols('y')
    ##############################
    y_prime = 4 * y - 4 * x
    x0 = 0
    xn = 1
    y0 = 0
    h = 0.5
    f = x + 1 / (1 - x)
    ############################
    solution = rk4(x0, y0, h, xn, y_prime, x, y)
    print("Solution using Runge Kutta 4 method:")
    for x_val, y_val in solution:
        print(f"x: {x_val}, y: {y_val}", end='; ')
        print(f"Error: {abs(y_val - f.subs(x, x_val).evalf())}")
    # Min, Max = error_bound(f, x, x0, xn)
    # Min = abs(h**2 / 2 * Min)
    # Max = abs(h**2 / 2 * Max)
    # print(f"Error Bound: [{Min}, {Max}]")
