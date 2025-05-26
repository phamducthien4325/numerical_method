import sympy as sp
import math
from runge_kutta4 import rk4

def adam_bashforth1(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> float:
    xi = x0
    yi = y0

    step = (xn - x0) / h - 0
    
    for _ in range(int(step)):
        k1 = y_prime.subs({x: xi, y: yi}).evalf()
        yi = yi + h * k1
        xi += h
    
    return yi

def adam_bashforth2(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> float:
    init = rk4(x0, y0, h, x0 + h, y_prime, x, y)
    yi = [init[0][1], init[1][1]]
    xi = x0 + h

    step = (xn - x0) / h - 1

    for _ in range(int(step)):
        y_tmp = yi[1] + h / 2 * (3 * y_prime.subs({x: xi, y: yi[1]}).evalf() -
                y_prime.subs({x: xi - h, y: yi[0]}).evalf())
        yi[0] = yi[1]
        yi[1] = y_tmp
        xi += h

    return yi[1]


def adam_bashforth3(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> float:
    init = rk4(x0, y0, h, x0 + 2 * h, y_prime, x, y)
    yi = [init[0][1], init[1][1], init[2][1]]
    xi = x0 + 2 * h

    step = (xn - x0) / h - 2

    for _ in range(int(step)):
        y_tmp = yi[2] + h / 12 * (23 * y_prime.subs({x: xi, y: yi[2]}).evalf() -
                16 * y_prime.subs({x: xi - h, y: yi[1]}).evalf() + 
                5 * y_prime.subs({x: xi - 2 * h, y: yi[0]}).evalf())
        yi[0] = yi[1]
        yi[1] = yi[2]
        yi[2] = y_tmp
        xi += h

    return yi[2]

def adam_bashforth4(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> float:
    init = rk4(x0, y0, h, x0 + 3 * h, y_prime, x, y)
    yi = [init[0][1], init[1][1], init[2][1], init[3][1]]
    xi = x0 + 3 * h

    step = (xn - x0) / h - 3

    for _ in range(int(step)):
        y_tmp = yi[3] + h / 24 * (55 * y_prime.subs({x: xi, y: yi[3]}).evalf() -
                59 * y_prime.subs({x: xi - h, y: yi[2]}).evalf() + 
                37 * y_prime.subs({x: xi - 2 * h, y: yi[1]}).evalf() - 
                9 * y_prime.subs({x: xi - 3 * h, y: yi[0]}).evalf())
        yi[0] = yi[1]
        yi[1] = yi[2]
        yi[2] = yi[3]
        yi[3] = y_tmp
        xi += h

    return yi[3]



if __name__ == '__main__':
    x = sp.symbols('x')
    y = sp.symbols('y')
    ##############################
    y_prime = y - x ** 2 + 1
    x0 = 0
    xn = 1
    y0 = 0.5
    h = 0.2
    # f = 1 / 5 * x * sp.exp(3 * x) - 1 / 25 * sp.exp(3 * x) + 1 / 25 * sp.exp(-2 * x)
    ############################
    ad1 = adam_bashforth1(x0, y0, h, xn, y_prime, x, y)
    ad2 = adam_bashforth2(x0, y0, h, xn, y_prime, x, y)
    ad3 = adam_bashforth3(x0, y0, h, xn, y_prime, x, y)
    ad4 = adam_bashforth4(x0, y0, h, xn, y_prime, x, y)
    print("Solution using Adam-Bashforth methods:")
    print(f"1st order: {ad1}, Error: {abs(ad1 - f.subs(x, xn).evalf())}")
    print(f"2nd order: {ad2}, Error: {abs(ad2 - f.subs(x, xn).evalf())}")
    print(f"3rd order: {ad3}, Error: {abs(ad3 - f.subs(x, xn).evalf())}")
    print(f"4th order: {ad4}, Error: {abs(ad4 - f.subs(x, xn).evalf())}")
    # print(f"1st order: {ad1}")
    # print(f"2nd order: {ad2}")
    # print(f"3rd order: {ad3}")
    # print(f"4th order: {ad4}")