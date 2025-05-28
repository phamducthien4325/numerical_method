import sympy as sp
import math
from runge_kutta4 import rk4


def adam_moulton4(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> float:
    init = rk4(x0, y0, h, x0 + 3 * h, y_prime, x, y)
    yi = [init[0][1], init[1][1], init[2][1], init[3][1]]
    xi = x0 + 3 * h

    step = round((xn - x0) / h - 3)

    for _ in range(int(step)):
        y_pred = yi[3] + h / 24 * (55 * y_prime.subs({x: xi, y: yi[3]}).evalf() -
                59 * y_prime.subs({x: xi - h, y: yi[2]}).evalf() + 
                37 * y_prime.subs({x: xi - 2 * h, y: yi[1]}).evalf() - 
                9 * y_prime.subs({x: xi - 3 * h, y: yi[0]}).evalf())
        y_corr = yi[3] + h / 24 * (9 * y_prime.subs({x: xi + h, y: y_pred}).evalf() +
                19 * y_prime.subs({x: xi, y: yi[3]}).evalf() -
                5 * y_prime.subs({x: xi - h, y: yi[2]}).evalf() +
                y_prime.subs({x: xi - 2 * h, y: yi[1]}).evalf())
        yi[0] = yi[1]
        yi[1] = yi[2]
        yi[2] = yi[3]
        yi[3] = y_corr
        xi += h

    return yi[3]



if __name__ == '__main__':
    x = sp.symbols('x')
    y = sp.symbols('y')
    ##############################
    y_prime = x * sp.exp(3 * x) - 2 * x
    x0 = 0
    xn = 1
    y0 = 0
    h = 0.2
    f = 1 / 5 * x * sp.exp(3 * x) - 1/ 25 * sp.exp(3 * x) + 1 / 25 * sp.exp(-2 * x)
    ############################
    res = adam_moulton4(x0, y0, h, xn, y_prime, x, y)
    print("Solution using Adam-Moulton methods:")
    print(f"4th order: {res}, Error: {abs(res - f.subs(x, xn).evalf())}")
    # print(f"1st order: {ad1}")
    # print(f"2nd order: {ad2}")
    # print(f"3rd order: {ad3}")
    # print(f"4th order: {ad4}")