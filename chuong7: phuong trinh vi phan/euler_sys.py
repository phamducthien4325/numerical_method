import sympy as sp
from scipy.optimize import minimize_scalar
import math

def euler_system_method(x0: float, Y0: sp.Matrix, h: float, xn: float, F, x, Ysym) -> list[tuple[float, float]]:
    xi = x0
    Yi = Y0
    result = [(xi, Yi)]

    step = round((xn - x0) / h)
    
    for _ in range(int(step)):
        Fi = F.subs({x: xi})
        for i in range(len(Y0)):
            Fi[i] = Fi[i].subs(dict(zip(Ysym, Yi)))
        Yi = Yi + h * Fi
        xi += h
        result.append((xi, Yi))
    
    return result

if __name__ == '__main__':
    x = sp.symbols('x')
    ###########################################
    Y = sp.Matrix(sp.symbols('y1 y2')) 
    F = sp.Matrix([3 * Y[0] + 2 * Y[1] - (2 * x ** 2 + 1) * sp.exp(2 * x), 
                   4 * Y[0] + Y[1] + (x ** 2 + 2 * x - 4) * sp.exp(2 * x)])
    x0 = 0
    xn = 1
    h = 0.2
    Y0 = sp.Matrix([1, 1])
    F = sp.Matrix([1 / 3 * sp.exp(5 * x) - 1 / 3 * sp.exp(-x) + sp.exp(2 * x),
                  1 / 3 * sp.exp(5 * x) + 2 / 3 * sp.exp(-x) + x ** 2 * sp.exp(2 * x)])
    ########################################
    result = euler_system_method(x0, Y0, h, xn, F, x, Y)
    for xi, Yi in result:
        print(f'x = {xi}, Y = {Yi}')
        print(f"    Error: {abs(Yi - F.subs({x: xi}).evalf())}")
