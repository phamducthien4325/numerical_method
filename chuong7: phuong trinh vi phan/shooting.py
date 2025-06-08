import sympy as sp
import numpy as np
import math
from runge_kutta_sys import rk4_system_method

def rk_sys(a: float, b: float, h: float, y0: float, yp0: float, Y, y, yp, x)-> list[float]:
    F = sp.Matrix([
        yp,
        Y
    ])
    Y0 = sp.Matrix([y0, yp0])
    result = rk4_system_method(a, Y0, h, b, F, x, sp.Matrix([y, yp]))
    ys = [float(Yi[0].evalf()) for _, Yi in result]
    return ys

def linear_shooting(a: float, b: float, y1_0: float, yp1_0: float, y2_0: float, yp2_0: float, h: float, yn: float, y1, y2, y, yp, x):
    y1_xi = np.array(rk_sys(a, b, h, y1_0, yp1_0, y1, y, yp, x))
    y2_xi = np.array(rk_sys(a, b, h, y2_0, yp2_0, y2, y, yp, x))
    w = y1_xi + ((yn - y1_xi[-1]) / y2_xi[-1]) * y2_xi
    return w.flatten().tolist()

if __name__ == '__main__':
    yp , y = sp.symbols("y' y")
    x = sp.symbols("x")
    ###########################################
    px = 0
    qx = 100
    rx = 0
    a = 0
    b = 1
    y0 = 1
    yn = math.exp(-10)
    h = 0.05
    fx = sp.exp(-10 * x)
    ########################################
    y1 = px * yp + qx * y + rx
    y1_0 = y0
    yp1_0 = 0
    y2 = px * yp + qx * y
    y2_0 = 0
    yp2_0 = 1
    result = linear_shooting(a, b, y1_0, yp1_0, y2_0, yp2_0, h, yn, y1, y2, y, yp, x)
    xi = a
    for w in result:
        print(f'w({xi}) = {w}')
        print(f'    Error: {abs(w - fx.subs({x: xi}).evalf())}')
        xi = round(xi + h, 2)
