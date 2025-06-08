# con goi la Linear Finite-Difference

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

def tiet_kiem(a: float, b: float, alpha: float, beta: float, h: float, px, qx, rx, y, yp, x):
    N = round((b - a) / h - 1)
    xi = a + h
    A = np.zeros((N, N)) 
    rhs = np.zeros(N)      
    a_ = np.zeros(N + 1)
    b_ = np.zeros(N + 1)
    c_ = np.zeros(N + 1)
    d_ = np.zeros(N + 1)
    a_[1] = 2 + h ** 2 * qx.subs({x: xi}).evalf()
    b_[1] = -1 + (h / 2) * px.subs({x: xi}).evalf()
    d_[1] = -h**2 * rx.subs({x: xi}).evalf() + (1 + h / 2 * px.subs({x: xi}).evalf()) * alpha
    for i in range(2, int(N)):
        xi = a + i * h
        a_[i] = 2 + h ** 2 * qx.subs({x: xi}).evalf()
        b_[i] = -1 + (h / 2) * px.subs({x: xi}).evalf()
        c_[i] = -1 - (h / 2) * px.subs({x: xi}).evalf()
        d_[i] = -h**2 * rx.subs({x: xi}).evalf()
    xi = b - h
    a_[N] = 2 + h ** 2 * qx.subs({x: xi}).evalf()
    c_[N] = -1 - (h / 2) * px.subs({x: xi}).evalf()
    d_[N] = -h**2 * rx.subs({x: xi}).evalf() + (1 - h / 2 * px.subs({x: xi}).evalf()) * beta

    l = np.zeros(N + 1)
    u = np.zeros(N + 1)
    z = np.zeros(N + 1)
    l[1] = a_[1]
    u[1] = b_[1] / a_[1]
    z[1] = d_[1] / l[1]

    for i in range(2, N):
        l[i] = a_[i] - c_[i] * u[i - 1]
        u[i] = b_[i] / l[i]
        z[i] = (d_[i] - c_[i] * z[i - 1]) / l[i]
    l[N] = a_[N] - c_[N] * u[N - 1]
    z[N] = (d_[N] - c_[N] * z[N - 1]) / l[N]
    w_full = np.zeros(N + 2)
    w_full[0] = alpha
    w_full[N + 1] = beta
    w_full[N] = z[N]
    for i in range(N - 1, 0, -1):
        w_full[i] = z[i] - u[i] * w_full[i + 1]

    return w_full

if __name__ == '__main__':
    yp , y = sp.symbols("y' y")
    x = sp.symbols("x")
    ###########################################
    px = -4 * x**-1
    qx = -2 * x**-2
    rx = 2 * x**-2 * sp.ln(x)
    a = 1
    b = 2
    y0 = -1/2
    yn = math.log(2)
    h = 0.05
    ########################################
    px, qx, rx = sp.sympify(px), sp.sympify(qx), sp.sympify(rx)
    ypp = px * yp + qx * y + rx
    result = tiet_kiem(a, b, y0, yn, h, px, qx, rx, y, yp, x)
    xi = a
    for w in result:
        print(f'w({xi}) = {w}')
        # print(f'    Error: {abs(w - fx.subs({x: xi}).evalf())}')
        xi = round(xi + h, 2)
