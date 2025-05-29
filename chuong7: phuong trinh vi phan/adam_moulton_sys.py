import sympy as sp
from runge_kutta_sys import rk4_system_method

def adam_moulton_system_method(x0: float, Y0: sp.Matrix, h: float, xn: float, F, x, Ysym) -> sp.Matrix:
    init = rk4_system_method(x0, Y0, h, x0 + 3 * h, F, x, Ysym)
    result = init[:]
    xi = x0 + 3 * h
    Yi = [init[0][1], init[1][1], init[2][1], init[3][1]]

    step = round((xn - x0) / h - 3)
    
    for _ in range(int(step)):
        subs = [dict(zip(Ysym, Yi[0])), dict(zip(Ysym, Yi[1])), dict(zip(Ysym, Yi[2])), dict(zip(Ysym, Yi[3]))]
        Y_pred = Yi[3] + h / 24 * (55 * F.subs({x: xi}).applyfunc(lambda f: f.subs(subs[3])) -
                59 * F.subs({x: xi - h}).applyfunc(lambda f: f.subs(subs[2])) + 
                37 * F.subs({x: xi - 2 * h}).applyfunc(lambda f: f.subs(subs[1])) - 
                9 * F.subs({x: xi - 3 * h}).applyfunc(lambda f: f.subs(subs[0])))
        Y_corr = Yi[3] + h / 24 * (9 * F.subs({x: xi + h}).applyfunc(lambda f: f.subs(dict(zip(Ysym, Y_pred)))) +
                19 * F.subs({x: xi}).applyfunc(lambda f: f.subs(subs[3])) -
                5 * F.subs({x: xi - h}).applyfunc(lambda f: f.subs(subs[2])) +
                F.subs({x: xi - 2 * h}).applyfunc(lambda f: f.subs(subs[1])))
        Yi[0] = Yi[1]
        Yi[1] = Yi[2]
        Yi[2] = Yi[3]
        Yi[3] = Y_corr
        xi += h
        result.append((xi, Y_corr))
    
    return result

if __name__ == '__main__':
    x = sp.symbols('x')
    ###########################################
    x = sp.symbols('x')
    Y = sp.Matrix(sp.symbols('y1 y2'))
    F = sp.Matrix([
        Y[0] - Y[1] + 2,
        -Y[0] + Y[1] + 4 * x
    ])
    Y0 = sp.Matrix([-1, 0])
    x0 = 0
    xn = 1
    h = 0.1
    actual = sp.Matrix([
        -1 / 2 * sp.exp(2 * x) + x**2 + 2 * x - 1 / 2,
        0.5 * sp.exp(2 * x) + x**2 - 0.5
    ])
    ########################################
    result = adam_moulton_system_method(x0, Y0, h, xn, F, x, Y)
    for xi, Yi in result:
        print(f'x = {xi}, Y = {Yi}')
        print(f"    Error: {abs(Yi - actual.subs({x: xi}).evalf())}")
