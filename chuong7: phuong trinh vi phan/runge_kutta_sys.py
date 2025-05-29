import sympy as sp

def rk4_system_method(x0: float, Y0: sp.Matrix, h: float, xn: float, F, x, Ysym) -> list[tuple[float, float]]:
    xi = x0
    Yi = Y0
    result = [(xi, Yi)]

    step = round((xn - x0) / h)
    
    for _ in range(int(step)):
        subs1 = dict(zip(Ysym, Yi))
        k1 = F.subs({x: xi}).applyfunc(lambda f: f.subs(subs1))

        subs2 = dict(zip(Ysym, Yi + h/2 * k1))
        k2 = F.subs({x: xi + h/2}).applyfunc(lambda f: f.subs(subs2))

        subs3 = dict(zip(Ysym, Yi + h/2 * k2))
        k3 = F.subs({x: xi + h/2}).applyfunc(lambda f: f.subs(subs3))

        subs4 = dict(zip(Ysym, Yi + h * k3))
        k4 = F.subs({x: xi + h}).applyfunc(lambda f: f.subs(subs4))
        Yi = Yi + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        xi += h
        result.append((xi, Yi))
    
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
    result = rk4_system_method(x0, Y0, h, xn, F, x, Y)
    for xi, Yi in result:
        print(f'x = {xi}, Y = {Yi}')
        print(f"    Error: {abs(Yi - actual.subs({x: xi}).evalf())}")
