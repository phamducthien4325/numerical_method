import sympy as sp

def gram_schmidt(n: int, a: float, b: float, w = 1, x = sp.symbols('x')) -> list:
    phi = []
    phi0 = 1
    phi.append(phi0)
    if n == 0:
        return phi
    b1 = sp.integrate(x * w, (x, a, b)).evalf() / sp.integrate(w, (x, a, b)).evalf()
    phi1 = x - b1
    phi.append(phi1)
    if n == 1:
        return phi
    for i in range(2, n + 1):
        bi = sp.integrate(x * w * phi[i - 1] ** 2, (x, a, b)).evalf() / sp.integrate(w * phi[i - 1] ** 2, (x, a, b)).evalf()
        ci = sp.integrate(x * w * phi[i - 1] * phi[i - 2], (x, a, b)).evalf() / sp.integrate(w * phi[i - 2] ** 2, (x, a, b)).evalf()
        phi_i = (x - bi) * phi[i - 1] - ci * phi[i - 2]
        phi.append(sp.simplify(phi_i))
    return phi


if __name__ == '__main__':
    n = 1
    a = 0
    b = sp.oo
    x = sp.symbols('x')
    w = sp.exp(-x)
    phi = gram_schmidt(n, a, b, w, x)
    for i in range(n + 1):
        if i == 0:
            print(f'φ_{i}(x) = {phi[i]}')
        else:
            print(f'φ_{i}(x) = {phi[i]}')