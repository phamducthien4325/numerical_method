from gram_schmidt import gram_schmidt
import sympy as sp

def orthogonal(n: int, a: float, b: float, f, phi, w = 1, x = sp.symbols('x')):
    px = 0

    for i in range(n + 1):
        ai = sp.integrate(w * phi[i] * f, (x, a, b)).evalf() / sp.integrate(w * phi[i] ** 2, (x, a, b)).evalf()
        px += ai * phi[i]
        
    return sp.simplify(px)

if __name__ == '__main__':
    x = sp.symbols('x')
    f = sp.exp(-2 * x)
    w = sp.exp(-x)
    a = 0
    b = sp.oo
    n = 1
    for n in range(1, 4):
        phi = gram_schmidt(n, a, b, w, x)
        px = orthogonal(n, a, b, f, phi, w, x)
        print(f'P{n}(x) = {px}')