from gram_schmidt import gram_schmidt
import sympy as sp

def orthogonal(n: int, a: float, b: float, f, phi, x = sp.symbols('x')):
    px = 0

    for i in range(n + 1):
        ai = sp.integrate(f * phi[i], (x, a, b)).evalf() / sp.integrate(phi[i] ** 2, (x, a, b)).evalf()
        px += ai * phi[i]
        
    return sp.simplify(px)

if __name__ == '__main__':
    x = sp.symbols('x')
    f = x ** 2 + 3 * x + 2
    n = 2
    a = 1
    b = 3
    phi = gram_schmidt(n, a, b)
    px = orthogonal(n, a, b, f, phi, x)
    print(f'P(x) = {px}')