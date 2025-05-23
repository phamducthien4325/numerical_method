#  Cong thuc Composite Trapezoidal Rule, khi n = 1 thi cong thuc tro thanh Trapezoidal Rule
#  ____________________________________________________________________________________
# |                      h                 n-1               (b - a) * h^2            | ξ ∈ (a, b)
# |∫[a to b] f(x) dx  = --- * [f(a) + 2 *   Σ f(xj)+ f(b)] - --------------- * f''(ξ) | h = (b - a) / n
# |                      2                 j=1                     12                 | xj = a + j * h
# |___________________________________________________________________________________| j = 0, 1, 2, ..., n
import sympy as sp

def trapezoidal(f, a: float, b: float, n: int) -> float:
    """
    :param f: function to be integrated
    :param a: lower limit of integration
    :param b: upper limit of integration
    :param n: number of subintervals
    :return: approximate value of the integral
    """
    h = (b - a) / n
    return h / 2 * (f(a) + 2 * sum(f(a + j * h) for j in range(1, n)) + f(b))

if __name__ == "__main__":
    x = sp.symbols('x')
    f = x ** 4
    a = 0.5
    b = 1
    n = 1
    h = (b - a) / n
    integ = trapezoidal(sp.lambdify(x, f), a, b, n)
    print(f"Trapezoidal Rule: {integ}")
    print(f"Error: {abs(sp.integrate(f, (x, a, b)) - integ)}")
