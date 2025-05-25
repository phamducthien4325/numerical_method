import sympy as sp

def euler_method(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> list[tuple[float, float]]:
    """
    Euler's method for solving ordinary differential equations.
    
    :param x0: initial x value
    :param y0: initial y value
    :param h: step size
    :param xn: final x value
    :param y_prime: the derivative of y with respect to x (dy/dx)
    :param x: sympy symbol for x
    :param y: sympy symbol for y
    :return: list of tuples (x, y) representing the solution
    """
    xi = x0
    yi = y0
    result = [(xi, yi)]
    
    while xi < xn:
        yi += h * y_prime.subs({x: xi, y: yi}).evalf()
        xi += h
        result.append((xi, yi))
    
    return result

if __name__ == '__main__':
    x = sp.symbols('x')
    y = sp.symbols('y')
    y_prime = 1 + (x - y) ** 2
    x0 = 2
    y0 = 1
    h = 0.5
    xn = 3
    solution = euler_method(x0, y0, h, xn, y_prime, x, y)
    print("Solution using Euler's method:")
    for x_val, y_val in solution:
        print(f"x: {x_val}, y: {y_val}")
