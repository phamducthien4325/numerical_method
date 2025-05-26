import sympy as sp
import math

def mid_point(x0: float, y0: float, h: float, xn: float, y_prime, x, y) -> list[tuple[float, float]]:
    xi = x0
    yi = y0
    result = [(xi, yi)]
    
    while xi < xn:
        y1_2 = yi + h / 2 * y_prime.subs({x: xi, y: yi}).evalf()
        x1_2 = xi + h / 2
        yi += h * y_prime.subs({x: x1_2, y: y1_2})
        xi += h
        result.append((xi, yi))
    
    return result

if __name__ == '__main__':
    x = sp.symbols('x')
    y = sp.symbols('y')
    y_prime = x * sp.exp(3 * x) - 2 * y
    x0 = 0
    xn = 1
    y0 = 0
    h = 0.5
    f = 1 / 5 * x * sp.exp(3 * x) - 1/25 * sp.exp(3 * x) + 1 / 25 * sp.exp(-2 * x)
    solution = mid_point(x0, y0, h, xn, y_prime, x, y)
    print("Solution using Modified Euler's method:")
    for x_val, y_val in solution:
        print(f"x: {x_val}, y: {y_val}", end='; ')
        print(f"Error: {abs(y_val - f.subs(x, x_val).evalf())}")
    # Min, Max = error_bound(f, x, x0, xn)
    # Min = abs(h**2 / 2 * Min)
    # Max = abs(h**2 / 2 * Max)
    # print(f"Error Bound: [{Min}, {Max}]")
