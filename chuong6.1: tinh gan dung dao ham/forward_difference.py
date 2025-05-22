#  ____________________________________________
# |          f(x0 + h) - f(x0)     h           |
# |f'(x0) = ------------------- - --- * f''(ξ) |
# |                  h             2           | # với ξ ∈ (x0, x0 + h)
# |____________________________________________|
# neu h < 0 thi goi la backward difference
def forward_difference(f, x0: float, h: float) -> float:
    return (f(x0 + h) - f(x0)) / h
    
def forward_difference_2(fx0: float, fxh: float, x0: float, h: float) -> float:
    return (fxh - fx0) / h

def fx(x):
    return x ** 2 + 3 * x + 2

if __name__ == '__main__':
    x0 = 0.6
    fx0 = 0.5646
    h = -0.1
    fxh = 0.4794
    f_prime = forward_difference_2(fx0, fxh, x0, h)
    print(f"f'({x0}) = {f_prime:.4f}")