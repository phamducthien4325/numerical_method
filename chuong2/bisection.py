def bisection(f, a, b, tol=1e-5, max_iter=100):
    if f(a) * f(b) > 0:
        print("Khoang a va b khong thoa man")
        return None

    iter_count = 0
    while iter_count < max_iter:
        c = a + (b - a) / 2
        if f(c) == 0:
            return c
        elif (b - a) / 2 < tol:
            return c
        elif f(c) * f(a) < 0:
            b = c
        else:
            a = c
        iter_count += 1
    print("khong tim thay nghiem sau {} vong lap".format(max_iter))
    return None

def func(x):
    return x**3 + 4 * x**2 - 10

if __name__ == "__main__":
    result = bisection(func, 1, 2)
    print(result)
