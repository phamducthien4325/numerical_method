def newton_inter(x: list[float],
                 fx: list[float],
                 x0: float) -> float:
    n = len(x)
    result = fx[0]
    mul = 1
    for i in range(1, n):
        mul *= (x0 - x[i - 1])
        for j in range(n - i):
            fx[j] = (fx[j + 1] - fx[j]) / (x[j + i] - x[j])
        result += mul * fx[0]
    return result

if __name__ == '__main__':
    n = int(input('Nhap bac n cua p(x): '))
    x = []
    fx = []
    for i in range(n + 1):
        x.append(float(input(f'x[{i}]: ')))
        fx.append(float(input(f'f(x[{i}]): ')))
    x0 = float(input('Nhap x0: '))
    print(f'f({x0}) = {newton_inter(x, fx, x0)}')