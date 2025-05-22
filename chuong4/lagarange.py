def lagarange(x: list[float],
              fx: list[float],
              x0: float) -> float:
    n = len(x)
    result = 0
    for i in range(n):
        temp = 1
        for j in range(n):
            if j != i:
                temp *= (x0 - x[j]) / (x[i] - x[j])
        result += fx[i] * temp
    return result

if __name__ == '__main__':
    n = int(input('Nhap bac n cua p(x): '))
    x = []
    fx = []
    for i in range(n + 1):
        x.append(float(input(f'x[{i}]: ')))
        fx.append(float(input(f'f(x[{i}]): ')))
    x0 = float(input('Nhap x0: '))
    print(f'f({x0}) = {lagarange(x, fx, x0)}')