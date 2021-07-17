from tabulate import tabulate


def benchmark_parameters(vert, hor, vert_name, hor_name, evalf, tablefmt="grid", num_experiments = 1):
    data = [[None for j in range(len(hor) + 1)] for i in range(len(vert) + 1)]

    for i in range(len(vert)):
        data[i + 1][0] = str(vert[i])

    for j in range(len(hor)):
        data[0][j + 1] = str(hor[j])

    data[0][0] = f"{vert_name} / {hor_name}"

    for _ in range(num_experiments):
        for i in range(1, len(vert) + 1):
            for j in range(1, len(hor) + 1):
                if data[i][j] is None:
                    data[i][j] = evalf(vert[i - 1], hor[j - 1])
                else:
                    data[i][j] = [sum(x) for x in zip(data[i][j], evalf(vert[i - 1], hor[j - 1]))]

    for i in range(1, len(vert) + 1):
        for j in range(1, len(hor) + 1):
            data[i][j] = ', '.join([str(elem) for elem in [round(x / num_experiments, 3) for x in data[i][j]]])

    return tabulate(data, tablefmt=tablefmt)


def compare_methods (rows, tablefmt="grid"):
    header = ['Method', 'RealTime', 'Correct', 'Wrong', 'Unmapped']
    data = []
    data.append(header)
    for row in rows:
        data.append(row)
    return tabulate(data, tablefmt=tablefmt)


def compare_real_methods (rows, tablefmt="grid"):
    header = ['Method', 'RealTime', 'Unmapped']
    data = []
    data.append(header)
    for row in rows:
        data.append(row)
    return tabulate(data, tablefmt=tablefmt)