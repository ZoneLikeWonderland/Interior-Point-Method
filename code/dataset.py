import numpy as np
import os


def read_data(data_folder):
    A, b, c, x_star = [
        np.matrix(np.loadtxt(os.path.join(data_folder, f"{name}.csv"), delimiter=","))
        for name in ["A", "b", "c", "x_star"]
    ]

    m, n = A.shape
    b = b.reshape((m, 1))
    c = c.reshape((n, 1))
    x_star = x_star.reshape((n, 1))

    return A, b, c, x_star


if __name__ == "__main__":
    A, b, c, x = read_data("./data/data1")
    print(A)
    print(b)
    print(c)
    print(x)
