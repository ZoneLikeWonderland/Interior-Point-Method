import numpy as np
import os
import scipy.io


def read_data(data_folder_or_path):
    if data_folder_or_path.endswith(".mat"):
        f = scipy.io.loadmat(data_folder_or_path)
        try:
            A = f["Problem"]["A"][0][0].todense()
            b = f["Problem"]["b"][0][0]
            c = f["Problem"]["c"][0][0]
        except:
            A = f["A"].todense()
            b = f["b"]
            c = f["c"]
        x_star = None
    else:
        A, b, c, x_star = [
            np.matrix(np.loadtxt(os.path.join(data_folder_or_path, f"{name}.csv"), delimiter=","))
            for name in ["A", "b", "c", "x_star"]
        ]

    m, n = A.shape
    b = b.reshape((m, 1))
    c = c.reshape((n, 1))
    if x_star is not None:
        x_star = x_star.reshape((n, 1))

    return A, b, c, x_star


if __name__ == "__main__":
    A, b, c, x = read_data("data/LP_MATLAB/afiro.mat")
    # print(A)
    # print(b)
    # print(c)
    # print(x)
