import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from dataset import *
import time
style.use("ggplot")


def IPM(A, b, c, detail=False):

    def f_primal(x):
        return c.T.dot(x)

    def F(x, lam, s):
        return np.vstack((
            A.T.dot(lam)+s-c,
            A.dot(x)-b,
            np.multiply(x, s)
        ))

    def gF(x, lam, s):
        g = np.zeros((n*2+m, n*2+m))
        g[0:n, n:n+m] = A.T
        g[n:n+m, 0:n] = A
        g[n+m:n+m+n, 0:n] = np.identity(n)
        g[0:n, n+m:n+m+n] = np.diag(s.flat)
        g[n+m:n+m+n, n+m:n+m+n] = np.diag(x.flat)
        return g

    def init_iterate():
        x_0 = A.T.dot(np.linalg.inv(A.dot(A.T)).dot(b))
        lam_0 = np.linalg.inv(A.dot(A.T)).dot(A.dot(c))
        s_0 = c-A.T.dot(lam_0)

        d_x, d_s = max(-np.min(x_0), 0), max(-np.min(s_0), 0)
        x_0 = x_0+d_x
        s_0 = s_0+d_s

        u = 1
        dh_x, dh_s = x_0.T.dot(s_0)/np.sum(s_0)*u, x_0.T.dot(s_0)/np.sum(x_0)*u
        x_0 = x_0+dh_x
        s_0 = s_0+dh_s

        return x_0, lam_0, s_0

    time_cost = {
        "init": 0,
        "gF&F": 0,
        "solve": 0,
        "line": 0,
        "update": 0
    }
    start = time.time()

    m, n = A.shape
    x_k = np.random.random((n, 1))
    lam_k = np.random.random((m, 1))
    s_k = np.random.random((n, 1))
    x_k, lam_k, s_k = init_iterate()

    eps = 1e-9
    gamma = 0.9
    c1 = 1e-4

    time_cost["init"] += time.time()-start

    F_k = F(x_k, lam_k, s_k)
    if detail:
        print(
            "k=-1",
            "f(x) = {:5.2f},".format(f_primal(x_k).item()),
            "|F(x)_d| = {:5.10f},".format(np.linalg.norm(F_k[0:n], np.inf)),
            "|F(x)_p| = {:5.10f},".format(np.linalg.norm(F_k[n:n+m], np.inf)),
            "|F(x)_0| = {:5.10f},".format(np.linalg.norm(F_k[n+m:n+m+n], np.inf)),
        )
    for k in range(1000):
        start = time.time()
        # gF_k = gF(x_k, lam_k, s_k)
        F_k = F(x_k, lam_k, s_k)

        time_cost["gF&F"] += time.time()-start
        start = time.time()

        mu = x_k.T.dot(s_k)
        sigma = 0.01
        toler = mu*sigma/n
        if False:
            rhs = -F_k
            rhs[n+m:n+m+n] += toler

            dz_k, resdual, rank, singular = np.linalg.lstsq(gF_k.T, rhs)
            dx_k = dz_k[0:n]
            dlam_k = dz_k[n:n+m]
            ds_k = dz_k[n+m:n+m+n]

        F_d = F_k[0:n]
        F_p = F_k[n:n+m]
        F_0 = F_k[n+m:n+m+n]-toler

        XS_inv = np.diag(np.multiply(x_k, 1/s_k).flat)
        S_invF_0 = F_0/s_k
        AXS_inv = A.dot(XS_inv)
        dlam_k = np.linalg.solve(AXS_inv.dot(A.T), -F_p-AXS_inv.dot(F_d)+A.dot(S_invF_0))
        ds_k = -F_d-A.T.dot(dlam_k)
        dx_k = -S_invF_0-XS_inv.dot(ds_k)
        dz_k = np.vstack((dx_k, dlam_k, ds_k))
        time_cost["solve"] += time.time()-start
        start = time.time()

        alpha = 1

        while (x_k+alpha*dx_k < 0).any() or (s_k+alpha*ds_k < 0).any():
            alpha *= gamma

        # while (F(x_k+alpha*dx_k, lam_k, s_k+alpha*ds_k < 0) >= F_k+c1*alpha*dz_k.T.dot(dz_k)).any():
        #     alpha *= gamma

        time_cost["line"] += time.time()-start
        start = time.time()

        x_k += alpha*dx_k
        lam_k += alpha*dlam_k
        s_k += alpha*ds_k

        time_cost["update"] += time.time()-start
        start = time.time()

        if detail:
            print(
                "k={:4d},".format(k),
                "f(x) = {:5.2f},".format(f_primal(x_k).item()),
                "|dz_d| = {:5.10f},".format(np.linalg.norm(dz_k[0:n], np.inf)),
                "|dz_p| = {:5.10f},".format(np.linalg.norm(dz_k[n:n+m], np.inf)),
                "|dz_0| = {:5.10f},".format(np.linalg.norm(dz_k[n+m:n+m+n], np.inf)),
                "|F(x)_d| = {:5.10f},".format(np.linalg.norm(F_k[0:n], np.inf)),
                "|F(x)_p| = {:5.10f},".format(np.linalg.norm(F_k[n:n+m], np.inf)),
                "|F(x)_0| = {:5.10f},".format(np.linalg.norm(F_k[n+m:n+m+n], np.inf)),
                "alpha = {:f},".format(alpha)
            )

        if np.linalg.norm(F_k, np.inf) < eps:
            break

    total_time = sum([i for _, i in time_cost.items()])
    if detail:
        print(time_cost)
        print("total time", total_time)

    return x_k, total_time


if __name__ == "__main__":

    A, b, c, x_star = read_data("./data/data1")
    A = A[:-1]
    b = b[:-1]

    t = 0
    for i in range(100):
        x, tc = IPM(A, b, c, detail=i == 0)
        t += tc
    print("avg time", t/100)
