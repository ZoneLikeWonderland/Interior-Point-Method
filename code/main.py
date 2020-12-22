import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from dataset import *
import time
import os
import glob
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
        AAT = A.dot(A.T)
        x_0 = A.T.dot(np.linalg.pinv(AAT).dot(b))
        lam_0 = np.linalg.pinv(AAT).dot(A.dot(c))
        s_0 = c-A.T.dot(lam_0)

        d_x, d_s = max(-np.min(x_0), 0), max(-np.min(s_0), 0)
        x_0 = x_0+d_x
        s_0 = s_0+d_s

        u = 0.1
        dh_x, dh_s = max(x_0.T.dot(s_0)/(np.sum(s_0)+u)*u, u), max(x_0.T.dot(s_0)/np.sum(x_0+u)*u, u)
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
    status = None

    rank = np.linalg.matrix_rank(A)
    if rank < A.shape[0]:
        print("WARNING: not full rank", rank, "<", A.shape[0])

    m, n = A.shape
    x_k = np.random.random((n, 1))
    lam_k = np.random.random((m, 1))
    s_k = np.random.random((n, 1))
    x_k, lam_k, s_k = init_iterate()

    eps = 1e-9
    gamma = 0.1
    c1 = 0.9

    time_cost["init"] += time.time()-start

    F_k = F(x_k, lam_k, s_k)
    if detail:
        print(
            "k=  -1",
            "f(x)={:5.2f},".format(f_primal(x_k).item()),
            "|F(x)_d|= {:5.6e},".format(np.linalg.norm(F_k[0:n], np.inf)),
            "|F(x)_p|= {:5.6e},".format(np.linalg.norm(F_k[n:n+m], np.inf)),
            "|F(x)_0|= {:5.6e},".format(np.linalg.norm(F_k[n+m:n+m+n], np.inf)),
        )
        old_dz_k = np.zeros((n+m+n, 1))

    for k in range(500):
        start = time.time()
        gF_k = gF(x_k, lam_k, s_k)
        F_k = F(x_k, lam_k, s_k)

        time_cost["gF&F"] += time.time()-start
        start = time.time()

        mu = x_k.T.dot(s_k)
        sigma = 0.01
        toler = mu*sigma/n
        # if False:
        #     rhs=-F_k
        #     rhs[n+m:n+m+n] += toler

        #     dz_k, resdual, rank, singular=np.linalg.lstsq(gF_k.T, rhs)
        #     dx_k=dz_k[0:n]
        #     dlam_k=dz_k[n:n+m]
        #     ds_k=dz_k[n+m:n+m+n]

        F_d = F_k[0:n]
        F_p = F_k[n:n+m]
        F_0 = F_k[n+m:n+m+n]-toler

        XS_inv = np.diag((x_k/s_k).flat)
        S_invF_0 = F_0/s_k
        AXS_inv = A.dot(XS_inv)
        dlam_k, _, _, _ = np.linalg.lstsq(AXS_inv.dot(A.T), -F_p-AXS_inv.dot(F_d)+A.dot(S_invF_0), rcond=None)
        ds_k = -F_d-A.T.dot(dlam_k)
        dx_k = -S_invF_0-XS_inv.dot(ds_k)
        dz_k = np.vstack((dx_k, dlam_k, ds_k))
        time_cost["solve"] += time.time()-start
        start = time.time()

        alpha = 1

        alpha = min(min(np.min((-x_k/dx_k)[np.where(dx_k < 0)]), np.min((-s_k/ds_k)[np.where(ds_k < 0)]))*0.9, 1)

        alpha_t = alpha
        while np.linalg.norm(F(x_k+alpha*dx_k, lam_k+alpha*dlam_k, s_k+alpha*ds_k), np.inf) >= np.linalg.norm(F_k+c1*alpha*gF_k.T.dot(dz_k), np.inf):
            alpha *= gamma
            if alpha*gamma == 0:
                print("WARNING: alpha=0")
                break
            if (x_k+alpha*dx_k).T.dot(s_k+alpha*ds_k) < gamma*mu:
                print("WARNING: too small")
                break

        time_cost["line"] += time.time()-start
        start = time.time()

        x_k += alpha*dx_k
        lam_k += alpha*dlam_k
        s_k += alpha*ds_k

        time_cost["update"] += time.time()-start
        start = time.time()

        if detail:
            print(
                "k={:3d},".format(k),
                "f(x)= {:5.3f},".format(f_primal(x_k).item()),
                "|dz_d|= {:5.6e},".format(np.linalg.norm(dz_k[0:n], np.inf)),
                "|dz_p|= {:5.6e},".format(np.linalg.norm(dz_k[n:n+m], np.inf)),
                "|dz_0|= {:5.6e},".format(np.linalg.norm(dz_k[n+m:n+m+n], np.inf)),
                "|F(x)_d|= {:5.6e},".format(np.linalg.norm(F_k[0:n], np.inf)),
                "|F(x)_p|= {:5.6e},".format(np.linalg.norm(F_k[n:n+m], np.inf)),
                "|F(x)_0|= {:5.6e},".format(np.linalg.norm(F_k[n+m:n+m+n], np.inf)),
                "alpha= {:.5e},".format(alpha),
                "theta= {:.5e}".format(
                    # np.linalg.norm(dz_k-old_dz_k)
                    np.linalg.norm(alpha*dz_k)
                )
            )
            old_dz_k = dz_k

        if np.linalg.norm(F_k, np.inf) < eps:
            status = "DONE"
            break

        if np.linalg.norm(dz_k, np.inf) == 0:
            raise Exception("FAILED to solve: update delta=0")
        if np.linalg.norm(dz_k, np.inf) > 1e31:
            raise Exception("FAILED to solve: update delta exploded")
        if np.linalg.norm(alpha*dz_k, np.inf) < eps:
            print("WARNING: predicted cannot converge")
            status = "HARD"
            break

    if status is None:
        status = "EXPIRED"
    x_k[np.where(np.abs(x_k) < eps)] = 0
    total_time = sum([i for _, i in time_cost.items()])
    if detail:
        print("x =", x_k.T.tolist())
        print(time_cost)
        print("total time", total_time)

    return x_k, f_primal(x_k).item(), np.linalg.norm(F_k, np.inf), total_time, status


if __name__ == "__main__":

    # TIMES=1
    # t=0
    # for i in range(TIMES):
    #     x, tc=IPM(A, b, c, detail=i == 0)
    #     t += tc
    # print("avg time", t/TIMES)

    # A, b, c, x_star = read_data(r"data/LP_MATLAB\beaconfd.mat")
    # x, tc, status = IPM(A, b, c, detail=True)
    # exit()

    for path in glob.glob("data/LP_MATLAB/*.mat"):
        print("try", path)
        A, b, c, x_star = read_data(path)
        try:
            x, tgt, res, tc, status = IPM(A, b, c, detail=False)
            print("\t"*4, status, tgt, res)
        except Exception as e:
            print(e)
