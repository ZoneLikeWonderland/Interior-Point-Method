import numpy as np
from enum import Enum
import scipy.sparse as ss
import scipy.linalg as sl
import os
import time
import argparse


class IPM:
    class Mat(np.ndarray):
        '''
        mxn
        '''
    class Vec(np.ndarray):
        '''
        mx1
        '''

    max_iterations = 500
    eps = 1e-9
    terminal = 1e-6
    gamma = 0.9
    c1 = 1e-4

    use_M = False

    class Status(Enum):
        TBD = 0
        DONE = 1
        HARD = 2
        EXPIRED = 3
        UNBOUNDED = 4
        VIOLATED = 5

        @staticmethod
        def spec(status):
            return [
                "TBD",
                "successfully solved",
                "hard to solve (cannot converge)",
                "max iterations exceeded",
                "feasible but unbounded",
                "infeasible"
            ][status.value]

    class SD_type(Enum):
        LLT = 0
        LS = 1

    def read_file(self, path: str, required: bool = True) -> np.ndarray:
        if not os.path.exists(path):
            if not required:
                return None
            print("[ ERROR ] could not find", path)
            exit(-1)
        mat = np.loadtxt(path, delimiter=",", ndmin=2)
        return mat

    def f_primal(self, x: Vec) -> float:
        return (x.T@self.c).item()

    def F(self, x: Vec, lam: Vec, s: Vec) -> Vec:
        return np.vstack((
            self.A.T@lam+s-self.c,
            self.A@x - self.b,
            x*s
        ))

    def gF(self, x: Vec, lam: Vec, s: Vec) -> None:
        m, n = self.m, self.n
        for i in range(n):
            self.gF_k[i, n+m+i] = s[i, 0]
            self.gF_k[n+m+i, n+m+i] = x[i, 0]

    def init_iterate(self) -> (Vec, Vec, Vec):
        AATinv = np.linalg.pinv(self.A_dense@self.A_dense.T)
        x_0 = self.A.T@(AATinv@self.b)
        lam_0 = AATinv@(self.A@self.c)
        s_0 = self.c - self.A.T@lam_0
        x_0 += max(-np.min(x_0), 0)
        s_0 += max(-np.min(s_0), 0)

        x_dot_s = (x_0.T@s_0).item()
        dh_x = x_dot_s / (s_0.sum()+self.eps)+self.eps
        dh_s = x_dot_s / (x_0.sum()+self.eps)+self.eps
        x_0 += dh_x
        s_0 += dh_s
        return x_0, lam_0, s_0

    def enforce_fullrank(self, A: Mat, b: Vec) -> (Mat, Vec):
        Araw = A.copy()
        AT = A.T
        nrows, ncols = AT.shape
        ind_index = []
        leadi = 0
        leadj = 0
        while leadi < nrows and leadj < ncols:
            max_index = leadi
            for r in range(leadi+1, nrows):
                if (abs(AT[r, leadj]) > abs(AT[max_index, leadj])):
                    max_index = r
            if (AT[max_index, leadj] == 0):
                leadj += 1
                continue

            AT[[leadi, max_index], :] = AT[[max_index, leadi], :]
            ind_index.append(leadj)
            AT[leadi] /= AT[leadi, leadj]
            for r in range(leadi+1, nrows):
                d = AT[leadi, leadj]
                m = AT[r, leadj] / d
                AT[r] -= AT[leadi] * m
            leadi += 1
            leadj += 1
        A_t = np.zeros((len(ind_index), A.shape[1]))
        b_t = np.zeros((len(ind_index), b.shape[1]))
        for i in range(len(ind_index)):
            A_t[i] = Araw[ind_index[i]]
            b_t[i] = b[ind_index[i]]
        return A_t, b_t

    def init(self, A_r: Mat, b_r: Vec, c: Vec, use_M: bool) -> None:
        A, b = self.enforce_fullrank(A_r, b_r)
        m, n = A.shape
        if use_M:
            M = max(np.max(np.abs(c))*np.mean(np.abs(b))*m, m)

            A = np.hstack((A, np.identity(m)))
            c = np.vstack((c, M*np.ones((m, 1))))
            m, n = A.shape
            self.use_M = True

        self.m, self.n = m, n
        self.A = ss.csr_matrix(A)
        self.A_dense = A
        self.b = b
        self.c = c

        self.gF_k = np.zeros((n+m+n, n+m+n))
        self.gF_k[0:n, n:n+m] = A.T
        self.gF_k[n:n+m, 0:n] = A
        self.gF_k[n+m:n+m+n, 0:n] = np.identity(n)
        self.gF_k = ss.lil_matrix(self.gF_k)

    def solve_newton(self, F_d: Vec, F_p: Vec,
                     F_0: Vec, first: bool) -> (Vec, Vec, Vec):
        if (first):
            self.XS_inv = np.diag((self.x_k/self.s_k).flat)
            AXS_inv = self.A@self.XS_inv
            self.AXS_invAT = AXS_inv@self.A.T
            self.rhs_0 = -F_p - AXS_inv@F_d

            try:
                self.L = sl.cholesky(self.AXS_invAT, lower=True)
            except:
                self.L = None

        S_invF_0 = F_0/self.s_k
        if self.L is not None:
            dlam_k = sl.cho_solve((self.L, True), self.rhs_0+self.A@S_invF_0)
        else:
            dlam_k = np.linalg.lstsq(self.AXS_invAT, self.rhs_0+self.A@S_invF_0, rcond=-1)[0]

        ds_k = -F_d - self.A.T @ dlam_k
        dx_k = -S_invF_0 - self.XS_inv @ ds_k
        return dx_k, dlam_k, ds_k

    def maximum_alpha(self, _k: Vec, d_k: Vec) -> float:
        _k = _k.copy()
        d_k = d_k.copy()
        select = np.where(d_k >= 0)
        _k[select] = np.inf
        d_k[select] = -1
        alpha = min(np.min(_k / -d_k), 1) * (1 - self.eps)
        return alpha

    def solve(self, args, is_presolve: bool = False, optimal: Mat = None) -> Status:
        total_start_time = time.time()
        status = self.Status.TBD
        m, n = self.m, self.n
        stay = 0
        self.x_k, self.lam_k, self.s_k = self.init_iterate()

        for k in range(1, self.max_iterations+1):
            F_k = self.F(self.x_k, self.lam_k, self.s_k)
            self.gF(self.x_k, self.lam_k, self.s_k)

            F_d = F_k[0:n]
            F_p = F_k[n:n+m]
            F_0 = F_k[n+m:n+m+n]

            dx_k, dlam_k, ds_k = self.solve_newton(F_d, F_p, F_0, True)

            alpha_p = self.maximum_alpha(self.x_k, dx_k)
            alpha_d = self.maximum_alpha(self.s_k, ds_k)

            mu = (self.x_k.T@self.s_k).item() / n
            sigma = (((self.x_k+alpha_p*dx_k).T@(self.s_k+alpha_d*ds_k)).item() / n / mu)**2
            toler = mu*sigma / n

            F_d = F_k[0:n]
            F_p = F_k[n:n+m]
            F_0 = F_k[n+m:n+m+n] - toler
            dx_k, dlam_k, ds_k = self.solve_newton(F_d, F_p, F_0, False)

            alpha_p = self.maximum_alpha(self.x_k, dx_k)
            alpha_d = self.maximum_alpha(self.s_k, ds_k)

            while (np.max((self.x_k+alpha_p*dx_k) *
                          (self.s_k+alpha_d*ds_k))) <= self.gamma * mu:
                alpha_p *= self.gamma
                alpha_d *= self.gamma
                if (alpha_d < self.eps**2 and alpha_p < self.eps**2):
                    break

            dz_k = np.vstack((alpha_p * dx_k, alpha_d * dlam_k, alpha_d * ds_k))
            desc = self.c1 * self.gF_k.T * dz_k

            while (np.linalg.norm(self.F(self.x_k+alpha_p * dx_k, self.lam_k+alpha_d * dlam_k,
                                         self.s_k+alpha_d * ds_k), np.inf) >=
                   np.linalg.norm((F_k+desc), np.inf)):
                alpha_p *= self.gamma
                alpha_d *= self.gamma
                desc *= self.gamma
                if (alpha_d < self.eps**2 and alpha_p < self.eps**2):
                    break

            if (args is not None and args.detail):
                print("k=%3d,"
                      "f(x)=%5.6f,"
                      "|dx|=%.3e,"
                      "|dlam|=%.3e,"
                      "|ds|=%.3e,"
                      "|F(x)_d|=%.3e,"
                      "|F(x)_p|=%.3e,"
                      "|F(x)_0|=%.3e,"
                      "alpha_p=%.3e,"
                      "alpha_d=%.3e,"
                      "" %
                      (k, self.f_primal(self.x_k), np.linalg.norm(dx_k, np.inf),
                       np.linalg.norm(dlam_k, np.inf),
                       np.linalg.norm(ds_k, np.inf),
                       np.linalg.norm(F_d, np.inf),
                       np.linalg.norm(F_p, np.inf),
                       np.linalg.norm(F_0, np.inf), alpha_p, alpha_d))

            if (np.linalg.norm(F_k, np.inf) < self.terminal):
                status = self.Status.DONE
                break

            if (alpha_p < self.eps and alpha_d < self.eps):
                stay += 1
                if stay > 3:
                    if (np.linalg.norm(F_k, np.inf) < 1e-3):
                        status = self.Status.HARD
                    elif (np.linalg.norm(dz_k, np.inf) > 1e3):
                        status = self.Status.UNBOUNDED
                    break

            self.x_k += alpha_p * dx_k
            self.lam_k += alpha_d * dlam_k
            self.s_k += alpha_d * ds_k

        if (is_presolve):
            return status
        else:
            if (status.value > self.Status.HARD.value):
                c_real = self.c
                self.c = np.zeros_like(self.c)
                presolve_status = self.solve(None, True)
                self.c = c_real
                if (presolve_status.value > self.Status.HARD.value):
                    status = self.Status.VIOLATED

        if self.use_M:
            if (status == self.Status.TBD or status == self.Status.DONE) and np.linalg.norm(self.x_k[-m:], np.inf) > self.terminal:
                status = self.Status.VIOLATED
        if (status == self.Status.TBD):
            status = self.Status.EXPIRED

        total_wall_time_cost = time.time()-total_start_time

        print("total solving wall time =", total_wall_time_cost, "sec")
        print("status                  =", status.name, ":", self.Status.spec(status))
        if (status == self.Status.DONE):
            print("primal-dual optimal solution:")
            print("optimal objective value = %.4f" % self.f_primal(self.x_k))
            print("numbers of iterations   =", k)
            sol = self.x_k
            if self.use_M:
                sol = sol[:-m]
            if args is not None and args.solution:
                print("solution x              = {")
                print(",\n".join([", ".join(["%10.4f" % j for j in sol.T[0][i:i+10]])
                                  for i in range(0, sol.shape[0], 10)]),
                      "\n}")
            if optimal is not None:
                print("given x*:||x_k-x*||_inf =", np.linalg.norm(sol-optimal, np.inf),
                      "(optimal solution may not be unique)")
        return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data_folder", help="a folder containing `A.csv`, `b.csv` and `c.csv` [and `x_star.csv`]")
    parser.add_argument("--detail", required=False, default=False,
                        help="show detailed information each iteration", action="store_true")
    parser.add_argument("--solution", required=False, default=False,
                        help="show solution if solved", action="store_true")
    parser.add_argument("-M", required=False, default=False,
                        help="use the big-M form", action="store_true")
    args = parser.parse_args()

    ipm = IPM()
    A0 = ipm.read_file(os.path.join(args.data_folder, "A.csv"))
    b0 = ipm.read_file(os.path.join(args.data_folder, "b.csv"))
    c0 = ipm.read_file(os.path.join(args.data_folder, "c.csv"))
    x0 = ipm.read_file(os.path.join(args.data_folder, "x_star.csv"), required=False)

    ipm.init(A0, b0, c0, use_M=args.M)
    ipm.solve(args, optimal=x0)

    print("solver terminated successfully")
