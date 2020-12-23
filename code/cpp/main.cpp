#include "Eigen/Core"
#include "Eigen/Dense"
#include <bits/stdc++.h>
using namespace std;
using namespace Eigen;

#define TIMER_START(name) auto start_time_##name = chrono::steady_clock::now();
#define TIMER_RECORD(name)                                                     \
    time_cost[#name] +=                                                        \
        (chrono::steady_clock::now() - start_time_##name).count() / 1e9;

template <class DataType> class IPM {
    typedef Matrix<DataType, -1, -1> Mat;
    typedef Matrix<DataType, -1, 1> Vec;
    Mat A;
    Vec b;
    Vec c;
    size_t m, n;

    DataType eps = 1e-6;
    DataType terminal = 1e-5;
    DataType gamma = 0.9;
    DataType c1 = 1e-4;

    Mat gF_k;
    Vec F_k;
    Vec x_k, lam_k, s_k;
    Vec dx_k, dlam_k, ds_k;
    Vec F_d, F_p, F_0;
    DataType alpha_p, alpha_d;
    enum Status { TBD, DONE, HARD };

  public:
    auto read_file(string path) {
        ifstream file(path);
        string row_str;
        vector<DataType> data;
        size_t rows = 0;
        while (file >> row_str) {
            istringstream row(row_str);
            DataType i;
            char ch;
            while (row >> i) {
                data.emplace_back(i);
                row >> ch;
            }
            rows += 1;
        }
        size_t cols = data.size() / rows;
        Mat map_matrix =
            Map<Matrix<DataType, -1, -1, RowMajor>>(data.data(), rows, cols);
        return map_matrix;
    }

    DataType f_primal(Vec x) { return (x.transpose() * c).value(); }

    auto F(Vec x, Vec lam, Vec s) {
        Vec ret(n + m + n, 1);
        ret << A.transpose() * lam + s - c, A * x - b, x.array() * s.array();
        return ret;
    }

    auto gF(Vec x, Vec lam, Vec s) {
        Mat g = Mat::Zero(n + m + n, n + m + n);
        g.block(0, n, n, m) = A.transpose();
        g.block(n, 0, m, n) = A;
        g.block(n + m, 0, n, n) = Mat::Identity(n, n);
        g.block(0, n + m, n, n) = DiagonalMatrix<DataType, -1>(s);
        g.block(n + m, n + m, n, n) = DiagonalMatrix<DataType, -1>(x);
        return g;
    }

    tuple<Vec, Vec, Vec> init_iterate() {
        Mat AATinv = (A * A.transpose())
                         .completeOrthogonalDecomposition()
                         .pseudoInverse();

        Vec x_0 = A.transpose() * (AATinv * b);
        Vec lam_0 = AATinv * (A * c);
        Vec s_0 = c - A.transpose() * lam_0;
        // cout << x_0;

        auto d_x = max(-(double)x_0.minCoeff(), 0.);
        auto d_s = max(-(double)s_0.minCoeff(), 0.);
        x_0.array() += d_x;
        s_0.array() += d_s;

        DataType u = 1;
        auto x_dot_s = (x_0.transpose() * s_0).value();
        auto dh_x = x_dot_s / (s_0.sum() + eps) * u + eps;
        auto dh_s = x_dot_s / (x_0.sum() + eps) * u + eps;
        x_0.array() += dh_x;
        s_0.array() += dh_s;

        if ((x_0.array() != x_0.array()).any() ||
            ((x_0 - x_0).array() != (x_0 - x_0).array()).any()) {
            cout << "x_0 nan\n";
            exit(-1);
        }
        if ((s_0.array() != s_0.array()).any() ||
            ((s_0 - s_0).array() != (s_0 - s_0).array()).any()) {
            cout << "s_0 nan\n";
            exit(-1);
        }
        // cout << s_0 << "\n";

        return {x_0, lam_0, s_0};
    }

    auto set_Abc(Mat A, Vec b, Vec c) {
        this->A = A;
        this->b = b;
        this->c = c;
        m = A.rows();
        n = A.cols();
    }

    tuple<Vec, Vec, Vec> solve_newton(Vec F_d, Vec F_p, Vec F_0) {
        auto XS_inv =
            DiagonalMatrix<DataType, -1>((x_k.array() / s_k.array()).matrix());
        // cout << x_k.array() / s_k.array() << "\n";
        // cout << XS_inv << "\n";
        auto S_invF_0 = (F_0.array() / s_k.array()).matrix();
        auto AXS_inv = A * XS_inv;
        // dlam_k, _, _,
        //     _ = np.linalg.lstsq(AXS_inv.dot(A.T),
        //                         -F_p - AXS_inv.dot(F_d) + A.dot(S_invF_0),
        //                         rcond = -1);
        Vec dlam_k = (AXS_inv * A.transpose())
                         .colPivHouseholderQr()
                         .solve(-F_p - AXS_inv * F_d + A * S_invF_0);
        Vec ds_k = -F_d - A.transpose() * dlam_k;
        Vec dx_k = -S_invF_0 - XS_inv * ds_k;
        // auto dz_k = np.vstack((dx_k, dlam_k, ds_k));
        return {dx_k, dlam_k, ds_k};
        // return dx_k, dlam_k, ds_k, dz_k;
        // return {};
    }

    auto maximum_alpha(Vec _k, Vec d_k) {
        auto select = d_k.array() < 0;
        _k = select.select(_k, 1);
        d_k = select.select(-d_k, 0);
        DataType alpha =
            min((DataType)((_k.array() / d_k.array()).minCoeff()), 1.) *
            (1 - eps);
        if (!isnormal(alpha)) {
            cout << "maximum_alpha fail\n";
            exit(-1);
        }
        return alpha;
    }

    auto solve(bool detail = false) {
        auto total_start_time = chrono::steady_clock::now();
        unordered_map<string, double> time_cost;
        TIMER_START(init);

        Status status = TBD;
        // auto r = f_primal(c);
        tie(x_k, lam_k, s_k) = init_iterate();
        TIMER_RECORD(init);

        // cout << x_k << "\n";
        // cout << lam_k << "\n";
        // cout << s_k << "\n";

        // if (detail) {
        //     F_k = F(x_k, lam_k, s_k);
        //     F_d = F_k.block(0, 0, n, 1);
        //     F_p = F_k.block(n, 0, m, 1);
        //     F_0 = F_k.block(n + m, 0, n, 1);
        //     printf("k= -1"
        //            "f(x)=%5.2f,"
        //            "|F(x)_d|=%.3e,"
        //            "|F(x)_p|=%.3e,"
        //            "|F(x)_0|=%.3e,",
        //            f_primal(x_k), F_d.template lpNorm<Infinity>(),
        //            F_p.template lpNorm<Infinity>(),
        //            F_0.template lpNorm<Infinity>());
        // }
        // cout << F_k << "\n";

        for (int k = 1; k <= 500; k++) {
            TIMER_START(calF);
            gF_k = gF(x_k, lam_k, s_k);
            F_k = F(x_k, lam_k, s_k);
            TIMER_RECORD(calF);

            // cout << F_k.rows() << " " << n << " " << m << "\n";

            // time_cost["gF&F"] += time.time()-start
            // start = time.time()
            // cout << F_k.rows() << " " << F_k.cols() << "\n";
            TIMER_START(solve);
            F_d = F_k.block(0, 0, n, 1);
            F_p = F_k.block(n, 0, m, 1);
            F_0 = F_k.block(n + m, 0, n, 1);
            // cout << F_k << "\n";

            tie(dx_k, dlam_k, ds_k) = solve_newton(F_d, F_p, F_0);

            alpha_p = maximum_alpha(x_k, dx_k);
            alpha_d = maximum_alpha(s_k, ds_k);
            // cout << alpha_p << " " << alpha_d << "\n";

            auto mu = (x_k.transpose() * s_k).value() / n;
            auto sigma =
                ((x_k + alpha_p * dx_k).transpose() * (s_k + alpha_d * ds_k))
                    .value() /
                n / mu;
            auto toler = mu * sigma / n;

            F_d = F_k.block(0, 0, n, 1);
            F_p = F_k.block(n, 0, m, 1);
            F_0 = F_k.block(n + m, 0, n, 1).array() - toler;
            tie(dx_k, dlam_k, ds_k) = solve_newton(F_d, F_p, F_0);

            // DataType yita = 0.9;
            alpha_p = maximum_alpha(x_k, dx_k);
            alpha_d = maximum_alpha(s_k, ds_k);
            TIMER_RECORD(solve);
            // cout << alpha_p << " " << alpha_d << "\n";

            TIMER_START(search);

            while ((DataType)(((x_k + alpha_p * dx_k).array() *
                               (s_k + alpha_d * ds_k).array())
                                  .maxCoeff()) <= gamma * mu) {
                alpha_p *= gamma;
                alpha_d *= gamma;
                if (alpha_d < eps * eps and alpha_p < eps * eps) {
                    cout << "WARNING: alpha too small\n";
                    break;
                }
                if (alpha_p * gamma == 0 or alpha_d * gamma == 0) {
                    cout << "WARNING: alpha=0\n";
                    break;
                }
            }

            Vec dz_k(n + m + n, 1);
            dz_k << alpha_p * dx_k, alpha_d * dlam_k, alpha_d * ds_k;
            Vec desc = c1 * gF_k.transpose() * dz_k;

            while (F(x_k + alpha_p * dx_k, lam_k + alpha_d * dlam_k,
                     s_k + alpha_d * ds_k)
                       .template lpNorm<Infinity>() >=
                   (F_k + desc).template lpNorm<Infinity>()) {
                alpha_p *= gamma;
                alpha_d *= gamma;
                desc *= gamma;
                if (alpha_d < eps * eps and alpha_p < eps * eps) {
                    cout << "WARNING: alpha too small\n";
                    break;
                }
                if (alpha_p * gamma == 0 or alpha_d * gamma == 0) {
                    cout << "WARNING: alpha=0\n";
                    break;
                }
            }

            TIMER_RECORD(search);

            if (detail) {
                printf("k=%3d,"
                       "f(x)=%5.6f,"
                       "|dx|=%.3e,"
                       "|dlam|=%.3e,"
                       "|ds|=%.3e,"
                       "|F(x)_d|=%.3e,"
                       "|F(x)_p|=%.3e,"
                       "|F(x)_0|=%.3e,"
                       "alpha_p=%.3e,"
                       "alpha_d=%.3e,"
                    //    "dot_x=%.5e"
                       "\n",
                       k, f_primal(x_k), dx_k.template lpNorm<Infinity>(),
                       dlam_k.template lpNorm<Infinity>(),
                       ds_k.template lpNorm<Infinity>(),
                       F_d.template lpNorm<Infinity>(),
                       F_p.template lpNorm<Infinity>(),
                       F_0.template lpNorm<Infinity>(), alpha_p, alpha_d);
                // cout << x_k << "\n";
            }

            if (F_k.template lpNorm<Infinity>() < eps) {
                status = DONE;
                break;
            }

            // if np.linalg.norm(dz_k, 1) == 0:
            //     raise Exception("FAILED to solve: update delta=0")
            if (dz_k.template lpNorm<Infinity>() > 1e31) {
                cout << "FAILED to solve: update delta exploded\n";
                exit(-1);
            }
            if ((alpha_p * dx_k).template lpNorm<Infinity>() < eps and
                (alpha_d * ds_k).template lpNorm<Infinity>() < eps) {
                if (F_k.template lpNorm<Infinity>() < terminal) {
                    status = DONE;
                    break;
                }
                cout << "WARNING: predicted cannot converge\n";
                status = HARD;
                break;
            }

            TIMER_START(update);
            x_k += alpha_p * dx_k;
            lam_k += alpha_d * dlam_k;
            s_k += alpha_d * ds_k;
            TIMER_RECORD(update);
        }

        auto total_wall_time_cost =
            (chrono::steady_clock::now() - total_start_time).count();
        cout << "total_wall_time_cost = " << total_wall_time_cost / 1e9 << "\n";
        for (auto item : time_cost) {
            cout << item.first << ":" << item.second << "\n";
        }
    }
};

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "usage: " << argv[0] << " <data_folder>\n";
        exit(-1);
    }
    string folder = argv[1];
    if (folder[folder.size() - 1] == '/')
        folder = folder.substr(0, folder.size() - 1);

    IPM<double> ipm;
    auto A = ipm.read_file(folder + "/A.csv");
    auto b = ipm.read_file(folder + "/b.csv");
    auto c = ipm.read_file(folder + "/c.csv");

    ipm.set_Abc(A, b, c);
    ipm.solve(true);

    cout << "executable finished\n";
    return 0;
}