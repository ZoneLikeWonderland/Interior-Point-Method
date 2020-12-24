#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <bits/stdc++.h>
using namespace std;
using namespace Eigen;

#define TIMER_START(name) auto start_time_##name = chrono::steady_clock::now();
#define TIMER_RECORD(name)                                                     \
    time_cost[#name] +=                                                        \
        (chrono::steady_clock::now() - start_time_##name).count() / 1e9;

template <class DataType> class IPM {
    typedef Matrix<DataType, -1, -1> Mat;
    typedef SparseMatrix<DataType> SMat;
    typedef Matrix<DataType, -1, 1> Vec;
    SMat A;
    Mat A_dense;
    Vec b;
    Vec c;
    size_t m, n;

    size_t max_iterations = 500;
    DataType eps = 1e-6;
    DataType terminal = 1e-5;
    DataType gamma = 0.9;
    DataType c1 = 1e-4;

    DiagonalMatrix<DataType, -1> XS_inv;
    SMat AXS_inv;
    SMat AXS_invAT;
    Vec rhs_0;

    Mat gF_k;
    Vec F_k;
    Vec x_k, lam_k, s_k;
    Vec dx_k, dlam_k, ds_k;
    Vec F_d, F_p, F_0;
    DataType alpha_p, alpha_d;
    enum Status { TBD, DONE, HARD, EXPIRED, UNBOUNDED, VIOLATED };
    string status_str[6] = {"TBD",     "DONE",      "HARD",
                            "EXPIRED", "UNBOUNDED", "VIOLATED"};

    SimplicialLLT<SMat> solver_LLT;
    enum SD_type { LLT, LS };
    SD_type solver_type;

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

    DataType f_primal(const Vec &x) { return (x.transpose() * c).value(); }

    auto F(const Vec &x, const Vec &lam, const Vec &s) {
        Vec ret(n + m + n, 1);
        ret << A.transpose() * lam + s - c, A * x - b, x.array() * s.array();
        return ret;
    }

    auto gF(const Vec &x, const Vec &lam, const Vec &s) {
        for (int i = 0; i < n; i++) {
            gF_k(i, n + m + i) = s(i, 0);
            gF_k(n + m + i, n + m + i) = x(i, 0);
        }
    }

    tuple<Vec, Vec, Vec> init_iterate() {
        Mat AATinv = (A_dense * A_dense.transpose())
                         .completeOrthogonalDecomposition()
                         .pseudoInverse();

        Vec x_0 = A.transpose() * (AATinv * b);
        Vec lam_0 = AATinv * (A * c);
        Vec s_0 = c - A.transpose() * lam_0;

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
        return {x_0, lam_0, s_0};
    }

    tuple<Mat, Vec> enforce_fullrank(const Mat &A, const Vec &b) {
        Mat AT = A.transpose();
        const int nrows = AT.rows();
        const int ncols = AT.cols();
        vector<int> ind_index;
        for (int leadi = 0, leadj = 0; leadi < nrows and leadj < ncols;
             leadi++, leadj++) {
            auto max_index = leadi;
            for (int r = leadi + 1; r < nrows; r++) {
                if (abs(AT(r, leadj)) > abs(AT(max_index, leadj)))
                    max_index = r;
            }
            if (AT(max_index, leadj) == 0) {
                leadi--;
                continue;
            }
            AT.row(leadi).swap(AT.row(max_index));
            ind_index.emplace_back(leadj);
            AT.row(leadi) /= AT(leadi, leadj);
            for (int r = leadi + 1; r < nrows; r++) {
                auto d = AT(leadi, leadj);
                auto m = AT(r, leadj) / d;
                AT.row(r) -= AT.row(leadi) * m;
            }
        }
        Mat A_t(ind_index.size(), A.cols());
        Vec b_t(ind_index.size(), b.cols());
        for (int i = 0; i < ind_index.size(); i++) {
            A_t.row(i) = A.row(ind_index[i]);
            b_t(i) = b(ind_index[i]);
        }
        return {A_t, b_t};
    }

    auto init(const Mat &A_r, const Vec &b_r, const Vec &c) {
        Mat A;
        Vec b;
        tie(A, b) = enforce_fullrank(A_r, b_r);

        this->A = A.sparseView();
        this->A_dense = A;
        this->b = b;
        this->c = c;
        m = A.rows();
        n = A.cols();

        gF_k = Mat::Zero(n + m + n, n + m + n);
        gF_k.block(0, n, n, m) = A.transpose();
        gF_k.block(n, 0, m, n) = A;
        gF_k.block(n + m, 0, n, n) = Mat::Identity(n, n);
    }

    tuple<Vec, Vec, Vec> solve_newton(const Vec &F_d, const Vec &F_p,
                                      const Vec &F_0, bool first) {
        if (first) {
            XS_inv = DiagonalMatrix<DataType, -1>(
                (x_k.array() / s_k.array()).matrix());
            AXS_inv = A * XS_inv;
            AXS_invAT = AXS_inv * A.transpose();
            rhs_0 = -F_p - AXS_inv * F_d;

            solver_LLT.compute(AXS_invAT);
            if (solver_LLT.info() == Success) {
                solver_type = LLT;
            } else {
                solver_type = LS;
            }
        }

        auto S_invF_0 = (F_0.array() / s_k.array()).matrix();

        Vec dlam_k;
        if (solver_type == LLT)
            dlam_k = solver_LLT.solve(rhs_0 + A * S_invF_0);
        else if (solver_type == LS)
            dlam_k = Mat(AXS_invAT).colPivHouseholderQr().solve(rhs_0 +
                                                                A * S_invF_0);

        Vec ds_k = -F_d - A.transpose() * dlam_k;
        Vec dx_k = -S_invF_0 - XS_inv * ds_k;
        return {dx_k, dlam_k, ds_k};
    }

    auto maximum_alpha(const Vec &_k, const Vec &d_k) {
        auto select = d_k.array() < 0;
        DataType alpha = min((DataType)((select.select(_k, 1).array() /
                                         select.select(-d_k, 0).array())
                                            .minCoeff()),
                             1.) *
                         (1 - eps);
        if (!isnormal(alpha)) {
            cout << "maximum_alpha fail\n";
            exit(-1);
        }
        return alpha;
    }

    Status solve(bool detail = false, bool is_presolve = false) {
        auto total_start_time = chrono::steady_clock::now();
        unordered_map<string, double> time_cost;
        Status status = TBD;

        TIMER_START(presolve);
        if (!is_presolve) {
            auto c_real = c;
            c = Mat::Zero(c.rows(), c.cols());
            auto presolve_status = solve(false, true);
            c = c_real;
            // cout << "presolve_status = " << status_str[presolve_status] << "\n";
            if (presolve_status > HARD) { status = VIOLATED; }
        }
        TIMER_RECORD(presolve);

        if (status == TBD) {

            TIMER_START(init);
            tie(x_k, lam_k, s_k) = init_iterate();
            TIMER_RECORD(init);

            for (int k = 1; k <= max_iterations; k++) {
                TIMER_START(calF);
                F_k = F(x_k, lam_k, s_k);
                TIMER_RECORD(calF);
                TIMER_START(calgF);
                gF(x_k, lam_k, s_k);
                TIMER_RECORD(calgF);

                TIMER_START(solve);
                F_d = F_k.block(0, 0, n, 1);
                F_p = F_k.block(n, 0, m, 1);
                F_0 = F_k.block(n + m, 0, n, 1);

                tie(dx_k, dlam_k, ds_k) = solve_newton(F_d, F_p, F_0, true);

                alpha_p = maximum_alpha(x_k, dx_k);
                alpha_d = maximum_alpha(s_k, ds_k);

                auto mu = (x_k.transpose() * s_k).value() / n;
                auto sigma = ((x_k + alpha_p * dx_k).transpose() *
                              (s_k + alpha_d * ds_k))
                                 .value() /
                             n / mu;
                auto toler = mu * sigma / n;

                F_d = F_k.block(0, 0, n, 1);
                F_p = F_k.block(n, 0, m, 1);
                F_0 = F_k.block(n + m, 0, n, 1).array() - toler;
                tie(dx_k, dlam_k, ds_k) = solve_newton(F_d, F_p, F_0, false);

                alpha_p = maximum_alpha(x_k, dx_k);
                alpha_d = maximum_alpha(s_k, ds_k);
                TIMER_RECORD(solve);

                TIMER_START(search);

                while ((DataType)(((x_k + alpha_p * dx_k).array() *
                                   (s_k + alpha_d * ds_k).array())
                                      .maxCoeff()) <= gamma * mu) {
                    alpha_p *= gamma;
                    alpha_d *= gamma;
                    if (alpha_d < eps * eps and alpha_p < eps * eps) { break; }
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
                    if (alpha_d < eps * eps and alpha_p < eps * eps) { break; }
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
                           "\n",
                           k, f_primal(x_k), dx_k.template lpNorm<Infinity>(),
                           dlam_k.template lpNorm<Infinity>(),
                           ds_k.template lpNorm<Infinity>(),
                           F_d.template lpNorm<Infinity>(),
                           F_p.template lpNorm<Infinity>(),
                           F_0.template lpNorm<Infinity>(), alpha_p, alpha_d);
                }

                if (F_k.template lpNorm<Infinity>() < terminal) {
                    status = DONE;
                    break;
                }

                if (dz_k.template lpNorm<Infinity>() > 1e31) {
                    status = UNBOUNDED;
                    break;
                }
                if ((alpha_p * dx_k).template lpNorm<Infinity>() < eps and
                    (alpha_d * ds_k).template lpNorm<Infinity>() < eps) {
                    if (F_k.template lpNorm<Infinity>() < sqrt(terminal))
                        status = HARD;
                    else
                        status = UNBOUNDED;
                    break;
                }

                TIMER_START(update);
                x_k += alpha_p * dx_k;
                lam_k += alpha_d * dlam_k;
                s_k += alpha_d * ds_k;
                TIMER_RECORD(update);
            }

            if (status == TBD) { status = EXPIRED; }
            if (is_presolve) return status;
        }

        auto total_wall_time_cost =
            (chrono::steady_clock::now() - total_start_time).count();
        cout << "total_wall_time_cost = " << total_wall_time_cost / 1e9 << "\n";
        double total_time_cost = 0;
        for (auto item : time_cost) {
            cout << item.first << ":" << item.second << "\n";
            total_time_cost += item.second;
        }
        cout << "total_time_cost = " << total_time_cost << "\n";
        cout << "status = " << status_str[status] << " tgt = " << f_primal(x_k)
             << " |F(x_k)| = " << F_k.template lpNorm<Infinity>() << "\n";
        return status;
    }
};

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "usage: " << argv[0] << " <data_folder>\n";
        exit(-1);
    }
    string folder = argv[1];
    if (folder[folder.size() - 1] == '/')
        folder = folder.substr(0, folder.size() - 1);

    auto detail = false;
    if (argc > 2 and string(argv[2]) == "--detail") detail = true;

    IPM<double> ipm;
    auto A = ipm.read_file(folder + "/A.csv");
    auto b = ipm.read_file(folder + "/b.csv");
    auto c = ipm.read_file(folder + "/c.csv");

    ipm.init(A, b, c);
    ipm.solve(detail);

    cout << "executable finished\n";
    return 0;
}