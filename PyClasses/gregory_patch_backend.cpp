#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <cmath>

namespace py = pybind11;

// Bernstein basis functions
double Bernstein(int n, int k, double t) {
    if (k < 0 || k > n) {
        return 0; 
    }
    // Calculate binomial coefficient C(n, i)
    double binom = 1;
    for (int i = 1; i <= k; ++i) {
        binom = binom * (n - i + 1) / i; 
    }
    return binom * std::pow(t, k) * std::pow(1 - t, n - k);
}


// Grg function
Eigen::Vector3d Grg(const Eigen::MatrixXd &CtrlPts, double u, double v, double eps) {
    Eigen::Vector3d p = Eigen::Vector3d::Zero();
    int n = 3; // Degree of Bernstein polynomial in u
    int m = 3; // Degree of Bernstein polynomial in v

    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            Eigen::Vector3d xij = Eigen::Vector3d::Zero(); // Initialize xij
            if (i >= 1 && i <= 2 && j >= 1 && j <= 2) {
                Eigen::Vector3d cp1;
                Eigen::Vector3d cp2;

                // These indices must match the flatCtrlPts order in Python
                if (i == 1 && j == 1) {
                    cp1 = CtrlPts.row(12);
                    cp2 = CtrlPts.row(13);
                } else if (i == 1 && j == 2) {
                    cp1 = CtrlPts.row(18);
                    cp2 = CtrlPts.row(19);
                } else if (i == 2 && j == 1) {
                    cp1 = CtrlPts.row(14);
                    cp2 = CtrlPts.row(15);
                } else { // i == 2 && j == 2
                    cp1 = CtrlPts.row(16);
                    cp2 = CtrlPts.row(17);
                }

                double den;
                if (i == 1 && j == 1) {
                    den = std::max(eps, u + v);
                    xij = (u * cp1 + v * cp2) / den;
                } else if (i == 1 && j == 2) {
                    den = std::max(eps, u + 1 - v);
                    xij = (u * cp1 + (1 - v) * cp2) / den;
                } else if (i == 2 && j == 1) {
                    den = std::max(eps, v + 1 - u);
                    xij = ((1 - u) * cp1 + v * cp2) / den;
                } else { // i == 2 && j == 2
                    den = std::max(eps, 2 - u - v);
                    xij = ((1 - u) * cp1 + (1 - v) * cp2) / den;
                }
            } else {
                // Direct access for boundary and corner nodes
                // These indices must match the flatCtrlPts order in Python
                if (i == 0 && j == 0) xij = CtrlPts.row(0);
                else if (i == 3 && j == 0) xij = CtrlPts.row(1);
                else if (i == 3 && j == 3) xij = CtrlPts.row(2);
                else if (i == 0 && j == 3) xij = CtrlPts.row(3);
                else if (i == 1 && j == 0) xij = CtrlPts.row(4);
                else if (i == 2 && j == 0) xij = CtrlPts.row(5);
                else if (i == 3 && j == 1) xij = CtrlPts.row(6);
                else if (i == 3 && j == 2) xij = CtrlPts.row(7);
                else if (i == 2 && j == 3) xij = CtrlPts.row(8);
                else if (i == 1 && j == 3) xij = CtrlPts.row(9);
                else if (i == 0 && j == 2) xij = CtrlPts.row(10);
                else if (i == 0 && j == 1) xij = CtrlPts.row(11);
            }
            
            double Bi = Bernstein(n, i, u);
            double Bj = Bernstein(m, j, v);
            p += Bi * Bj * xij;
        }
    }
    return p;
}

// Grg_derivs function
py::tuple Grg_derivs(const Eigen::MatrixXd &CtrlPts, double u, double v, double eps) {
    Eigen::Vector3d p = Eigen::Vector3d::Zero();
    Eigen::Vector3d D1p = Eigen::Vector3d::Zero();
    Eigen::Vector3d D2p = Eigen::Vector3d::Zero();
    int n = 3;
    int m = 3;

    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            Eigen::Vector3d xij = Eigen::Vector3d::Zero();
            Eigen::Vector3d D1xij = Eigen::Vector3d::Zero();
            Eigen::Vector3d D2xij = Eigen::Vector3d::Zero();

            if (i >= 1 && i <= 2 && j >= 1 && j <= 2) {
                Eigen::Vector3d cp1;
                Eigen::Vector3d cp2;

                if (i == 1 && j == 1) {
                    cp1 = CtrlPts.row(12);
                    cp2 = CtrlPts.row(13);
                    double den = std::max(eps, u + v);
                    xij = (u * cp1 + v * cp2) / den;
                    D1xij = cp1 / den - (u * cp1 + v * cp2) / (den * den);
                    D2xij = cp2 / den - (u * cp1 + v * cp2) / (den * den);
                } else if (i == 1 && j == 2) {
                    cp1 = CtrlPts.row(18);
                    cp2 = CtrlPts.row(19);
                    double den = std::max(eps, u + 1.0 - v);
                    xij = (u * cp1 + (1.0 - v) * cp2) / den;
                    D1xij = cp1 / den - (u * cp1 + (1.0 - v) * cp2) / (den * den);
                    D2xij = -cp2 / den + (u * cp1 + (1.0 - v) * cp2) / (den * den);
                } else if (i == 2 && j == 1) {
                    cp1 = CtrlPts.row(14);
                    cp2 = CtrlPts.row(15);
                    double den = std::max(eps, v + 1.0 - u);
                    xij = ((1.0 - u) * cp1 + v * cp2) / den;
                    D1xij = -cp1 / den + ((1.0 - u) * cp1 + v * cp2) / (den * den);
                    D2xij = cp2 / den - ((1.0 - u) * cp1 + v * cp2) / (den * den);
                } else { // i == 2 && j == 2
                    cp1 = CtrlPts.row(16);
                    cp2 = CtrlPts.row(17);
                    double den = std::max(eps, 2.0 - u - v);
                    xij = ((1.0 - u) * cp1 + (1.0 - v) * cp2) / den;
                    D1xij = -cp1 / den + ((1.0 - u) * cp1 + (1.0 - v) * cp2) / (den * den);
                    D2xij = -cp2 / den + ((1.0 - u) * cp1 + (1.0 - v) * cp2) / (den * den);
                }
            } else {
                if (i == 0 && j == 0) xij = CtrlPts.row(0);
                else if (i == 3 && j == 0) xij = CtrlPts.row(1);
                else if (i == 3 && j == 3) xij = CtrlPts.row(2);
                else if (i == 0 && j == 3) xij = CtrlPts.row(3);
                else if (i == 1 && j == 0) xij = CtrlPts.row(4);
                else if (i == 2 && j == 0) xij = CtrlPts.row(5);
                else if (i == 3 && j == 1) xij = CtrlPts.row(6);
                else if (i == 3 && j == 2) xij = CtrlPts.row(7);
                else if (i == 2 && j == 3) xij = CtrlPts.row(8);
                else if (i == 1 && j == 3) xij = CtrlPts.row(9);
                else if (i == 0 && j == 2) xij = CtrlPts.row(10);
                else if (i == 0 && j == 1) xij = CtrlPts.row(11);
            }

            double Bi = Bernstein(n, i, u);
            double Bj = Bernstein(m, j, v);
            double D1Bi = n * Bernstein(n - 1, i - 1, u) - n * Bernstein(n - 1, i, u);
            double D2Bj = m * Bernstein(m - 1, j - 1, v) - m * Bernstein(m - 1, j, v);

            p += Bi * Bj * xij;
            D1p += D1Bi * Bj * xij + Bi * Bj * D1xij;
            D2p += Bi * D2Bj * xij + Bi * Bj * D2xij;
        }
    }
    return py::make_tuple(p, D1p, D2p);
}

// Grg_derivs2 function
py::tuple Grg_derivs2(const Eigen::MatrixXd &CtrlPts, double u, double v, double eps) {
    Eigen::Vector3d p = Eigen::Vector3d::Zero();
    Eigen::Vector3d D1p = Eigen::Vector3d::Zero();
    Eigen::Vector3d D2p = Eigen::Vector3d::Zero();
    Eigen::Vector3d D1D1p = Eigen::Vector3d::Zero();
    Eigen::Vector3d D1D2p = Eigen::Vector3d::Zero();
    Eigen::Vector3d D2D2p = Eigen::Vector3d::Zero();
    int n = 3;
    int m = 3;

    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            Eigen::Vector3d xij = Eigen::Vector3d::Zero();
            Eigen::Vector3d D1xij = Eigen::Vector3d::Zero();
            Eigen::Vector3d D2xij = Eigen::Vector3d::Zero();
            Eigen::Vector3d D11xij = Eigen::Vector3d::Zero();
            Eigen::Vector3d D12xij = Eigen::Vector3d::Zero();
            Eigen::Vector3d D22xij = Eigen::Vector3d::Zero();

            if (i >= 1 && i <= 2 && j >= 1 && j <= 2) {
                Eigen::Vector3d cp1;
                Eigen::Vector3d cp2;

                if (i == 1 && j == 1) {
                    cp1 = CtrlPts.row(12);
                    cp2 = CtrlPts.row(13);
                    double den = std::max(eps, u + v);
                    xij = (u * cp1 + v * cp2) / den;
                    D1xij = cp1 / den - (u * cp1 + v * cp2) / (den * den);
                    D2xij = cp2 / den - (u * cp1 + v * cp2) / (den * den);
                    D11xij = -2 * v * (cp1 - cp2) / (den * den * den);
                    D12xij = (u - v) * (cp1 - cp2) / (den * den * den);
                    D22xij = 2 * u * (cp1 - cp2) / (den * den * den);
                } else if (i == 1 && j == 2) {
                    cp1 = CtrlPts.row(18);
                    cp2 = CtrlPts.row(19);
                    double den = std::max(eps, u + 1.0 - v);
                    xij = (u * cp1 + (1.0 - v) * cp2) / den;
                    D1xij = cp1 / den - (u * cp1 + (1.0 - v) * cp2) / (den * den);
                    D2xij = -cp2 / den + (u * cp1 + (1.0 - v) * cp2) / (den * den);
                    D11xij = (2 * (-1 + v) * (cp1 - cp2)) / (den * den * den);
                    D12xij = -((-1 + u + v) * (cp1 - cp2)) / (den * den * den);
                    D22xij = (2 * u * (cp1 - cp2)) / (den * den * den);
                } else if (i == 2 && j == 1) {
                    cp1 = CtrlPts.row(14);
                    cp2 = CtrlPts.row(15);
                    double den = std::max(eps, v + 1.0 - u);
                    xij = ((1.0 - u) * cp1 + v * cp2) / den;
                    D1xij = -cp1 / den + ((1.0 - u) * cp1 + v * cp2) / (den * den);
                    D2xij = cp2 / den - ((1.0 - u) * cp1 + v * cp2) / (den * den);
                    D11xij = (2 * v * (cp1 - cp2)) / (-den * den * den);
                    D12xij = ((-1 + u + v) * (cp1 - cp2)) / (den * den * den);
                    D22xij = (2 * (-1 + u) * (cp1 - cp2)) / (-den * den * den);
                } else { // i == 2 && j == 2
                    cp1 = CtrlPts.row(16);
                    cp2 = CtrlPts.row(17);
                    double den = std::max(eps, 2.0 - u - v);
                    xij = ((1.0 - u) * cp1 + (1.0 - v) * cp2) / den;
                    D1xij = -cp1 / den + ((1.0 - u) * cp1 + (1.0 - v) * cp2) / (den * den);
                    D2xij = -cp2 / den + ((1.0 - u) * cp1 + (1.0 - v) * cp2) / (den * den);
                    D11xij = -((2 * (-1 + v) * (cp1 - cp2)) / (-den * den * den));
                    D12xij = ((u - v) * (cp1 - cp2)) / (-den * den * den);
                    D22xij = (2 * (-1 + u) * (cp1 - cp2)) / (-den * den * den);
                }
            } else {
                if (i == 0 && j == 0) xij = CtrlPts.row(0);
                else if (i == 3 && j == 0) xij = CtrlPts.row(1);
                else if (i == 3 && j == 3) xij = CtrlPts.row(2);
                else if (i == 0 && j == 3) xij = CtrlPts.row(3);
                else if (i == 1 && j == 0) xij = CtrlPts.row(4);
                else if (i == 2 && j == 0) xij = CtrlPts.row(5);
                else if (i == 3 && j == 1) xij = CtrlPts.row(6);
                else if (i == 3 && j == 2) xij = CtrlPts.row(7);
                else if (i == 2 && j == 3) xij = CtrlPts.row(8);
                else if (i == 1 && j == 3) xij = CtrlPts.row(9);
                else if (i == 0 && j == 2) xij = CtrlPts.row(10);
                else if (i == 0 && j == 1) xij = CtrlPts.row(11);
            }

            double Bi = Bernstein(n, i, u);
            double Bj = Bernstein(m, j, v);
            double D1Bi = n * Bernstein(n - 1, i - 1, u) - n * Bernstein(n - 1, i, u);
            double D2Bj = m * Bernstein(m - 1, j - 1, v) - m * Bernstein(m - 1, j, v);
            double DD1Bi = n * (n - 1) * Bernstein(n - 2, i - 2, u) - 2 * n * (n - 1) * Bernstein(n - 2, i - 1, u) + n * (n - 1) * Bernstein(n - 2, i, u);
            double DD2Bj = m * (m - 1) * Bernstein(m - 2, j - 2, v) - 2 * m * (m - 1) * Bernstein(m - 2, j - 1, v) + m * (m - 1) * Bernstein(m - 2, j, v);

            p += Bi * Bj * xij;
            D1p += D1Bi * Bj * xij + Bi * Bj * D1xij;
            D2p += Bi * D2Bj * xij + Bi * Bj * D2xij;
            D1D1p += (DD1Bi * xij + 2 * D1Bi * D1xij + Bi * D11xij) * Bj;
            D1D2p += D1Bi * D2Bj * xij + D1Bi * Bj * D2xij + Bi * D2Bj * D1xij + Bi * Bj * D12xij;
            D2D2p += (DD2Bj * xij + 2 * D2Bj * D2xij + Bj * D22xij) * Bi;
        }
    }
    return py::make_tuple(p, D1p, D2p, D1D1p, D1D2p, D2D2p);
}

// MinDist function
py::tuple MinDist(const Eigen::MatrixXd &CtrlPts, const Eigen::Vector3d &x, int seeding, double eps) {
    double umin = 0.0;
    double vmin = 0.0;
    Eigen::Vector3d initial_point = CtrlPts.row(0);
    double dmin = (x - initial_point).norm();

    for (int i = 0; i <= seeding; ++i) {
        double u = static_cast<double>(i) / seeding;
        for (int j = 0; j <= seeding; ++j) {
            double v = static_cast<double>(j) / seeding;
            Eigen::Vector3d p = Grg(CtrlPts, u, v, eps);
            double d = (x - p).norm();
            if (d < dmin) {
                dmin = d;
                umin = u;
                vmin = v;
            }
        }
    }
    return py::make_tuple(umin, vmin);
}

// Helper function for D3Grg
Eigen::Vector3d D3Grg_helper(const Eigen::MatrixXd &CtrlPts, double u, double v, double eps, bool normalize) {
    py::tuple derivs = Grg_derivs(CtrlPts, u, v, eps);
    Eigen::Vector3d D1p = derivs[1].cast<Eigen::Vector3d>();
    Eigen::Vector3d D2p = derivs[2].cast<Eigen::Vector3d>();
    Eigen::Vector3d D3p = D1p.cross(D2p);
    if (normalize) {
        double norm_D3p = D3p.norm();
        if (norm_D3p > 1e-12) { // Avoid division by zero
            D3p /= norm_D3p;
        }
    }
    return D3p;
}

// find_projection function
py::tuple find_projection(const Eigen::MatrixXd &CtrlPts, const Eigen::Vector3d &xs, py::tuple t_py, double bs_r, double eps) {
    Eigen::Vector2d t(t_py[0].cast<double>(), t_py[1].cast<double>());

    double tol = 1e-15;
    double res = 1.0 + tol;
    int niter = 0;
    Eigen::Vector2d tcandidate = t;

    Eigen::Vector3d xc_candidate = Grg(CtrlPts, tcandidate.x(), tcandidate.y(), eps);
    double dist = (xs - xc_candidate).norm();

    double opa = 1e-2;

    while (res > tol && (t.x() >= -opa && t.x() <= 1.0 + opa) && (t.y() >= -opa && t.y() <= 1.0 + opa)) {
        py::tuple derivs2 = Grg_derivs2(CtrlPts, t.x(), t.y(), eps);
        Eigen::Vector3d xc = derivs2[0].cast<Eigen::Vector3d>();
        Eigen::Vector3d D1p = derivs2[1].cast<Eigen::Vector3d>();
        Eigen::Vector3d D2p = derivs2[2].cast<Eigen::Vector3d>();
        Eigen::Vector3d D1D1p = derivs2[3].cast<Eigen::Vector3d>();
        Eigen::Vector3d D1D2p = derivs2[4].cast<Eigen::Vector3d>();
        Eigen::Vector3d D2D2p = derivs2[5].cast<Eigen::Vector3d>();

        Eigen::Matrix<double, 3, 2> dxcdt;
        dxcdt.col(0) = D1p;
        dxcdt.col(1) = D2p;

        Eigen::Vector2d f = -2 * dxcdt.transpose() * (xs - xc);

        Eigen::Matrix2d K;
        K(0,0) = 2.0 * (-(xs - xc).dot(D1D1p) + D1p.dot(D1p));
        K(0,1) = 2.0 * (-(xs - xc).dot(D1D2p) + D1p.dot(D2p));
        K(1,0) = K(0,1);
        K(1,1) = 2.0 * (-(xs - xc).dot(D2D2p) + D2p.dot(D2p));

        Eigen::Vector2d dt = -K.inverse() * f;
        t += dt;

        res = dt.norm();

        if (res < std::sqrt(tol) && !(t.x() > 0 && t.x() < 1 && t.y() > 0 && t.y() < 1)) {
            return py::make_tuple(-1.0, -1.0);
        }

        niter++;
        if (niter > 10) {
            Eigen::Vector3d xc_new = Grg(CtrlPts, t.x(), t.y(), eps);
            double dist_new = (xs - xc_new).norm();
            if (dist_new < dist) {
                dist = dist_new;
                tcandidate = t;
            }
            if (niter > 13) {
                // proj_final_check for tcandidate
                if (!(tcandidate.x() > 0 && tcandidate.x() < 1 && tcandidate.y() > 0 && tcandidate.y() < 1)) {
                    double t1 = std::min(std::max(0.0, tcandidate.x()), 1.0);
                    double t2 = std::min(std::max(0.0, tcandidate.y()), 1.0);
                    Eigen::Vector3d xc0 = Grg(CtrlPts, t1, t2, eps);
                    Eigen::Vector3d nor0 = D3Grg_helper(CtrlPts, t1, t2, eps, true);
                    Eigen::Vector3d x_tang = (xs - xc0) - (xs - xc0).dot(nor0) * nor0;
                    if (x_tang.norm() > 2.0 * bs_r / 100.0) {
                        return py::make_tuple(-1.0, -1.0);
                    }
                }
                return py::make_tuple(tcandidate.x(), tcandidate.y());
            }
        }
    }

    // proj_final_check for t
    if (!(t.x() > 0 && t.x() < 1 && t.y() > 0 && t.y() < 1)) {
        double t1 = std::min(std::max(0.0, t.x()), 1.0);
        double t2 = std::min(std::max(0.0, t.y()), 1.0);
        Eigen::Vector3d xc0 = Grg(CtrlPts, t1, t2, eps);
        Eigen::Vector3d nor0 = D3Grg_helper(CtrlPts, t1, t2, eps, true);
        Eigen::Vector3d x_tang = (xs - xc0) - (xs - xc0).dot(nor0) * nor0;
        if (x_tang.norm() > 2.0 * bs_r / 100.0) {
            return py::make_tuple(-1.0, -1.0);
        }
    }

    return py::make_tuple(t.x(), t.y());
}


PYBIND11_MODULE(gregory_patch_backend, m) {
    m.doc() = "C++ backend for Gregory patch calculations";
    m.def("Grg", &Grg, "A function that calculates a point on a Gregory patch");
    m.def("Grg_derivs", &Grg_derivs, "A function that calculates first derivatives of Grg");
    m.def("Grg_derivs2", &Grg_derivs2, "A function that calculates second derivatives of Grg");
    m.def("MinDist", &MinDist, "A function that calculates the minimum distance from a point to a Gregory patch");
    m.def("find_projection", &find_projection, "A function that finds the projection of a point onto a Gregory patch");
}