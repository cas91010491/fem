#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
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

PYBIND11_MODULE(gregory_patch_backend, m) {
    m.doc() = "C++ backend for Gregory patch calculations";
    m.def("Grg", &Grg, "A function that calculates a point on a Gregory patch");
}