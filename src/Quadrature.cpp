// Quadrature.cpp
#include "Quadrature.hpp"
#include <stdexcept>
#include <algorithm>
#include "Logger.hpp"

namespace Quadrature {

double legendre(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;

    return ((2.0 * n - 1.0) * x * legendre(n - 1, x) - (n - 1.0) * legendre(n - 2, x)) / n;
}

double dLegendre(int n, double x) {
    if (n == 0) return 0.0;
    return n * (x * legendre(n, x) - legendre(n - 1, x)) / (x * x - 1.0);
}

std::vector<double> legendreRoots(int n) {
    std::vector<double> roots(n);
    double tol = 1e-10;
    int maxIter = 100;
    for (int i = 0; i < n; ++i) {
        double x0 = std::cos(M_PI * (i + 0.75) / (n + 0.5));
        double x = x0;
        double dx = 0.0;
        for (int iter = 0; iter < maxIter; ++iter) {
            double Pn = legendre(n, x);
            double dPn = dLegendre(n, x);
            dx = -Pn / dPn;
            x += dx;
            if (std::abs(dx) < tol) break;
        }
        double Pn_final = legendre(n, x);
        if (std::abs(Pn_final) > tol) {
            std::cerr << "Warning: Gauss-Legendre root " << i+1 << " did not converge.\n";
        }
        roots[i] = x;
    }
    std::sort(roots.begin(), roots.end());
    return roots;
}

std::vector<double> legendreWeights(const std::vector<double>& roots, int n) {
    std::vector<double> weights(roots.size());

    for (size_t i = 0; i < roots.size(); ++i) {
        double x = roots[i];
        double dPn = dLegendre(n, x);
        weights[i] = 2.0 / ((1.0 - x * x) * dPn * dPn);
    }

    return weights;
}

std::pair<std::vector<double>, std::vector<double>> gaussChebyshev(int n) {
    std::vector<double> points(n);
    std::vector<double> weights(n, 2.0 * M_PI / n); // Poids constants pour Gauss-Chebyshev

    for (int i = 1; i <= n; ++i) {
        points[i-1] = M_PI * (2.0 * i - 1.0) / (n);
    }

    return {points, weights};
}

} // namespace Quadrature
