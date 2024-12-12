// AngularQuadrature.cpp
#include "AngularQuadrature.hpp"
#include "Quadrature.hpp"
#include "Logger.hpp"
#include <cmath>

namespace SNSolver {

AngularQuadrature::AngularQuadrature(int theta_order, int phi_order) : thetaOrder_(theta_order), phiOrder_(phi_order) {
    generateQuadrature();
}

const std::vector<Direction>& AngularQuadrature::getDirections() const {
    return directions_;
}

void AngularQuadrature::generateQuadrature() {
    int n_theta = thetaOrder_;
    if (n_theta < 1) {
        // Logger::error("Invalid theta order: " + std::to_string(n_theta));
        throw std::invalid_argument("Theta order must be at least 1.");
    }
    int n_phi = phiOrder_;
    if (n_phi < 1) {
        // Logger::error("Invalid phi order: " + std::to_string(n_phi));
        throw std::invalid_argument("Phi order must be at least 1.");
    }

    // Generate gauss-legendre for theta
    std::vector<double> roots = Quadrature::legendreRoots(n_theta);
    std::vector<double> weights_legendre = Quadrature::legendreWeights(roots, n_theta);

    // Generate gauss-chebyshev for phi
    auto [phi_points, phi_weights] = Quadrature::gaussChebyshev(n_phi);

    // Generate directions with a tensor product of the two quadratures
    for(int i = 0; i < n_theta; ++i) {
        for(int j = 0; j < n_phi; ++j) {
            double phi = phi_points[j];
            double weight = weights_legendre[i] * phi_weights[j];
            directions_.emplace_back(Direction{roots[i], phi, weight});
        }
    }
    // Logger::info("Angular quadrature generated successfully. Number of directions: " + std::to_string(directions_.size()));

    // Vérifier la normalisation des poids
    double total_weight = 0.0;
    for(const auto& dir : directions_) {
        total_weight += dir.weight;
    }
    // Logger::info("Sum of weights: " + std::to_string(total_weight) + " (should be close to 4*pi ≈ " + std::to_string(4.0 * M_PI) + ")");
}

} // namespace SNSolver
