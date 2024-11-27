// AngularQuadrature.cpp
#include "AngularQuadrature.hpp"
#include "Logger.hpp" // Assure-toi que Logger est implémenté
#include <cmath>

namespace SNSolver {

AngularQuadrature::AngularQuadrature(int snOrder) : snOrder_(snOrder) {
    generateQuadrature();
}

const std::vector<Direction>& AngularQuadrature::getDirections() const {
    return directions_;
}

void AngularQuadrature::generateQuadrature() {
    int n_theta = snOrder_;
    int n_phi = snOrder_;

    // Generate gauss-legendre for theta
    std::vector<double> roots = Quadrature::legendreRoots(n_theta);
    std::vector<double> weights_legendre = Quadrature::legendreWeights(roots, n_theta);

    // Generate gauss-chebyshev for phi
    auto [phi_points, phi_weights] = Quadrature::gaussChebyshev(n_phi);

    // Generate directions with a tensor product of the two quadratures
    for(int i = 0; i < n_theta; ++i) {
        double mu = std::cos(roots[i]); // mu = cos(theta)
        for(int j = 0; j < n_phi; ++j) {
            double phi = phi_points[j];
            double weight = weights_legendre[i] * phi_weights[j];
            directions_.emplace_back(Direction{mu, phi, weight});
        }
    }
    Logger::info("Angular quadrature generated successfully. Number of directions: " + std::to_string(directions_.size()));

    // Vérifier la normalisation des poids
    double total_weight = 0.0;
    for(const auto& dir : directions_) {
        total_weight += dir.weight;
    }
    Logger::info("Sum of weights: " + std::to_string(total_weight) + " (should be close to 4*pi ≈ " + std::to_string(4.0 * M_PI) + ")");
}

} // namespace SNSolver
