// tests/test_AngularQuadrature.cpp
#include "AngularQuadrature.hpp"
#include "Quadrature.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace SNSolver {

// angular quadrature (theta_order = 2, phi_order = 4)
TEST(AngularQuadratureTest, GenerateQuadratureSnOrder2) {
    int theta_order = 2;
    int phi_order = 4;
    AngularQuadrature aq(theta_order, phi_order);
    const std::vector<Direction>& directions = aq.getDirections();
    
    // Calcul attendu : theta_order = 2, phi_order = 4 => 2 * 4 = 8 directions
    EXPECT_EQ(directions.size(), 8);

    // Vérifier les mu et phi des directions
    // Pour snOrder = 2, les racines Gauss-Legendre sont ±sqrt(1/3)
    std::vector<double> expected_mu = { -std::sqrt(1.0/3.0), std::sqrt(1.0/3.0) };
    std::vector<double> expected_phi = { M_PI/4, M_PI/2, 3*M_PI/4, M_PI }; // Si n_phi = 4

    // Puisque phi_order = 4, les points de Gauss-Chebyshev pour n_phi = 4 sont pi*(2k +1)/(2n)
    std::vector<double> expected_phi_calculated = { M_PI * 1.0 / 8.0, M_PI * 3.0 / 8.0, M_PI * 5.0 / 8.0, M_PI * 7.0 / 8.0 };

    // Vérifier que chaque combinaison mu x phi est présente
    for(int i = 0; i < snOrder_; ++i) { // snOrder_ = 2
        for(int j = 0; j < snOrder_ * 2; ++j) { // phi_order = 4
            int index = i * (snOrder_ * 2) + j;
            EXPECT_NEAR(directions[index].mu, expected_mu[i], 1e-10);
            EXPECT_NEAR(directions[index].phi, expected_phi_calculated[j], 1e-10);
        }
    }

    // Vérifier la somme des poids
    double total_weight = 0.0;
    for(const auto& dir : directions) {
        total_weight += dir.weight;
    }
    EXPECT_NEAR(total_weight, 4.0 * M_PI, 1e-6);
}

//n_theta = 1, n_phi = 2)
TEST(AngularQuadratureTest, GenerateQuadratureSnOrder1) {
    int n_theta = 1;
    int n_phi = 2;
    AngularQuadrature aq(n_theta, n_phi);
    const std::vector<Direction>& directions = aq.getDirections();
    
    // Calcul attendu : theta_order = 1, phi_order = 2 => 1 * 2 = 2 directions
    EXPECT_EQ(directions.size(), 2);

    // Vérifier les mu et phi des directions
    // Pour snOrder = 1, les racines Gauss-Legendre sont cos(0) = 1
    std::vector<double> expected_mu = { std::cos(0.0) };
    std::vector<double> expected_phi_calculated = { M_PI/2, 3*M_PI/2 }; // Gauss-Chebyshev n_phi = 2

    // Vérifier chaque direction
    EXPECT_NEAR(directions[0].mu, expected_mu[0], 1e-10);
    EXPECT_NEAR(directions[0].phi, expected_phi_calculated[0], 1e-10);
    
    EXPECT_NEAR(directions[1].mu, expected_mu[0], 1e-10);
    EXPECT_NEAR(directions[1].phi, expected_phi_calculated[1], 1e-10);

    // Vérifier la somme des poids
    double total_weight = 0.0;
    for(const auto& dir : directions) {
        total_weight += dir.weight;
    }
    EXPECT_NEAR(total_weight, 4.0 * M_PI, 1e-6);
}

// testing normalization of weights for n_theta = 2, n_phi = 4
TEST(AngularQuadratureTest, WeightNormalizationSnOrder3) {
    int n_theta = 2;
    int n_phi = 4;
    AngularQuadrature aq(n_theta, n_phi);
    const std::vector<Direction>& directions = aq.getDirections();
    
    // Vérifier la somme des poids
    double total_weight = 0.0;
    for(const auto& dir : directions) {
        total_weight += dir.weight;
    }
    EXPECT_NEAR(total_weight, 4.0 * M_PI, 1e-6);
}

// Test de consistance entre Quadrature et AngularQuadrature pour snOrder = 4
TEST(AngularQuadratureTest, ConsistencyWithQuadratureSnOrder4) {
    int snOrder = 4;
    AngularQuadrature aq(snOrder);
    const std::vector<Direction>& directions = aq.getDirections();

    // Vérifier que la somme des poids est proche de 4*pi
    double total_weight = 0.0;
    for(const auto& dir : directions) {
        total_weight += dir.weight;
    }
    EXPECT_NEAR(total_weight, 4.0 * M_PI, 1e-6);
}

// Test de gestion des ordres invalides (snOrder < 1)
TEST(AngularQuadratureTest, InvalidSnOrder) {
    EXPECT_THROW(SNSolver::AngularQuadrature aq(0), std::invalid_argument);
    EXPECT_THROW(SNSolver::AngularQuadrature aq(-2), std::invalid_argument);
}

} // namespace SNSolver

// main.cpp pour les tests (si nécessaire)
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
