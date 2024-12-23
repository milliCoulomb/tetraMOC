// tests/test_AngularQuadrature.cpp
#include "AngularQuadrature.hpp"
#include "Quadrature.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "Vector3D.hpp"

// angular quadrature (theta_order = 2, phi_order = 4)
TEST(AngularQuadratureTest, GenerateQuadratureSnOrder2) {
    int theta_order = 2;
    int phi_order = 4;
    AngularQuadrature aq(theta_order, phi_order);
    const std::vector<Direction>& directions = aq.getDirections();
    
    EXPECT_EQ(directions.size(), 8);

    std::vector<double> expected_mu = {-std::sqrt(1.0/3.0), std::sqrt(1.0/3.0)};
    std::vector<double> expected_phi = { M_PI/4.0, 3.0 * M_PI/4.0, 5.0*M_PI/4.0, 7.0 / 4.0 * M_PI }; // Si n_phi = 4

    for(int i = 0; i < theta_order; ++i) { // snOrder_ = 2
    
        for(int j = 0; j < phi_order; ++j) { // phi_order = 4
            int index = i * phi_order + j;
            EXPECT_NEAR(directions[index].mu, expected_mu[i], 1e-10);
            EXPECT_NEAR(directions[index].phi, expected_phi[j], 1e-10);
        }
    }

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
    std::vector<double> expected_mu = { 0.0 };
    std::vector<double> expected_phi_calculated = { M_PI/2.0, 3.0*M_PI/2.0 }; // Gauss-Chebyshev n_phi = 2

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
    int n_theta = 4;
    int n_phi = 4;
    AngularQuadrature aq(n_theta, n_phi);
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
    EXPECT_THROW(AngularQuadrature aq(0, 4), std::invalid_argument);
    EXPECT_THROW(AngularQuadrature aq(-2, 4), std::invalid_argument);
}

// we need to check if the angular quadrature is symmetric
TEST(AngularQuadratureTest, Symmetry) {
    int n_theta = 2;
    int n_phi = 2;
    AngularQuadrature aq(n_theta, n_phi);
    const std::vector<Direction>& directions = aq.getDirections();
    std::vector<Vector3D> directions_3D;
    // Vérifier la symétrie par rapport à l'axe z
    for(int i = 0; i < n_theta; ++i) {
        for(int j = 0; j < n_phi; ++j) {
            int index = i * n_phi + j;
            // build the 3D direction vector
            Vector3D direction(std::sqrt(1 - directions[index].mu * directions[index].mu) * std::cos(directions[index].phi),
                                 std::sqrt(1 - directions[index].mu * directions[index].mu) * std::sin(directions[index].phi),
                                 directions[index].mu);
            directions_3D.push_back(direction);
        }
    }
    // then we loop on half of the directions
    for(int i = 0; i < directions_3D.size() / 2; ++i) {
        // we check if the direction is the opposite of the other
        EXPECT_NEAR(directions_3D[i].x, -directions_3D[directions_3D.size() - i - 1].x, 1e-10);
        EXPECT_NEAR(directions_3D[i].y, -directions_3D[directions_3D.size() - i - 1].y, 1e-10);
        EXPECT_NEAR(directions_3D[i].z, -directions_3D[directions_3D.size() - i - 1].z, 1e-10);
    }
}

// add a unit test for the addDirection method
TEST(AngularQuadratureTest, AddDirection) {
    AngularQuadrature aq(2, 4);
    Direction d = {0.5, 0.5, 0.5};
    aq.addDirection(d);
    const std::vector<Direction>& directions = aq.getDirections();
    EXPECT_EQ(directions.size(), 9);
    EXPECT_EQ(directions[8].mu, 0.5);
    EXPECT_EQ(directions[8].phi, 0.5);
    EXPECT_EQ(directions[8].weight, 0.5);
}

// test the overloaded constructor with a vector of directions
TEST(AngularQuadratureTest, OverloadedConstructor) {
    std::vector<Direction> directions = {{0.5, 0.5, 0.5}, {0.3, 0.3, 0.3}};
    AngularQuadrature aq(directions);
    const std::vector<Direction>& directions_ = aq.getDirections();
    EXPECT_EQ(directions_.size(), 2);
    EXPECT_EQ(directions_[0].mu, 0.5);
    EXPECT_EQ(directions_[0].phi, 0.5);
    EXPECT_EQ(directions_[0].weight, 0.5);
    EXPECT_EQ(directions_[1].mu, 0.3);
    EXPECT_EQ(directions_[1].phi, 0.3);
    EXPECT_EQ(directions_[1].weight, 0.3);
}

// test the getTotalWeight method
TEST(AngularQuadratureTest, GetTotalWeight) {
    AngularQuadrature aq(4, 8);
    double total_weight = aq.getTotalWeight();
    EXPECT_NEAR(total_weight, 4.0 * M_PI, 1e-6);
}

// // main.cpp pour les tests (si nécessaire)
// int main(int argc, char **argv) {
//     ::testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }
