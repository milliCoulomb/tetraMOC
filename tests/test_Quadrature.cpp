// tests/test_Quadrature.cpp
#include "Quadrature.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace Quadrature {

// roots for n = 2
TEST(QuadratureTest, GaussLegendreRootsN2) {
    int n = 2;
    std::vector<double> expected_roots = { -std::sqrt(1.0/3.0), std::sqrt(1.0/3.0) };
    std::vector<double> roots = legendreRoots(n);
    ASSERT_EQ(roots.size(), expected_roots.size());
    for(size_t i = 0; i < roots.size(); ++i) {
        EXPECT_NEAR(roots[i], expected_roots[i], 1e-10);
    }
}

// weights for n = 2
TEST(QuadratureTest, GaussLegendreWeightsN2) {
    int n = 2;
    std::vector<double> roots = legendreRoots(n);
    std::vector<double> weights = legendreWeights(roots, n);
    std::vector<double> expected_weights = { 1.0, 1.0 };
    ASSERT_EQ(weights.size(), expected_weights.size());
    for(size_t i = 0; i < weights.size(); ++i) {
        EXPECT_NEAR(weights[i], expected_weights[i], 1e-10);
    }
}

// points for n = 3
TEST(QuadratureTest, GaussChebyshevPointsN3) {
    int n = 3;
    std::vector<double> expected_points = { M_PI/6, M_PI/2, 5*M_PI/6 };
    auto [points, weights] = gaussChebyshev(n);
    ASSERT_EQ(points.size(), expected_points.size());
    for(size_t i = 0; i < points.size(); ++i) {
        EXPECT_NEAR(points[i], expected_points[i], 1e-10);
    }
}

// weights for n = 3
TEST(QuadratureTest, GaussChebyshevWeightsN3) {
    int n = 3;
    std::vector<double> expected_weights = { M_PI/3, M_PI/3, M_PI/3 };
    auto [points, weights] = gaussChebyshev(n);
    ASSERT_EQ(weights.size(), expected_weights.size());
    for(size_t i = 0; i < weights.size(); ++i) {
        EXPECT_NEAR(weights[i], expected_weights[i], 1e-10);
    }
}

// assess roots for n = 3
TEST(QuadratureTest, GaussLegendreRootsN3) {
    int n = 3;
    std::vector<double> expected_roots = { -std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0) };
    std::vector<double> roots = legendreRoots(n);
    ASSERT_EQ(roots.size(), expected_roots.size());
    for(size_t i = 0; i < roots.size(); ++i) {
        EXPECT_NEAR(roots[i], expected_roots[i], 1e-10);
    }
}

// assess weights for n = 3
TEST(QuadratureTest, GaussLegendreWeightsSumN3) {
    int n = 3;
    std::vector<double> roots = legendreRoots(n);
    std::vector<double> weights = legendreWeights(roots, n);
    double sum_weights = 0.0;
    for(auto w : weights) sum_weights += w;
    EXPECT_NEAR(sum_weights, 2.0, 1e-10); // Intégrale de P0(x) = 2
}

} // namespace Quadrature

// // main.cpp pour les tests (si nécessaire)
// int main(int argc, char **argv) {
//     ::testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }
