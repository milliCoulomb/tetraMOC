// tests/test_Vector3D.cpp
#include "Vector3D.hpp"
#include <gtest/gtest.h>

// testing initialization of a vector
TEST(Vector3DTest, DefaultConstructor) {
    Vector3D v;
    EXPECT_DOUBLE_EQ(v.x, 0.0);
    EXPECT_DOUBLE_EQ(v.y, 0.0);
    EXPECT_DOUBLE_EQ(v.z, 0.0);
}

// testing initialization of a vector with parameters
TEST(Vector3DTest, ParameterizedConstructor) {
    Vector3D v(1.0, 2.0, 3.0);
    EXPECT_DOUBLE_EQ(v.x, 1.0);
    EXPECT_DOUBLE_EQ(v.y, 2.0);
    EXPECT_DOUBLE_EQ(v.z, 3.0);
}

// testing the addition of vectors
TEST(Vector3DTest, Addition) {
    Vector3D v1(1.0, 2.0, 3.0);
    Vector3D v2(4.0, 5.0, 6.0);
    Vector3D v3 = v1 + v2;
    EXPECT_DOUBLE_EQ(v3.x, 5.0);
    EXPECT_DOUBLE_EQ(v3.y, 7.0);
    EXPECT_DOUBLE_EQ(v3.z, 9.0);
}

// test of the subtraction of vectors
TEST(Vector3DTest, Subtraction) {
    Vector3D v1(5.0, 7.0, 9.0);
    Vector3D v2(1.0, 2.0, 3.0);
    Vector3D v3 = v1 - v2;
    EXPECT_DOUBLE_EQ(v3.x, 4.0);
    EXPECT_DOUBLE_EQ(v3.y, 5.0);
    EXPECT_DOUBLE_EQ(v3.z, 6.0);
}

// test of the multiplication of a vector by a scalar
TEST(Vector3DTest, ScalarMultiplication) {
    Vector3D v(1.0, -2.0, 3.0);
    Vector3D v2 = v * 2.5;
    EXPECT_DOUBLE_EQ(v2.x, 2.5);
    EXPECT_DOUBLE_EQ(v2.y, -5.0);
    EXPECT_DOUBLE_EQ(v2.z, 7.5);
}

// test of the scalar product
TEST(Vector3DTest, DotProduct) {
    Vector3D v1(1.0, 3.0, -5.0);
    Vector3D v2(4.0, -2.0, -1.0);
    double dot = v1.dot(v2);
    EXPECT_DOUBLE_EQ(dot, 1.0 * 4.0 + 3.0 * (-2.0) + (-5.0) * (-1.0)); // 4 - 6 + 5 = 3
    EXPECT_DOUBLE_EQ(dot, 3.0);
}

// test of the cross product
TEST(Vector3DTest, CrossProduct) {
    Vector3D v1(2.0, 3.0, 4.0);
    Vector3D v2(5.0, 6.0, 7.0);
    Vector3D cross = v1.cross(v2);
    EXPECT_DOUBLE_EQ(cross.x, (3.0 * 7.0 - 4.0 * 6.0)); // 21 - 24 = -3
    EXPECT_DOUBLE_EQ(cross.y, (4.0 * 5.0 - 2.0 * 7.0)); // 20 - 14 = 6
    EXPECT_DOUBLE_EQ(cross.z, (2.0 * 6.0 - 3.0 * 5.0)); // 12 - 15 = -3
    EXPECT_DOUBLE_EQ(cross.x, -3.0);
    EXPECT_DOUBLE_EQ(cross.y, 6.0);
    EXPECT_DOUBLE_EQ(cross.z, -3.0);
}

// test of normalization of a vector
TEST(Vector3DTest, Normalization) {
    Vector3D v(0.0, 3.0, 4.0);
    Vector3D normalized = v.normalized();
    double expectedNorm = 1.0;
    double actualNorm = normalized.norm();
    EXPECT_NEAR(actualNorm, expectedNorm, 1e-9);
}

// test of normalization of a zero vector
TEST(Vector3DTest, NormalizationOfZeroVector) {
    Vector3D v(0.0, 0.0, 0.0);
    Vector3D normalized = v.normalized();
    EXPECT_DOUBLE_EQ(normalized.x, 0.0);
    EXPECT_DOUBLE_EQ(normalized.y, 0.0);
    EXPECT_DOUBLE_EQ(normalized.z, 0.0);
}

// test the equality of vectors
TEST(Vector3DTest, Equality) {
    Vector3D v1(1.0, 2.0, 3.0);
    Vector3D v2(1.0, 2.0, 3.0);
    EXPECT_TRUE(v1 == v2);
}

// test the inequality of vectors
TEST(Vector3DTest, Inequality) {
    Vector3D v1(1.0, 2.0, 3.0);
    Vector3D v2(1.0, 2.1, 3.0);
    EXPECT_TRUE(v1 != v2);
}

// // main.cpp for tests (optional if using gtest_main)
// int main(int argc, char **argv) {
//     ::testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }
