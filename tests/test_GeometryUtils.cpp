// tests/test_GeometryUtils.cpp
#include <gtest/gtest.h>
#include "GeometryUtils.hpp"
#include "Vector3D.hpp"

namespace SNSolver {

// Helper function to compare two vectors with a tolerance
bool vectorsAlmostEqual(const Vector3D& v1,
                        const Vector3D& v2,
                        double tol = 1e-6) {
    return (std::abs(v1.x - v2.x) < tol) &&
           (std::abs(v1.y - v2.y) < tol) &&
           (std::abs(v1.z - v2.z) < tol);
}

TEST(GeometryUtilsTest, ComputeFaceNormal_StandardTriangle) {
    // Define a right-angled triangle in the XY plane
    Vector3D n0(0.0, 0.0, 0.0);
    Vector3D n1(1.0, 0.0, 0.0);
    Vector3D n2(0.0, 1.0, 0.0);
    
    std::array<Vector3D, 3> triangle = {n0, n1, n2};
    
    Vector3D normal = computeFaceNormal(triangle, Vector3D(0.0, 0.0, 0.0));
    
    // Expected normal is along +Z axis
    Vector3D expected_normal(0.0, 0.0, 1.0);
    
    EXPECT_TRUE(vectorsAlmostEqual(normal, expected_normal)) << "Normal of standard triangle should be along +Z axis.";
}

TEST(GeometryUtilsTest, ComputeFaceNormal_DegenerateTriangle) {
    // Define a degenerate triangle where all points are colinear
    Vector3D n0(0.0, 0.0, 0.0);
    Vector3D n1(1.0, 1.0, 1.0);
    Vector3D n2(2.0, 2.0, 2.0);
    
    std::array<Vector3D, 3> triangle = {n0, n1, n2};
    
    Vector3D normal = computeFaceNormal(triangle, Vector3D(0.0, 0.0, 0.0));
    
    // Expected normal is zero vector
    Vector3D expected_normal(0.0, 0.0, 0.0);
    
    EXPECT_TRUE(vectorsAlmostEqual(normal, expected_normal)) << "Normal of degenerate triangle should be zero vector.";
}

TEST(GeometryUtilsTest, SamplePointOnTriangle_PointWithinTriangle) {
    // Define a triangle
    Vector3D n0(0.0, 0.0, 0.0);
    Vector3D n1(1.0, 0.0, 0.0);
    Vector3D n2(0.0, 1.0, 0.0);
    
    std::array<Vector3D, 3> triangle = {n0, n1, n2};
    
    // Sample multiple points and verify they lie within the triangle
    for(int i = 0; i < 1000; ++i) {
        Vector3D point = samplePointOnTriangle(triangle);
        
        double u = point.x;
        double v = point.y;
        double w = 1.0 - u - v;
        
        // Check barycentric coordinates
        EXPECT_GE(u, 0.0) << "Sampled u should be >= 0";
        EXPECT_GE(v, 0.0) << "Sampled v should be >= 0";
        EXPECT_GE(w, 0.0) << "Sampled w should be >= 0";
        EXPECT_LE(u + v, 1.0) << "Sum of u and v should be <= 1";
    }
}

TEST(GeometryUtilsTest, SamplePointOnTriangle_DegenerateTriangle) {
    // Define a degenerate triangle where all points are colinear
    Vector3D n0(0.0, 0.0, 0.0);
    Vector3D n1(1.0, 1.0, 1.0);
    Vector3D n2(2.0, 2.0, 2.0);
    
    std::array<Vector3D, 3> triangle = {n0, n1, n2};
    
    // Sample a point
    Vector3D point = samplePointOnTriangle(triangle);
    
    // Since the triangle is degenerate, the point should lie on the line
    // Check if (x, y, z) satisfies y = x and z = x
    EXPECT_NEAR(point.y, point.x, 1e-6) << "For degenerate triangle, y should equal x.";
    EXPECT_NEAR(point.z, point.x, 1e-6) << "For degenerate triangle, z should equal x.";
}

} // namespace SNSolver
