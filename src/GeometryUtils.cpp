// src/GeometryUtils.cpp
#include "GeometryUtils.hpp"
#include "Vector3D.hpp"
#include <cmath>
#include <random>
#include <thread>

// Modified computeFaceNormal to ensure outward pointing normals
Vector3D computeFaceNormal(const std::array<Vector3D, 3>& triangle, const Vector3D& cell_center) {
    Vector3D u = triangle[1] - triangle[0];
    Vector3D v = triangle[2] - triangle[0];
    Vector3D n = u.cross(v);

    double norm = n.norm();
    if(norm == 0) return {0.0, 0.0, 0.0}; // Degenerate triangle

    n = n.normalized();

    // Compute a vector from the cell center to a point on the face
    Vector3D vector_to_face = triangle[0] - cell_center;

    // If the dot product is positive, the normal points outward; if negative, flip it
    if(n.dot(vector_to_face) < 0) {
        n = -n;
    }

    return n;
}


Vector3D samplePointOnTriangle(const std::array<Vector3D, 3>& triangle) {
    // Thread-local random number generator for thread safety
    thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    double u = distribution(generator);
    double v = distribution(generator);

    // Ensure the point is inside the triangle
    if(u + v > 1.0) {
        u = 1.0 - u;
        v = 1.0 - v;
    }

    double x = triangle[0].x + u * (triangle[1].x - triangle[0].x) + v * (triangle[2].x - triangle[0].x);
    double y = triangle[0].y + u * (triangle[1].y - triangle[0].y) + v * (triangle[2].y - triangle[0].y);
    double z = triangle[0].z + u * (triangle[1].z - triangle[0].z) + v * (triangle[2].z - triangle[0].z);

    return Vector3D(x, y, z);
}