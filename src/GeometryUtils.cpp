// src/GeometryUtils.cpp
#include "GeometryUtils.hpp"
#include <cmath>
#include <random>

namespace SNSolver {

std::array<double, 3> computeFaceNormal(const std::array<Node, 3>& triangle) {
    // Vectors along two edges
    double ux = triangle[1].x - triangle[0].x;
    double uy = triangle[1].y - triangle[0].y;
    double uz = triangle[1].z - triangle[0].z;

    double vx = triangle[2].x - triangle[0].x;
    double vy = triangle[2].y - triangle[0].y;
    double vz = triangle[2].z - triangle[0].z;

    // Cross product to get the normal
    double nx = uy * vz - uz * vy;
    double ny = uz * vx - ux * vz;
    double nz = ux * vy - uy * vx;

    // Normalize the normal vector
    double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    if (norm < 1e-8) { // Avoid division by zero
        return {0.0, 0.0, 0.0};
    }
    return {nx / norm, ny / norm, nz / norm};
}

std::array<double, 3> samplePointOnTriangle(const std::array<Node, 3>& triangle) {
    static thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double u = dis(gen);
    double v = dis(gen);

    if (u + v > 1.0) {
        u = 1.0 - u;
        v = 1.0 - v;
    }

    double x = u * triangle[0].x + v * triangle[1].x + (1.0 - u - v) * triangle[2].x;
    double y = u * triangle[0].y + v * triangle[1].y + (1.0 - u - v) * triangle[2].y;
    double z = u * triangle[0].z + v * triangle[1].z + (1.0 - u - v) * triangle[2].z;

    return {x, y, z};
}

} // namespace SNSolver