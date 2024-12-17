// src/GeometryUtils.cpp
#include "GeometryUtils.hpp"
#include <cmath>
#include <random>
#include <thread>

std::array<double, 3> computeFaceNormal(const std::array<Node, 3>& triangle) {
    // Vectors from the first node
    double ux = triangle[1].x - triangle[0].x;
    double uy = triangle[1].y - triangle[0].y;
    double uz = triangle[1].z - triangle[0].z;

    double vx = triangle[2].x - triangle[0].x;
    double vy = triangle[2].y - triangle[0].y;
    double vz = triangle[2].z - triangle[0].z;

    // Cross product to get normal
    double nx = uy * vz - uz * vy;
    double ny = uz * vx - ux * vz;
    double nz = ux * vy - uy * vx;

    // Normalize the normal vector
    double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
    if(norm == 0) return {0.0, 0.0, 0.0}; // Degenerate triangle

    return {nx / norm, ny / norm, nz / norm};
}

std::array<double, 3> samplePointOnTriangle(const std::array<Node, 3>& triangle) {
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

    return {x, y, z};
}