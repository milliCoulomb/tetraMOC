// include/GeometryUtils.hpp
#ifndef GEOMETRY_UTILS_HPP
#define GEOMETRY_UTILS_HPP

#include <array>
#include "MeshHandler.hpp" // Include MeshHandler to access the Node struct

// Function to compute the normal of a triangle face
std::array<double, 3> computeFaceNormal(const std::array<Node, 3>& triangle);

// Function to sample a point uniformly within a triangle
std::array<double, 3> samplePointOnTriangle(const std::array<Node, 3>& triangle);

#endif // GEOMETRY_UTILS_HPP 