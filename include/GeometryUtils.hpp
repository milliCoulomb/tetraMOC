// include/GeometryUtils.hpp
#ifndef GEOMETRY_UTILS_HPP
#define GEOMETRY_UTILS_HPP

#include <array>
#include "MeshHandler.hpp" // Include MeshHandler to access the Node struct
#include "Vector3D.hpp"    // Include Vector3D to access the Vector3D class

// Function to compute the normal of a triangle face
Vector3D computeFaceNormal(const std::array<Vector3D, 3>& triangle);

// Function to sample a point uniformly within a triangle
Vector3D samplePointOnTriangle(const std::array<Vector3D, 3>& triangle);

#endif // GEOMETRY_UTILS_HPP 