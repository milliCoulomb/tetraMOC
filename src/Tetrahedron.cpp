// src/Tetrahedron.cpp

#include "Tetrahedron.hpp"
#include "GeometryUtils.hpp"
#include <cmath>
#include <limits>
#include <iostream>

// Define constants for better readability and maintainability
const int NUM_FACES = 4;
const std::array<std::array<int, 3>, NUM_FACES> FACE_VERTEX_INDICES = {{
    {0, 1, 2},
    {0, 1, 3},
    {0, 2, 3},
    {1, 2, 3}
}};

// Constructor for Tetrahedron
Tetrahedron::Tetrahedron(const TetraCell& cell, const std::vector<Vector3D>& nodes, const CellVectorField& field) {
    // Validate that the cell has exactly four node IDs
    if(cell.node_ids.size() != static_cast<size_t>(NUM_FACES)) {
        Logger::error("TetraCell must have exactly four node IDs.");
        throw std::invalid_argument("Invalid number of node IDs in TetraCell.");
    }

    // Initialize the vertices and center of mass
    CenterOfMass = Vector3D(0.0, 0.0, 0.0);
    for(int i = 0; i < NUM_FACES; ++i) {
        int node_id = cell.node_ids[i];
        // Validate node_id
        if(node_id < 0 || node_id >= static_cast<int>(nodes.size())) {
            Logger::error("Node ID " + std::to_string(node_id) + " is out of bounds.");
            throw std::out_of_range("Node ID out of bounds.");
        }
        vertices[i] = nodes[node_id];
        CenterOfMass += vertices[i];
    }
    CenterOfMass = CenterOfMass / 4.0; // Calculate the average

    // Initialize the velocity vector
    velocity = field;
}

bool Tetrahedron::findExit(const Vector3D& x0, const Vector3D& v, double& t_exit, Vector3D& x_exit, int& exit_face_id) const {
    t_exit = std::numeric_limits<double>::max();
    exit_face_id = -1;

    for(int i = 0; i < NUM_FACES; ++i) {
        // Retrieve the three vertices of the current face
        std::array<Vector3D, 3> face_vertices = {
            vertices[FACE_VERTEX_INDICES[i][0]],
            vertices[FACE_VERTEX_INDICES[i][1]],
            vertices[FACE_VERTEX_INDICES[i][2]]
        };
        Vector3D p1 = face_vertices[0];
        Vector3D p2 = face_vertices[1];
        Vector3D p3 = face_vertices[2];

        try {
            // Compute the normal of the face
            Vector3D normal = computeFaceNormal(face_vertices, CenterOfMass);

            double D = normal.dot(p1);
            double denom = normal.dot(v);

            if(std::abs(denom) < EPSILON) {
                continue; // Ray parallel to face
            }

            double t = (D - normal.dot(x0)) / denom;
            if(t <= EPSILON) {
                continue; // Intersection behind the start point or too close
            }

            // Compute the intersection point
            Vector3D intersect = x0 + v * t;

            // Check if the point is inside the triangle using barycentric coordinates
            Vector3D edge1 = p2 - p1;
            Vector3D edge2 = p3 - p1;
            Vector3D v2 = intersect - p1;

            double dot00 = edge1.dot(edge1);
            double dot01 = edge1.dot(edge2);
            double dot02 = edge1.dot(v2);
            double dot11 = edge2.dot(edge2);
            double dot12 = edge2.dot(v2);

            double denom_bary = dot00 * dot11 - dot01 * dot01;
            if(denom_bary < EPSILON) {
                Logger::warning("Degenerate triangle detected on face " + std::to_string(i));
                continue; // Degenerate triangle
            }

            double u = (dot11 * dot02 - dot01 * dot12) / denom_bary;
            double v_bary = (dot00 * dot12 - dot01 * dot02) / denom_bary;

            // Introduce tolerance to account for floating-point precision
            const double BARY_TOLERANCE = 1e-8;
            if(u >= -BARY_TOLERANCE && v_bary >= -BARY_TOLERANCE && (u + v_bary) <= 1.0 + BARY_TOLERANCE) {
                if(t < t_exit) {
                    t_exit = t;
                    x_exit = intersect;
                    exit_face_id = i;
                }
            }
        } catch (const std::exception& e) {
            Logger::error("Exception in findExit on face " + std::to_string(i) + ": " + std::string(e.what()));
            continue; // Skip this face and continue with others
        }
    }

    if(exit_face_id != -1) {
        return true;
    }
    return false;
}