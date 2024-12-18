// src/Tetrahedron.cpp

#include "Tetrahedron.hpp"
#include <cmath>
#include <limits>
#include <iostream>

// Removed: using namespace SNSolver;

Tetrahedron::Tetrahedron(const TetraCell& cell, const std::vector<Vector3D>& nodes, const CellVectorField& field) {
    for(int i = 0; i < 4; ++i) {
        vertices[i] = Vector3D(nodes[cell.node_ids[i]].x, nodes[cell.node_ids[i]].y, nodes[cell.node_ids[i]].z);
    }
    // Initialize the velocity vector
    velocity = Vector3D(field.x, field.y, field.z);
}

bool Tetrahedron::findExit(const Vector3D& x0, const Vector3D& v, double& t_exit, Vector3D& x_exit, int& exit_face_id) const {
    t_exit = std::numeric_limits<double>::max();
    exit_face_id = -1;

    for(int i = 0; i < 4; ++i) {
        // Define the three vertices of the face
        std::array<Vector3D, 3> face_vertices;
        switch(i) {
            case 0:
                face_vertices = { vertices[0], vertices[1], vertices[2] };
                break;
            case 1:
                face_vertices = { vertices[0], vertices[1], vertices[3] };
                break;
            case 2:
                face_vertices = { vertices[0], vertices[2], vertices[3] };
                break;
            case 3:
                face_vertices = { vertices[1], vertices[2], vertices[3] };
                break;
            default:
                continue;
        }

        // Compute the plane of the face
        Vector3D p1 = face_vertices[0];
        Vector3D p2 = face_vertices[1];
        Vector3D p3 = face_vertices[2];

        Vector3D edge1 = p2 - p1;
        Vector3D edge2 = p3 - p1;
        Vector3D normal = edge1.cross(edge2).normalized();

        double D = normal.dot(p1);

        double denom = normal.dot(v);
        if(std::abs(denom) < 1e-12) {
            continue; // Ray parallel to face
        }

        double t = (D - normal.dot(x0)) / denom;
        if(t <= 1e-10) {
            continue; // Intersection behind the start point or too close
        }

        // Compute the intersection point
        Vector3D intersect = x0 + v * t;

        // Check if the point is inside the triangle using barycentric coordinates
        Vector3D v0 = p3 - p1;
        Vector3D v1 = p2 - p1;
        Vector3D v2 = intersect - p1;

        double dot00 = v0.dot(v0);
        double dot01 = v0.dot(v1);
        double dot02 = v0.dot(v2);
        double dot11 = v1.dot(v1);
        double dot12 = v1.dot(v2);

        double denom_bary = dot00 * dot11 - dot01 * dot01;
        if(denom_bary == 0.0) {
            continue; // Degenerate triangle
        }

        double u = (dot11 * dot02 - dot01 * dot12) / denom_bary;
        double v_bary = (dot00 * dot12 - dot01 * dot02) / denom_bary;

        if(u >= 0.0 && v_bary >= 0.0 && (u + v_bary) <= 1.0) {
            if(t < t_exit) {
                t_exit = t;
                x_exit = { intersect.x, intersect.y, intersect.z };
                exit_face_id = i;
            }
        }
    }

    if(exit_face_id != -1) {
        return true;
    }
    return false;
}