// src/RayTracer.cpp

#include "RayTracer.hpp"
#include "Tetrahedron.hpp" // Ensure this includes Tetrahedron class with findExit method
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>

// Implementation of RayTracer constructor
RayTracer::RayTracer(const MeshHandler& mesh_handler, const Field& field_handler)
    : mesh(mesh_handler), field(field_handler) {
    // Initialization if needed
}

// Implementation of traceRay
std::vector<CellTrace> RayTracer::traceRay(int start_cell_id, const Vector3D& start_point, int max_iter) const {
    std::vector<CellTrace> pathline;
    int current_cell_id = start_cell_id;
    Vector3D current_point = start_point;
    double total_time = 0.0;

    for(int iter = 0; iter < max_iter; ++iter) {
        // Validate current_cell_id
        if(current_cell_id < 0 || current_cell_id >= static_cast<int>(mesh.getCells().size())) {
            std::cerr << "Error: Invalid current cell ID: " << current_cell_id << std::endl;
            break;
        }

        // Retrieve the current cell and its associated vector field
        const TetraCell& cell = mesh.getCells()[current_cell_id];
        const CellVectorField& field_val = field.getVectorFields()[current_cell_id];

        // Initialize the Tetrahedron with current cell data
        Tetrahedron tetra(cell, mesh.getNodes(), field_val);

        // Define the velocity vector
        Vector3D v(field_val.vx, field_val.vy, field_val.vz);

        double t_exit;
        Vector3D x_exit;
        int exit_face_id;

        // Find the exit point and corresponding face
        bool has_exit = tetra.findExit(current_point, v, t_exit, x_exit, exit_face_id);
        if(!has_exit) {
            std::cerr << "Warning: No exit found for cell " << current_cell_id << " at iteration " << iter << std::endl;
            break;
        }

        // Record the traversal segment
        CellTrace segment;
        segment.cell_id = current_cell_id;
        segment.time_spent = t_exit;
        segment.start_point = current_point;
        segment.end_point = x_exit;
        pathline.push_back(segment);
        total_time += t_exit;

        // Retrieve the neighboring cell using MeshHandler's method
        int neighbor_cell_id = mesh.getNeighborCell(current_cell_id, exit_face_id);
        if(neighbor_cell_id == -1) {
            // Ray has exited the domain
            std::cout << "Info: Ray exited the domain at iteration " << iter << std::endl;
            break;
        }

        // Update for the next iteration
        current_cell_id = neighbor_cell_id;
        current_point = x_exit;
    }

    std::cout << "Info: Ray tracing completed. Total time elapsed: " << total_time << std::endl;
    return pathline;
}
