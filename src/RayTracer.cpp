// src/RayTracer.cpp

#include "RayTracer.hpp"
#include "Tetrahedron.hpp" // Ensure this includes Tetrahedron class with findExit method
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cassert>

// Constructor for variable direction tracing
RayTracer::RayTracer(const MeshHandler& mesh_handler, const Field& field_handler)
    : mesh_(mesh_handler),
      mode_(RayTracerMode::VARIABLE_DIRECTION),
      field_ptr_(&field_handler),
      fixed_direction_(0.0, 0.0, 0.0) // Initialize to zero vector (unused)
{
    // Initialization if needed
}

// Constructor for constant direction tracing
RayTracer::RayTracer(const MeshHandler& mesh_handler, const Vector3D& fixed_direction)
    : mesh_(mesh_handler),
      mode_(RayTracerMode::CONSTANT_DIRECTION),
      field_ptr_(nullptr),
      fixed_direction_(fixed_direction.normalized()) // Normalize the fixed direction
{
    // Ensure fixed_direction is not zero
    assert(!(fixed_direction_.x == 0.0 && fixed_direction_.y == 0.0 && fixed_direction_.z == 0.0) && "Fixed direction cannot be zero.");
}


std::vector<CellTrace> RayTracer::traceRay(int start_cell_id, const Vector3D& start_point, int max_iter) const {
    std::vector<CellTrace> pathline;
    int current_cell_id = start_cell_id;
    Vector3D current_point = start_point;
    // double total_time = 0.0;

    for(int iter = 0; iter < max_iter; ++iter) {
        // Validate current_cell_id
        if(current_cell_id < 0 || current_cell_id >= static_cast<int>(mesh_.getCells().size())) {
            // std::cerr << "Error: Invalid current cell ID: " << current_cell_id << std::endl;
            Logger::error("Invalid current cell ID: " + std::to_string(current_cell_id));
            break;
        }

        // Retrieve the current cell
        const TetraCell& cell = mesh_.getCells()[current_cell_id];

        // Determine the velocity vector based on mode
        Vector3D v;
        if(mode_ == RayTracerMode::VARIABLE_DIRECTION) {
            // Retrieve the velocity from the field's vector field
            const Vector3D& field_val = field_ptr_->getVectorFields()[current_cell_id];
            v = field_val;
        } else { // CONSTANT_DIRECTION
            v = fixed_direction_;
        }

        // Initialize the Tetrahedron with current cell data
        Tetrahedron tetra(cell, mesh_.getNodes(), v);

        double t_exit;
        Vector3D x_exit;
        int exit_face_id;

        // Find the exit point and corresponding face
        bool has_exit = tetra.findExit(current_point, v, t_exit, x_exit, exit_face_id);
        if(!has_exit) {
            // std::cerr << "Warning: No exit found for cell " << current_cell_id << " at iteration " << iter << std::endl;
            // also cout the current point and velocity
            Logger::error("No exit found for cell " + std::to_string(current_cell_id) + " at iteration " + std::to_string(iter));
            // Logger::info("Current point " + std::to_string(current_point));
            // std::cout << "Current point: " << current_point << std::endl;
            // std::cout << "Velocity: " << v << std::endl;
            break;
        }

        // Record the traversal segment
        CellTrace segment;
        segment.cell_id = current_cell_id;
        segment.time_spent = t_exit;
        segment.start_point = current_point;
        segment.end_point = x_exit;
        pathline.push_back(segment);
        // total_time += t_exit;

        // Retrieve the neighboring cell using MeshHandler's method
        int neighbor_cell_id = mesh_.getNeighborCell(current_cell_id, exit_face_id);
        if(neighbor_cell_id == -1) {
            // Ray has exited the domain
            // std::cout << "Info: Ray exited the domain at iteration " << iter << std::endl;
            break;
        }

        // Update for the next iteration
        current_cell_id = neighbor_cell_id;
        current_point = x_exit;
    }

    // std::cout << "Info: Ray tracing completed. Total time elapsed: " << total_time << std::endl;
    return pathline;
}