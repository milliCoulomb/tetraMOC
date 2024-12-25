// src/RayTracer.cpp

#include "RayTracer.hpp"
#include "Tetrahedron.hpp" // Ensure this includes Tetrahedron class with findExit method
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>
#include <cassert>

const Vector3D ZERO_VECTOR(0.0, 0.0, 0.0);

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
RayTracer::RayTracer(const MeshHandler& mesh_handler, const Vector3D& fixed_direction, const double direction_weight)
    : mesh_(mesh_handler),
      mode_(RayTracerMode::CONSTANT_DIRECTION),
      field_ptr_(nullptr),
      fixed_direction_(fixed_direction.normalized()), // Normalize the fixed direction
      direction_weight_(direction_weight)
{
    // Ensure fixed_direction is not zero
    if (fixed_direction_.isAlmostEqual(ZERO_VECTOR)) {
        Logger::error("Fixed direction cannot be zero.");
        throw std::invalid_argument("Fixed direction cannot be zero.");
    }
}


std::vector<CellTrace> RayTracer::traceRay(int start_cell_id, const Vector3D& start_point, int max_iter) const {
    const int MAX_ITERATIONS = 1000;
    std::vector<CellTrace> pathline;
    pathline.reserve(max_iter); // Preallocate memory

    int current_cell_id = start_cell_id;
    Vector3D current_point = start_point;

    for(int iter = 0; iter < max_iter && iter < MAX_ITERATIONS; ++iter) {
        // Validate current_cell_id
        if(current_cell_id < 0 || current_cell_id >= static_cast<int>(mesh_.getCells().size())) {
            Logger::error("Invalid current cell ID: " + std::to_string(current_cell_id));
            break;
        }

        // Retrieve the current cell
        const TetraCell& cell = mesh_.getCells()[current_cell_id];

        // Determine the velocity vector based on mode
        Vector3D v;
        if(mode_ == RayTracerMode::VARIABLE_DIRECTION) {
            assert(field_ptr_ != nullptr && "field_ptr_ is null in VARIABLE_DIRECTION mode.");
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
            break; // reached the end of the domain
        }

        // Record the traversal segment using emplace_back
        pathline.emplace_back(CellTrace{current_cell_id, t_exit, current_point, x_exit});

        // Validate exit_face_id
        if(exit_face_id < 0 || exit_face_id >= mesh_.getTotalFaces()) {
            Logger::error("Invalid exit_face_id: " + std::to_string(exit_face_id) + " at iteration " + std::to_string(iter));
            break;
        }

        // Retrieve the neighboring cell using MeshHandler's method
        int neighbor_cell_id = mesh_.getNeighborCell(current_cell_id, exit_face_id);
        if(neighbor_cell_id == -1) {
            Logger::info("Ray exited the domain at iteration " + std::to_string(iter));
            break;
        }

        // Validate the neighbor cell ID
        if(neighbor_cell_id < 0 || neighbor_cell_id >= static_cast<int>(mesh_.getCells().size())) {
            Logger::error("Invalid neighbor_cell_id: " + std::to_string(neighbor_cell_id) + " at iteration " + std::to_string(iter));
            break;
        }

        // Update for the next iteration
        current_cell_id = neighbor_cell_id;
        current_point = x_exit;
    }

    if(pathline.size() >= static_cast<size_t>(max_iter)) {
        Logger::warning("Ray tracing reached maximum iterations (" + std::to_string(max_iter) + ").");
    }

    return pathline;
}