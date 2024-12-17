// src/RayTracerManager.cpp
#include "RayTracerManager.hpp"
#include "GeometryUtils.hpp"
#include <cmath>
#include <omp.h>
#include <iostream> // For logging purposes

RayTracerManager::RayTracerManager(const MeshHandler& mesh,
                                   const Field& field,
                                   RayTracer& tracer,
                                   AngularQuadrature& angular_quadrature)
    : mesh_(mesh),
      field_(field),
      tracer_(tracer),
      angular_quadrature_(angular_quadrature) {}

bool RayTracerManager::isValidDirection(const std::array<double, 3>& face_normal,
                                        const Direction& dir,
                                        double angle_threshold_deg) const {
    // Convert mu and phi to Cartesian coordinates
    double theta = std::acos(dir.mu);
    double x = std::sin(theta) * std::cos(dir.phi);
    double y = std::sin(theta) * std::sin(dir.phi);
    double z = dir.mu;

    // Normalize face normal and direction vector
    double norm_face = std::sqrt(face_normal[0]*face_normal[0] +
                                 face_normal[1]*face_normal[1] +
                                 face_normal[2]*face_normal[2]);
    double norm_dir = std::sqrt(x*x + y*y + z*z);
    if(norm_face == 0 || norm_dir == 0) {
        std::cerr << "Warning: Zero-length normal or direction vector." << std::endl;
        return false; // Invalid normals
    }

    double dot = (face_normal[0] * x + face_normal[1] * y + face_normal[2] * z) / (norm_face * norm_dir);

    // Clamp dot to valid range to avoid numerical issues
    dot = std::max(-1.0, std::min(1.0, dot));

    // Compute angle in radians
    double angle_rad = std::acos(dot);

    // Define a forward cone: angle less than 90 - threshold
    return (dot > std::cos((90.0 - angle_threshold_deg) * M_PI / 180.0));
}

void RayTracerManager::generateTrackingData(int rays_per_face) {
    // Retrieve boundary faces
    std::vector<Face> boundary_faces = mesh_.getBoundaryFaces();
    if(boundary_faces.empty()) {
        std::cerr << "Warning: No boundary faces found." << std::endl;
        return;
    }

    // Retrieve directions from angular quadrature
    const std::vector<Direction>& directions = angular_quadrature_.getDirections();
    if(directions.empty()) {
        std::cerr << "Error: No directions available from AngularQuadrature." << std::endl;
        return;
    }

    // Utilize only half of the directions due to symmetry
    int half_directions = directions.size() / 2;
    if(half_directions == 0) {
        std::cerr << "Error: Insufficient directions for ray tracing." << std::endl;
        return;
    }

    // Determine the number of threads
    int num_threads = omp_get_max_threads();

    // Precompute total expected rays for memory reservation
    size_t total_rays = boundary_faces.size() * rays_per_face * half_directions;
    tracking_data_.reserve(total_rays);

    // Initialize thread-local tracking data
    std::vector<std::vector<TrackingData>> thread_local_data(num_threads);
    for(auto& local_data : thread_local_data) {
        local_data.reserve(rays_per_face * half_directions * (boundary_faces.size() / num_threads + 1));
    }

    // Parallelize over boundary faces
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::vector<TrackingData>& local_tracking_data = thread_local_data[thread_id];

        #pragma omp for schedule(dynamic)
        for(int i = 0; i < static_cast<int>(boundary_faces.size()); ++i) {
            const Face& face = boundary_faces[i];

            // Validate face node indices
            if(face.n0 >= mesh_.getNodes().size() ||
               face.n1 >= mesh_.getNodes().size() ||
               face.n2 >= mesh_.getNodes().size()) {
                std::cerr << "Warning: Face node index out of bounds. Face ID: " << i << std::endl;
                continue; // Skip invalid face
            }

            std::array<Node, 3> face_nodes = {
                mesh_.getNodes()[face.n0],
                mesh_.getNodes()[face.n1],
                mesh_.getNodes()[face.n2]
            };

            // Compute face normal
            Vector3D face_normal = computeFaceNormal(face_nodes);
            double norm = face_normal.norm();
            if(norm == 0) {
                std::cerr << "Warning: Degenerate face encountered. Face ID: " << i << std::endl;
                continue; // Skip degenerate face
            }

            // Determine the starting cell ID (only one adjacent cell for boundary faces)
            int start_cell_id = mesh_.getFaceAdjacentCell(face, /*only_one=*/true);
            if(start_cell_id == -1) {
                std::cerr << "Warning: Unable to determine starting cell for face ID: " << i << std::endl;
                continue; // Invalid face
            }

            // Sample starting points on the face
            for(int ray = 0; ray < rays_per_face; ++ray) {
                Vector3D start_point = samplePointOnTriangle(face_nodes);

                // Iterate over half directions
                for(int d = 0; d < half_directions; ++d) {
                    const Direction& dir = directions[d];

                    // Check if direction is valid (not parallel and outward)
                    if(isValidDirection(face_normal, dir, /*threshold=*/1.0)) { // 1 degree threshold
                        // Convert to Cartesian direction vector
                        double theta = std::acos(dir.mu);
                        double x_dir = std::sin(theta) * std::cos(dir.phi);
                        double y_dir = std::sin(theta) * std::sin(dir.phi);
                        double z_dir = dir.mu;

                        // Normalize direction vector
                        double norm_dir = std::sqrt(x_dir * x_dir + y_dir * y_dir + z_dir * z_dir);
                        if(norm_dir == 0) {
                            std::cerr << "Warning: Zero-length direction vector. Ray ID will not be unique." << std::endl;
                            continue; // Skip invalid direction
                        }
                        std::array<double, 3> direction = {x_dir / norm_dir, y_dir / norm_dir, z_dir / norm_dir};

                        // Trace the ray
                        std::vector<CellTrace> cell_traces = tracer_.traceRay(start_cell_id, start_point, 100);
                        if(cell_traces.empty()) {
                            std::cerr << "Warning: Ray tracing returned empty cell traces. Ray ID: " << ray_counter_ << std::endl;
                            continue; // Skip rays with no traversal
                        }

                        // Populate TrackingData
                        TrackingData data;
                        // Assign a unique ray_id based on thread and local counter
                        data.ray_id = thread_id * (rays_per_face * half_directions) + ray * half_directions + d;
                        data.direction = direction; // Assuming TrackingData has direction field
                        data.cell_traces = std::move(cell_traces);

                        // Add to local tracking data
                        local_tracking_data.push_back(std::move(data));
                    }
                }
            }
        }
    }

    // Merge all thread-local data into the main tracking_data_ vector
    for(auto& local_data : thread_local_data) {
        tracking_data_.insert(tracking_data_.end(),
                              std::make_move_iterator(local_data.begin()),
                              std::make_move_iterator(local_data.end()));
    }
}