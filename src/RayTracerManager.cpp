// src/RayTracerManager.cpp

#include "RayTracerManager.hpp"
#include "GeometryUtils.hpp"
#include <cmath>
#include <omp.h>
#include <iostream> // For logging purposes

RayTracerManager::RayTracerManager(const MeshHandler& mesh,
                                   const Field& base_field,
                                   AngularQuadrature& angular_quadrature)
    : mesh_(mesh),
      base_field_(base_field),
      angular_quadrature_(angular_quadrature)
{
    initializeRayTracers();
}

RayTracerManager::RayTracerManager(const MeshHandler& mesh,
                                   const Field& base_field,
                                   AngularQuadrature& angular_quadrature,
                                   const std::vector<Vector3D>& constant_directions)
    : mesh_(mesh),
      base_field_(base_field),
      angular_quadrature_(angular_quadrature)
{
    initializeRayTracers();
    initializeConstantDirectionRayTracers(constant_directions);
}

void RayTracerManager::initializeRayTracers()
{
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

    // Initialize RayTracer instances for variable direction mode
    for(int d = 0; d < half_directions; ++d) {
        const Direction& dir = directions[d];
        
        // Create a Vector3D direction vector
        Vector3D direction_vector(
            std::sqrt(1 - dir.mu * dir.mu) * std::cos(dir.phi),
            std::sqrt(1 - dir.mu * dir.mu) * std::sin(dir.phi),
            dir.mu
        );
        
        // Create a new Field instance sharing the same CellVectorField but with unique direction
        // Note: Since RayTracer in VARIABLE_DIRECTION mode uses Field's vector field, direction_vector can be arbitrary or set as needed
        // For this example, we assume direction_vector is used elsewhere or integrated into Field as per previous refactoring
        // Therefore, we initialize RayTracer in VARIABLE_DIRECTION mode without modifying the Field's direction
        
        // Create a RayTracer in VARIABLE_DIRECTION mode
        ray_tracers_.emplace_back(std::make_unique<RayTracer>(mesh_, base_field_));
    }
}

void RayTracerManager::initializeConstantDirectionRayTracers(const std::vector<Vector3D>& constant_directions)
{
    for(const auto& dir : constant_directions) {
        // Normalize the direction vector to ensure consistency
        Vector3D normalized_dir = dir.normalized();

        // Create a RayTracer in CONSTANT_DIRECTION mode
        ray_tracers_.emplace_back(std::make_unique<RayTracer>(mesh_, normalized_dir));
    }
}

bool RayTracerManager::isValidDirection(const Vector3D& face_normal,
                                        const Direction& dir,
                                        double angle_threshold_deg) const
{
    Vector3D DirVector(std::sqrt(1 - dir.mu * dir.mu) * std::cos(dir.phi),
                       std::sqrt(1 - dir.mu * dir.mu) * std::sin(dir.phi),
                       dir.mu);

    double norm_face = face_normal.norm();
    double norm_dir = DirVector.norm();
    if (norm_face == 0 || norm_dir == 0)
    {
        std::cerr << "Warning: Zero-length normal or direction vector." << std::endl;
        return false; // Invalid normals
    }

    double dot = face_normal.dot(DirVector) / (norm_face * norm_dir);

    // Clamp dot to valid range to avoid numerical issues
    dot = std::max(-1.0, std::min(1.0, dot));

    // Compute angle in radians
    double angle_rad = std::acos(dot);

    // Define a forward cone: angle less than 90 - threshold
    return (dot > std::cos((90.0 - angle_threshold_deg) * M_PI / 180.0));
}

void RayTracerManager::generateTrackingData(int rays_per_face)
{
    // Retrieve boundary faces
    std::vector<Face> boundary_faces = mesh_.getBoundaryFaces();
    if (boundary_faces.empty())
    {
        std::cerr << "Warning: No boundary faces found." << std::endl;
        return;
    }

    // Determine the number of RayTracer instances
    int num_ray_tracers = ray_tracers_.size();
    if(num_ray_tracers == 0) {
        std::cerr << "Error: No RayTracer instances initialized." << std::endl;
        return;
    }

    // Determine the number of threads
    int num_threads = omp_get_max_threads();

    // Precompute total expected rays for memory reservation
    size_t total_rays = boundary_faces.size() * rays_per_face * num_ray_tracers;
    tracking_data_.reserve(total_rays);

    // Initialize thread-local tracking data
    std::vector<std::vector<TrackingData>> thread_local_data(num_threads);
    for(auto &local_data : thread_local_data)
    {
        local_data.reserve(rays_per_face * num_ray_tracers * (boundary_faces.size() / num_threads + 1));
    }

    // Parallelize over boundary faces
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::vector<TrackingData> &local_tracking_data = thread_local_data[thread_id];

        #pragma omp for schedule(dynamic)
        for(int i = 0; i < static_cast<int>(boundary_faces.size()); ++i)
        {
            const Face &face = boundary_faces[i];

            // Validate face node indices
            if (face.n0 >= mesh_.getNodes().size() ||
                face.n1 >= mesh_.getNodes().size() ||
                face.n2 >= mesh_.getNodes().size())
            {
                std::cerr << "Warning: Face node index out of bounds. Face ID: " << i << std::endl;
                continue; // Skip invalid face
            }

            std::array<Vector3D, 3> face_nodes = {
                mesh_.getNodes()[face.n0],
                mesh_.getNodes()[face.n1],
                mesh_.getNodes()[face.n2]};

            // Compute face normal
            Vector3D face_normal = computeFaceNormal(face_nodes);
            double norm = face_normal.norm();
            if (norm == 0)
            {
                std::cerr << "Warning: Degenerate face encountered. Face ID: " << i << std::endl;
                continue; // Skip degenerate face
            }

            // Determine the starting cell ID (only one adjacent cell for boundary faces)
            int start_cell_id = mesh_.getFaceAdjacentCell(face, /*only_one=*/true);
            if (start_cell_id == -1)
            {
                std::cerr << "Warning: Unable to determine starting cell for face ID: " << i << std::endl;
                continue; // Invalid face
            }

            // Sample starting points on the face
            for (int ray = 0; ray < rays_per_face; ++ray)
            {
                Vector3D start_point = samplePointOnTriangle(face_nodes);

                // Iterate over RayTracer instances
                for(int d = 0; d < num_ray_tracers; ++d)
                {
                    const RayTracer& ray_tracer = *(ray_tracers_[d]);

                    // Determine direction based on RayTracer mode
                    Vector3D direction;
                    if(ray_tracer.getMode() == RayTracerMode::VARIABLE_DIRECTION) {
                        direction = ray_tracer.getField().getVectorFields()[start_cell_id];
                    } else { // CONSTANT_DIRECTION
                        direction = ray_tracer.getFixedDirection();
                    }

                    // Check if direction is valid (not parallel and outward)
                    // For VARIABLE_DIRECTION mode, retrieve corresponding Direction struct
                    // For CONSTANT_DIRECTION mode, direction is already defined
                    bool is_valid = false;
                    if(ray_tracer.getMode() == RayTracerMode::VARIABLE_DIRECTION) {
                        // Assuming the order in ray_tracers_ matches the directions used to initialize them
                        // Otherwise, store Direction information within RayTracer
                        const std::vector<Direction>& directions = angular_quadrature_.getDirections();
                        int direction_index = d < (directions.size() / 2) ? d : d - (directions.size() / 2);
                        if(direction_index < directions.size()) {
                            is_valid = isValidDirection(face_normal, directions[direction_index], /*threshold=*/1.0);
                        }
                    } else { // CONSTANT_DIRECTION
                        // Compute the angle between face_normal and direction
                        double dot = face_normal.normalized().dot(direction.normalized());
                        is_valid = (dot > std::cos((90.0 - 1.0) * M_PI / 180.0)); // 1 degree threshold
                    }

                    if(is_valid)
                    {
                        // Trace the ray using the corresponding RayTracer
                        std::vector<CellTrace> cell_traces = ray_tracer.traceRay(start_cell_id, start_point, 100);
                        if (cell_traces.empty())
                        {
                            std::cerr << "Warning: Ray tracing returned empty cell traces. Ray ID: " << i << std::endl;
                            continue; // Skip rays with no traversal
                        }

                        // Populate TrackingData
                        TrackingData data;
                        // Assign a unique ray_id based on thread and local counter
                        data.ray_id = thread_id * (rays_per_face * num_ray_tracers) + ray * num_ray_tracers + d;
                        data.direction = (ray_tracer.getMode() == RayTracerMode::VARIABLE_DIRECTION) ?
                                         direction :
                                         direction; // Both modes assign the direction used
                        data.cell_traces = std::move(cell_traces);

                        // Add to local tracking data
                        local_tracking_data.push_back(std::move(data));
                    }
                }
            }
        }
    }

    // Merge all thread-local data into the main tracking_data_ vector
    for(auto &local_data : thread_local_data)
    {
        tracking_data_.insert(tracking_data_.end(),
                              std::make_move_iterator(local_data.begin()),
                              std::make_move_iterator(local_data.end()));
    }
}