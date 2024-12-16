// src/RayTracerManager.cpp
#include "RayTracerManager.hpp"
#include "GeometryUtils.hpp"
#include <cmath>
#include <omp.h>

namespace SNSolver {

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

    // Compute dot product with face normal
    double dot = face_normal[0] * x + face_normal[1] * y + face_normal[2] * z;

    // Compute angle in degrees
    double angle_deg = std::acos(dot) * (180.0 / M_PI);

    // Valid if angle is greater than threshold and less than (180 - threshold)
    return (angle_deg > angle_threshold_deg) && (angle_deg < (180.0 - angle_threshold_deg));
}

void RayTracerManager::generateTrackingData(int rays_per_face) {
    // Retrieve boundary faces
    std::vector<Face> boundary_faces = tracer_.getBoundaryFaces();

    // Retrieve directions from angular quadrature
    const std::vector<Direction>& directions = angular_quadrature_.getDirections();

    // Utilize only half of the directions due to symmetry
    int half_directions = directions.size() / 2;

    // Prepare thread-local storage for tracking data
    int num_threads = omp_get_max_threads();
    std::vector<std::vector<TrackingData>> thread_local_data(num_threads);

    // Parallelize over boundary faces
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::vector<TrackingData>& local_tracking_data = thread_local_data[thread_id];

        #pragma omp for schedule(dynamic)
        for(int i = 0; i < static_cast<int>(boundary_faces.size()); ++i) {
            const Face& face = boundary_faces[i];

            // Retrieve the nodes of the face
            std::array<Node, 3> face_nodes = {
                mesh_.getNodes()[face.n0],
                mesh_.getNodes()[face.n1],
                mesh_.getNodes()[face.n2]
            };

            // Compute face normal
            std::array<double, 3> face_normal = computeFaceNormal(face_nodes);

            // Determine the starting cell ID (only one adjacent cell for boundary faces)
            int start_cell_id = tracer_.getFaceAdjacentCell(face, /*only_one=*/true);
            if(start_cell_id == -1) continue; // Invalid face

            // Sample starting points on the face
            for(int ray = 0; ray < rays_per_face; ++ray) {
                std::array<double, 3> start_point = samplePointOnTriangle(face_nodes);

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

                        // Create a direction vector
                        std::array<double, 3> direction = {x_dir, y_dir, z_dir};

                        // Trace the ray
                        std::vector<CellTrace> cell_traces = tracer_.traceRay(start_cell_id, start_point, 100);

                        // Populate TrackingData
                        TrackingData data;
                        // Assign a unique ray_id using atomic counter
                        data.ray_id = ray_counter_++;
                        data.cell_traces = cell_traces;

                        // Add to local tracking data
                        local_tracking_data.push_back(data);
                    }
                }
            }
        }
    }

    // Merge all thread-local data into the main tracking_data_ vector
    for(const auto& local_data : thread_local_data) {
        tracking_data_.insert(tracking_data_.end(),
                              local_data.begin(),
                              local_data.end());
    }
}

} // namespace SNSolver