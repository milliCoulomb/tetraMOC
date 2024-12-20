// src/RayTracerManager.cpp

#include "RayTracerManager.hpp"
#include "GeometryUtils.hpp"

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <atomic>

// Constructor for managing variable direction RayTracers
RayTracerManager::RayTracerManager(const MeshHandler& mesh,
                                   const Field& base_field,
                                   AngularQuadrature& angular_quadrature)
    : mesh_(mesh),
      base_field_(base_field),
      angular_quadrature_(angular_quadrature)
{
    initializeRayTracers();
}

// Overloaded Constructor for managing both variable and constant direction RayTracers
// The boolean flag indicates whether to use half of the quadrature for constant directions
RayTracerManager::RayTracerManager(const MeshHandler& mesh,
                                   const Field& base_field,
                                   AngularQuadrature& angular_quadrature,
                                   bool use_half_quadrature_for_constant)
    : mesh_(mesh),
      base_field_(base_field),
      angular_quadrature_(angular_quadrature)
{
    initializeRayTracers();

    if (use_half_quadrature_for_constant)
    {
        // Get all directions from AngularQuadrature
        const std::vector<Direction>& all_directions = angular_quadrature_.getDirections();

        // Initialize constant direction RayTracers using half of the directions (e.g., mu >= 0)
        initializeConstantDirectionRayTracers(all_directions);
    }
}

// Helper method to initialize RayTracers for variable directions
void RayTracerManager::initializeRayTracers()
{
    // Create a single RayTracer in VARIABLE_DIRECTION mode
    ray_tracers_.emplace_back(std::make_unique<RayTracer>(mesh_, base_field_));
}

// Helper method to initialize RayTracers with constant directions based on hemisphere
void RayTracerManager::initializeConstantDirectionRayTracers(const std::vector<Direction>& quadrature_directions)
{
    // Reserve additional space to avoid reallocations
    ray_tracers_.reserve(ray_tracers_.size() + (quadrature_directions.size() / 2));

    for (const auto& dir : quadrature_directions)
    {
        // Select directions from one hemisphere (e.g., mu >= 0)
        if (dir.mu >= 0.0)
        {
            // Convert Direction to Vector3D
            Vector3D vector_dir(
                std::sqrt(1.0 - dir.mu * dir.mu) * std::cos(dir.phi),
                std::sqrt(1.0 - dir.mu * dir.mu) * std::sin(dir.phi),
                dir.mu
            );

            // Normalize the direction vector
            Vector3D normalized_dir = vector_dir.normalized();

            // Instantiate a RayTracer in CONSTANT_DIRECTION mode
            ray_tracers_.emplace_back(std::make_unique<RayTracer>(mesh_, normalized_dir));
        }
    }
}

// Method to check if a direction is valid (incoming) for a given face normal
bool RayTracerManager::isValidDirection(const Vector3D& face_normal, const Vector3D& direction, double threshold) const
{
    // Compute the dot product between face normal and direction
    double dot_product = face_normal.dot(direction);

    // Check if the dot product is less than -threshold (incoming ray)
    return dot_product < -threshold;
}

// Method to generate tracking data by tracing rays
void RayTracerManager::generateTrackingData(int rays_per_face)
{
    // Clear previous tracking data
    tracking_data_.clear();

    // Retrieve all boundary faces
    const auto& boundary_faces = mesh_.getBoundaryFaces();

    // Initialize an atomic counter for unique ray IDs
    std::atomic<int> ray_id_counter(0);

    // Parallelize over RayTracers to utilize multiple directions
    #pragma omp parallel
    {
        std::vector<TrackingData> local_tracking_data;
        // Optionally reserve space to improve performance
        // local_tracking_data.reserve(boundary_faces.size() * rays_per_face);

        // Iterate over each RayTracer
        #pragma omp for nowait
        for (size_t i = 0; i < ray_tracers_.size(); ++i)
        {
            RayTracer* tracer = ray_tracers_[i].get();
            RayTracerMode mode = tracer->getMode();

            // Iterate over each boundary face
            for (const auto& face : boundary_faces)
            {
                // Retrieve the node indices for the face
                int n0 = face.n0;
                int n1 = face.n1;
                int n2 = face.n2;

                // Retrieve the node coordinates from MeshHandler
                const Vector3D& v0 = mesh_.getNodes()[n0];
                const Vector3D& v1 = mesh_.getNodes()[n1];
                const Vector3D& v2 = mesh_.getNodes()[n2];

                // Compute the face normal using GeometryUtils
                std::array<Vector3D, 3> triangle = { v0, v1, v2 };
                Vector3D face_normal = computeFaceNormal(triangle);

                // Determine the direction based on RayTracer mode
                Vector3D direction;
                if (mode == RayTracerMode::VARIABLE_DIRECTION)
                {
                    // Retrieve the direction from the Field based on the adjacent cell
                    int adjacent_cell_id = mesh_.getFaceAdjacentCell(face, true); // true for boundary face
                    Vector3D field_direction = base_field_.getVectorFields()[adjacent_cell_id];
                    if (field_direction.x == 0.0 && field_direction.y == 0.0 && field_direction.z == 0.0)
                    {
                        // Skip rays with zero direction
                        continue;
                    }
                    direction = field_direction.normalized();
                }
                else if (mode == RayTracerMode::CONSTANT_DIRECTION)
                {
                    Vector3D fixed_direction = tracer->getFixedDirection();
                    if (fixed_direction.x == 0.0 && fixed_direction.y == 0.0 && fixed_direction.z == 0.0)
                    {
                        // Skip rays with zero direction
                        continue;
                    }
                    direction = fixed_direction.normalized();
                }
                else
                {
                    continue; // Unknown mode
                }

                // Check if the direction is valid (incoming)
                if (isValidDirection(face_normal, direction, 1e-6))
                {
                    // For each ray per face
                    for (int ray = 0; ray < rays_per_face; ++ray)
                    {
                        // Sample a starting point on the face
                        Vector3D start_point = samplePointOnTriangle(triangle);

                        // Get the adjacent cell ID (assuming rays start from the boundary face into the adjacent cell)
                        int start_cell_id = mesh_.getFaceAdjacentCell(face, true); // true for boundary face

                        // Trace the ray through the mesh
                        std::vector<CellTrace> cell_traces = tracer->traceRay(start_cell_id, start_point, 100);

                        // Assign a unique ray_id
                        int ray_id = ray_id_counter.fetch_add(1, std::memory_order_relaxed);

                        // Populate TrackingData
                        TrackingData data;
                        data.ray_id = ray_id;
                        data.direction = direction;

                        // Populate cell_traces
                        data.cell_traces = cell_traces;

                        local_tracking_data.push_back(data);
                    }
                }
            }
        }

        // Protect shared resource with a critical section
        #pragma omp critical
        {
            tracking_data_.insert(tracking_data_.end(),
                                  local_tracking_data.begin(),
                                  local_tracking_data.end());
        }
    }
}

// Getter for tracking data (already provided in the header)
// const std::vector<TrackingData>& RayTracerManager::getTrackingData() const { return tracking_data_; }