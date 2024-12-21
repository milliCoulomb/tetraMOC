// src/RayTracerManager.cpp

#include "RayTracerManager.hpp"
#include "GeometryUtils.hpp"

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <atomic>
#include <iostream> // For debugging purposes

const Vector3D ZERO_VECTOR(0.0, 0.0, 0.0);

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

// Overloaded Constructor for managing constant direction RayTracers only
RayTracerManager::RayTracerManager(const MeshHandler& mesh,
                                   const Field& base_field,
                                   AngularQuadrature& angular_quadrature,
                                   bool constant_directions,
                                   bool use_half_quadrature_for_constant)
    : mesh_(mesh),
      base_field_(base_field),
      angular_quadrature_(angular_quadrature)
{
    if (!constant_directions)
    {
        // If not using constant directions, initialize variable direction RayTracer
        initializeRayTracers();
    }
    else
    {
        if (use_half_quadrature_for_constant)
        {
            // Get all directions from AngularQuadrature
            const std::vector<Direction>& all_directions = angular_quadrature_.getDirections();

            // Initialize constant direction RayTracers using half of the directions
            initializeConstantDirectionRayTracers(all_directions, use_half_quadrature_for_constant);
        }
        else
        {
            // if not, we use all directions
            const std::vector<Direction>& all_directions = angular_quadrature_.getDirections();
            initializeConstantDirectionRayTracers(all_directions, use_half_quadrature_for_constant);
        }
    }
    // if (use_half_quadrature_for_constant)
    // {
    //     // Get all directions from AngularQuadrature
    //     const std::vector<Direction>& all_directions = angular_quadrature_.getDirections();

    //     // Initialize constant direction RayTracers using half of the directions
    //     initializeConstantDirectionRayTracers(all_directions);
    // }
    // else
    // {
    //     // If not using half quadrature for constant directions, initialize variable direction RayTracer
    //     initializeRayTracers();
    // }
}

// Helper method to initialize RayTracers for variable directions
void RayTracerManager::initializeRayTracers()
{
    // Create a single RayTracer in VARIABLE_DIRECTION mode
    ray_tracers_.emplace_back(std::make_unique<RayTracer>(mesh_, base_field_));
}

// Helper method to initialize RayTracers with constant directions
void RayTracerManager::initializeConstantDirectionRayTracers(const std::vector<Direction>& quadrature_directions, bool use_half_quadrature_for_constant)
{
    // check if we use half or full quadrature
    size_t half_size = quadrature_directions.size();
    if (use_half_quadrature_for_constant)
    {
        half_size = quadrature_directions.size() / 2;
    }
    // Determine half of the quadrature directions
    size_t added_tracers = 0;

    for (size_t i = 0; i < half_size; ++i)
    {
        const auto& dir = quadrature_directions[i];

        // Convert Direction to Vector3D
        double sqrt_term = std::sqrt(1.0 - dir.mu * dir.mu);
        double x = sqrt_term * std::cos(dir.phi);
        double y = sqrt_term * std::sin(dir.phi);
        double z = dir.mu;

        Vector3D vector_dir(x, y, z);

        // Check if vector_dir is zero
        if (vector_dir.isAlmostEqual(ZERO_VECTOR))
        {
            // std::cerr << "Warning: Zero direction vector encountered. Skipping this direction." << std::endl;
            Logger::warning("Zero direction vector encountered. Skipping this direction.");
            continue; // Skip adding this RayTracer
        }

        // Normalize the direction vector
        Vector3D normalized_dir = vector_dir.normalized();

        // Instantiate a RayTracer in CONSTANT_DIRECTION mode
        ray_tracers_.emplace_back(std::make_unique<RayTracer>(mesh_, normalized_dir));
        added_tracers++;
    }
    Logger::info("Initialized " + std::to_string(added_tracers) + " constant direction RayTracers.");
    // std::cout << "Initialized " << added_tracers << " constant direction RayTracers." << std::endl;
}

// Method to check if a direction is valid (incoming) for a given face normal
bool RayTracerManager::isValidDirection(const Vector3D& face_normal, const Vector3D& direction, double threshold) const
{
    // Compute the dot product between face normal and direction
    double dot_product = face_normal.dot(direction);

    // Check if the dot product is less than -threshold (incoming ray)
    return dot_product < -threshold;
}

void RayTracerManager::generateTrackingData(int rays_per_face)
{
    // Set number of threads to 1 for debugging
    omp_set_num_threads(1);

    // Clear previous tracking data
    tracking_data_.clear();
    // std::cout << "Tracking data cleared." << std::endl;

    // Retrieve all boundary faces
    const auto& boundary_faces = mesh_.getBoundaryFaces();
    // std::cout << "Boundary faces: " << boundary_faces.size() << std::endl;
    // std::cout << "Number of RayTracers: " << ray_tracers_.size() << std::endl;

    // Showcase direction of all RayTracers
    for (size_t i = 0; i < ray_tracers_.size(); ++i)
    {
        RayTracer* tracer = ray_tracers_[i].get();
        RayTracerMode mode = tracer->getMode();
        Logger::info("RayTracer" + std::to_string(i) + " Mode: " + (mode == RayTracerMode::VARIABLE_DIRECTION ? "VARIABLE_DIRECTION" : "CONSTANT_DIRECTION"));
        // std::cout << "RayTracer " << i << " Mode: " << (mode == RayTracerMode::VARIABLE_DIRECTION ? "VARIABLE_DIRECTION" : "CONSTANT_DIRECTION") << std::endl;
        if(mode == RayTracerMode::CONSTANT_DIRECTION) {
            Vector3D dir = tracer->getFixedDirection();
            Logger::info("Direction: (" + std::to_string(dir.x) + ", " + std::to_string(dir.y) + ", " + std::to_string(dir.z) + ")");
            // std::cout << "Direction: (" << dir.x << ", " << dir.y << ", " << dir.z << ")" << std::endl;
        }
        else if(mode == RayTracerMode::VARIABLE_DIRECTION) {
            // std::cout << "Direction: (Variable Direction)" << std::endl;
            Logger::info("Direction: (Variable Direction)");
        }
    }

    // Initialize an atomic counter for unique ray IDs
    std::atomic<int> ray_id_counter(0);

    // Parallelize over RayTracers to utilize multiple directions
    #pragma omp parallel
    {
        std::vector<TrackingData> local_tracking_data;

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

                // Compute the cell center
                int adjacent_cell_id = mesh_.getFaceAdjacentCell(face, true); // true for boundary face
                // std::cout << "Adjacent Cell ID: " << adjacent_cell_id << std::endl;
                Vector3D cell_center = mesh_.getCellCenter(adjacent_cell_id);

                // Compute the face normal using GeometryUtils with cell center
                std::array<Vector3D, 3> triangle = { v0, v1, v2 };
                Vector3D face_normal = computeFaceNormal(triangle, cell_center);

                // Determine the direction based on RayTracer mode
                Vector3D direction;
                if (mode == RayTracerMode::VARIABLE_DIRECTION)
                {
                    // Retrieve the direction from the Field based on the adjacent cell
                    Vector3D field_direction = base_field_.getVectorFields()[adjacent_cell_id];
                    if (field_direction.isAlmostEqual(ZERO_VECTOR))
                    {
                        // Skip rays with zero direction
                        // std::cout << "Skipping RayTracer " << i << " for face due to zero direction." << std::endl;
                        Logger::warning("Skipping RayTracer " + std::to_string(i) + " for face due to zero direction.");
                        continue;
                    }
                    direction = field_direction.normalized();
                }
                else if (mode == RayTracerMode::CONSTANT_DIRECTION)
                {
                    Vector3D fixed_direction = tracer->getFixedDirection();
                    if (fixed_direction.isAlmostEqual(ZERO_VECTOR))
                    {
                        // Skip rays with zero direction
                        // std::cout << "Skipping RayTracer " << i << " for face due to zero direction." << std::endl;
                        Logger::warning("Skipping RayTracer " + std::to_string(i) + " for face due to zero direction.");
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
                        // Already have adjacent_cell_id

                        // Trace the ray through the mesh
                        std::vector<CellTrace> cell_traces = tracer->traceRay(adjacent_cell_id, start_point, 100);

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
