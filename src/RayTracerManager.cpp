// src/RayTracerManager.cpp

#include "RayTracerManager.hpp"
#include "GeometryUtils.hpp"

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <atomic>
#include <iostream> // For debugging purposes
#include <vector>
#include <array>
#include <mutex>

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
        double direction_weight = dir.weight;

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
        ray_tracers_.emplace_back(std::make_unique<RayTracer>(mesh_, normalized_dir, direction_weight));
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

void RayTracerManager::generateTrackingData(int rays_per_face, int max_ray_length)
{
    // Clear previous tracking data
    tracking_data_.clear();

    // Retrieve all boundary faces
    const auto& boundary_faces = mesh_.getBoundaryFaces();
    Logger::info("Boundary faces: " + std::to_string(boundary_faces.size()));

    // Showcase direction of all RayTracers
    for (size_t i = 0; i < ray_tracers_.size(); ++i)
    {
        RayTracer* tracer = ray_tracers_[i].get();
        RayTracerMode mode = tracer->getMode();
        Logger::info("RayTracer " + std::to_string(i) + " Mode: " + (mode == RayTracerMode::VARIABLE_DIRECTION ? "VARIABLE_DIRECTION" : "CONSTANT_DIRECTION"));
        // std::cout << "RayTracer " << i << " Mode: " << (mode == RayTracerMode::VARIABLE_DIRECTION ? "VARIABLE_DIRECTION" : "CONSTANT_DIRECTION") << std::endl;
        if(mode == RayTracerMode::CONSTANT_DIRECTION) {
            Vector3D dir = tracer->getFixedDirection();
            Logger::info("Direction: (" + std::to_string(dir.x) + ", " + std::to_string(dir.y) + ", " + std::to_string(dir.z) + ")");
        }
        else if(mode == RayTracerMode::VARIABLE_DIRECTION) {
            // std::cout << "Direction: (Variable Direction)" << std::endl;
            Logger::info("Direction: (Variable Direction)");
        }
    }

    // Determine the number of threads available
    int num_threads = omp_get_max_threads();

    // Initialize per-thread tracking data containers
    std::vector<std::vector<TrackingData>> thread_tracking_data(num_threads);

    // Initialize a vector to hold per-thread ray ID offsets
    std::vector<int> thread_ray_id_offsets(num_threads, 0);

    // First parallel region: Generate tracking data without ray IDs
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        auto& local_tracking_data = thread_tracking_data[thread_id];
        local_tracking_data.reserve(rays_per_face * boundary_faces.size() / ray_tracers_.size());

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
                Vector3D cell_center = mesh_.getCellCenter(adjacent_cell_id);

                // Compute the face normal using GeometryUtils with cell center
                std::array<Vector3D, 3> triangle = { v0, v1, v2 };
                Vector3D face_normal = computeFaceNormal(triangle, cell_center);

                // Determine the direction based on RayTracer mode
                Vector3D direction;
                double direction_weight = 1.0; // Default weight
                if (mode == RayTracerMode::VARIABLE_DIRECTION)
                {
                    // Retrieve the direction from the Field based on the adjacent cell
                    Vector3D field_direction = base_field_.getVectorFields()[adjacent_cell_id];
                    if (field_direction.isAlmostEqual(ZERO_VECTOR))
                    {
                        // Skip rays with zero direction
                        Logger::warning("Skipping RayTracer " + std::to_string(i) + " for face due to zero direction.");
                        continue;
                    }
                    direction = field_direction.normalized();
                }
                else if (mode == RayTracerMode::CONSTANT_DIRECTION)
                {
                    Vector3D fixed_direction = tracer->getFixedDirection();
                    direction_weight = tracer->getDirectionWeight();
                    if (fixed_direction.isAlmostEqual(ZERO_VECTOR))
                    {
                        // Skip rays with zero direction
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
                if (isValidDirection(face_normal, direction, 1e-5))
                {
                    // For each ray per face
                    for (int ray = 0; ray < rays_per_face; ++ray)
                    {
                        // Sample a starting point on the face
                        Vector3D start_point = samplePointOnTriangle(triangle);

                        // Trace the ray through the mesh
                        std::vector<CellTrace> cell_traces = tracer->traceRay(adjacent_cell_id, start_point, max_ray_length);

                        // If cell_traces is empty, discard the ray
                        if (cell_traces.empty()) {
                            continue;
                        }

                        // Populate TrackingData without ray_id for now
                        TrackingData data;
                        data.direction = direction;
                        data.direction_weight = direction_weight;
                        data.cell_traces = std::move(cell_traces);

                        local_tracking_data.emplace_back(std::move(data));
                    }
                }
            }
        }
    }

    // Compute ray ID offsets
    size_t total_rays = 0;
    for (int i = 0; i < num_threads; ++i)
    {
        thread_ray_id_offsets[i] = static_cast<int>(total_rays);
        total_rays += thread_tracking_data[i].size();
    }

    // Assign unique ray IDs
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int ray_offset = thread_ray_id_offsets[thread_id];
        auto& local_tracking_data = thread_tracking_data[thread_id];

        #pragma omp for nowait
        for (size_t i = 0; i < local_tracking_data.size(); ++i)
        {
            local_tracking_data[i].ray_id = ray_offset + static_cast<int>(i);
        }
    }

    // Merge all thread-local tracking data into the global tracking_data_
    // Reserve space to avoid multiple reallocations
    tracking_data_.reserve(total_rays);
    for (const auto& local_data : thread_tracking_data)
    {
        tracking_data_.insert(tracking_data_.end(), local_data.begin(), local_data.end());
    }

    Logger::info("Generated " + std::to_string(total_rays) + " tracking data entries.");
}

void RayTracerManager::doubleTrackingDataByReversing() {
    std::vector<TrackingData> reversed_tracking_data;
    reversed_tracking_data.reserve(tracking_data_.size());

    // Determine the current maximum ray_id to ensure uniqueness
    int max_ray_id = 0;
    for (const auto& ray : tracking_data_) {
        if (ray.ray_id > max_ray_id) {
            max_ray_id = ray.ray_id;
        }
    }

    // Start assigning new ray_ids from max_ray_id + 1
    int new_ray_id = max_ray_id + 1;

    for (const auto& original_ray : tracking_data_) {
        TrackingData reversed_ray = original_ray; // Copy existing data

        // Assign a new unique ray_id
        reversed_ray.ray_id = new_ray_id++;

        // Reverse the direction
        reversed_ray.direction = reversed_ray.direction * -1.0;

        // Reverse the order of cell_traces
        std::reverse(reversed_ray.cell_traces.begin(), reversed_ray.cell_traces.end());

        // // Swap start_point and end_point in each CellTrace
        // for (auto& trace : reversed_ray.cell_traces) {
        //     std::swap(trace.start_point, trace.end_point);
        // }

        // Optionally, adjust direction_weight if necessary
        // For this example, we'll keep it the same
        // reversed_ray.direction_weight = reversed_ray.direction_weight;

        reversed_tracking_data.push_back(reversed_ray);
    }

    // Append the reversed tracking data to the original tracking_data_
    tracking_data_.insert(tracking_data_.end(),
                          reversed_tracking_data.begin(),
                          reversed_tracking_data.end());

    Logger::info("Tracking data doubled by reversing each ray.");
}

// void RayTracerManager::symmetrizeTrackingData()
// {
//     std::vector<TrackingData> symmetrized_data;
//     symmetrized_data.reserve(tracking_data_.size() * 2);

//     // Atomic counter for unique ray IDs
//     std::atomic<int> new_ray_id_counter(tracking_data_.size());

//     for (const auto& data : tracking_data_)
//     {
//         // Original ray
//         symmetrized_data.push_back(data);

//         // Create reversed ray
//         TrackingData reversed_data;
//         reversed_data.ray_id = new_ray_id_counter.fetch_add(1);
//         reversed_data.direction = -data.direction;
//         reversed_data.cell_traces = data.cell_traces;
//         std::reverse(reversed_data.cell_traces.begin(), reversed_data.cell_traces.end());

//         // Reverse time_spent if it exists
//         if (!data.time_spent.empty())
//         {
//             reversed_data.time_spent = data.time_spent;
//             std::reverse(reversed_data.time_spent.begin(), reversed_data.time_spent.end());
//         }

//         // Swap start and end points
//         reversed_data.start_point = data.end_point;
//         reversed_data.end_point = data.start_point;

//         symmetrized_data.push_back(reversed_data);
//     }

//     tracking_data_ = std::move(symmetrized_data);
// }

