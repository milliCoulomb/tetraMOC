// include/RayTracerManager.hpp
#ifndef RAYTRACERMANAGER_HPP
#define RAYTRACERMANAGER_HPP

#include "RayTracer.hpp"
#include "AngularQuadrature.hpp"
#include "TrackingData.hpp"
#include <vector>
#include <array>
#include <mutex>
#include <atomic>

class RayTracerManager {
public:
    RayTracerManager(const MeshHandler& mesh,
                    const Field& field,
                    RayTracer& tracer,
                    AngularQuadrature& angular_quadrature);

    // Generate tracking data with a specified number of rays per boundary face
    void generateTrackingData(int rays_per_face);

    // Retrieve the tracking data
    const std::vector<TrackingData>& getTrackingData() const { return tracking_data_; }

private:
    const MeshHandler& mesh_;
    const Field& field_;
    RayTracer& tracer_;
    AngularQuadrature& angular_quadrature_;

    std::vector<TrackingData> tracking_data_;
    std::mutex data_mutex_;
    std::atomic<int> ray_counter_{0};

    // Function to determine if a direction is valid (not parallel)
    bool isValidDirection(const std::array<double, 3>& face_normal,
                          const Direction& dir,
                          double angle_threshold_deg) const;
};
#endif // RAYTRACERMANAGER_HPP