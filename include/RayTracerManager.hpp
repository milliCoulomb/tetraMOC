// include/RayTracerManager.hpp

#ifndef RAY_TRACER_MANAGER_HPP
#define RAY_TRACER_MANAGER_HPP

#include "MeshHandler.hpp"
#include "Field.hpp"
#include "RayTracer.hpp"
#include "AngularQuadrature.hpp"
#include "TrackingData.hpp"
#include "Vector3D.hpp"

#include <vector>
#include <memory>
#include <array>
#include <omp.h>

class RayTracerManager {
public:
    // Constructor for managing variable direction RayTracers
    RayTracerManager(const MeshHandler& mesh,
                    const Field& base_field, // Pass a base Field
                    AngularQuadrature& angular_quadrature);

    // Overloaded Constructor for managing both variable and constant direction RayTracers
    RayTracerManager(const MeshHandler& mesh,
                    const Field& base_field, // Pass a base Field
                    AngularQuadrature& angular_quadrature,
                    const std::vector<Vector3D>& constant_directions); // Additional constant directions

    void generateTrackingData(int rays_per_face);
    const std::vector<TrackingData>& getTrackingData() const { return tracking_data_; }

private:
    const MeshHandler& mesh_;
    const Field& base_field_;
    AngularQuadrature& angular_quadrature_;

    // Vector of RayTracers, both variable and constant direction
    std::vector<std::unique_ptr<RayTracer>> ray_tracers_;

    // Aggregated tracking data
    std::vector<TrackingData> tracking_data_;

    bool isValidDirection(const Vector3D& face_normal, const Direction& direction, double threshold) const;

    // Helper method to initialize RayTracers
    void initializeRayTracers();

    // Helper method to initialize RayTracers with constant directions
    void initializeConstantDirectionRayTracers(const std::vector<Vector3D>& constant_directions);
};

#endif // RAY_TRACER_MANAGER_HPP