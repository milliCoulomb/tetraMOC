// include/TrackingData.hpp
#ifndef TRACKINGDATA_HPP
#define TRACKINGDATA_HPP

#include "Vector3D.hpp"

struct CellTrace {
    int cell_id;
    double time_spent;
    Vector3D start_point;
    Vector3D end_point;
};

struct TrackingData {
    int ray_id = -1; // Initialize ray_id to a default value
    Vector3D direction;
    double direction_weight;
    std::vector<CellTrace> cell_traces;

    // Default constructor
    TrackingData() = default;

    // Copy constructor
    TrackingData(const TrackingData& other) = default;

    // Move constructor
    TrackingData(TrackingData&& other) noexcept = default;

    // Copy assignment operator
    TrackingData& operator=(const TrackingData& other) = default;

    // Move assignment operator
    TrackingData& operator=(TrackingData&& other) noexcept = default;
};

#endif // TRACKINGDATA_HPP