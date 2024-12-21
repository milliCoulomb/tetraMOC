// include/TrackingData.hpp
#ifndef TRACKINGDATA_HPP
#define TRACKINGDATA_HPP

#include <vector>
#include "Vector3D.hpp"

// Structure to store information about a single cell traversal
struct CellTrace {
    int cell_id;
    double time_spent;
    Vector3D start_point;
    Vector3D end_point;
};

// Structure to store tracking data for a single ray
struct TrackingData {
    int ray_id;
    Vector3D direction; // Added direction for better traceability
    double direction_weight; // Added direction weight for better traceability
    std::vector<CellTrace> cell_traces;
};

#endif // TRACKINGDATA_HPP