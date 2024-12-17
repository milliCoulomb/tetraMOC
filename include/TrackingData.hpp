// include/TrackingData.hpp
#ifndef TRACKINGDATA_HPP
#define TRACKINGDATA_HPP

#include <vector>
#include <array>
// Structure to store information about a single cell traversal
struct CellTrace {
    int cell_id;
    double time_spent;
    std::array<double, 3> start_point;
    std::array<double, 3> end_point;
};

// Structure to store tracking data for a single ray
struct TrackingData {
    int ray_id;
    std::array<double, 3> direction;
    std::vector<CellTrace> cell_traces;
};

#endif // TRACKINGDATA_HPP