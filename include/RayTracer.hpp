// include/ray_tracer.h
#ifndef RAY_TRACER_H
#define RAY_TRACER_H

#include "MeshHandler.hpp"
#include "Field.hpp"
#include "Tetrahedron.hpp"
#include "Vector3D.hpp" // Your existing Vector3D class
#include <vector>
#include <array>
#include <unordered_map>
#include <tuple>

// Custom hash function for tuple<int, int, int>
struct TupleHash {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        std::size_t seed = 0;
        std::hash<int> hasher;
        seed ^= hasher(std::get<0>(t)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hasher(std::get<1>(t)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hasher(std::get<2>(t)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }
};

class RayTracer {
public:
    RayTracer(const MeshHandler& mesh, const Field& field);
    bool loadFaceConnectivity(const std::string& filename);
    std::vector<DirectionData> traceRay(int start_cell_id, const std::array<double, 3>& start_point, int max_iter = 100);
    
private:
    const MeshHandler& mesh;
    const Field& field;
    
    // Map of sorted face (n0, n1, n2) to adjacent cell IDs
    std::unordered_map<std::tuple<int, int, int>, std::vector<int>, TupleHash> face_to_cells;
    
    // Find the neighboring cell given current cell and exit face
    int getNeighborCell(int current_cell_id, int exit_face_id) const;
};

#endif // RAY_TRACER_H