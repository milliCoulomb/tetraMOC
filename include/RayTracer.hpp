// include/RayTracer.hpp
#ifndef RAY_TRACER_HPP
#define RAY_TRACER_HPP

#include "MeshHandler.hpp"
#include "Field.hpp"
#include "Tetrahedron.hpp"
#include "Vector3D.hpp" // Ensure this defines Vector3D
#include <TrackingData.hpp>
#include <vector>
#include <array>
#include <unordered_map>
#include <tuple>


class RayTracer {
public:
    RayTracer(const MeshHandler& mesh, const Field& field);
    std::vector<CellTrace> traceRay(int start_cell_id, const Vector3D& start_point, int max_iter = 100) const;
    
private:
    const MeshHandler& mesh;
    const Field& field;
};

#endif // RAY_TRACER_HPP