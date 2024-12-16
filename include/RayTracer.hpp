// include/RayTracer.hpp
#ifndef RAY_TRACER_HPP
#define RAY_TRACER_HPP

#include "MeshHandler.hpp"
#include "Field.hpp"
#include "Tetrahedron.hpp"
#include "Vector3D.hpp" // Ensure this defines SNSolver::Vector3D
#include <vector>
#include <array>
#include <unordered_map>
#include <tuple>


class RayTracer {
public:
    RayTracer(const MeshHandler& mesh, const Field& field);
    std::vector<DirectionData> traceRay(int start_cell_id, const std::array<double, 3>& start_point, int max_iter = 100) const;
    
private:
    const MeshHandler& mesh;
    const Field& field;
};

#endif // RAY_TRACER_HPP