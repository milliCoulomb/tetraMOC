// include/Tetrahedron.hpp

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "MeshHandler.hpp"
#include "Vector3D.hpp" // Use your existing Vector3D class
#include "Field.hpp" // Include the Field class
#include <array>

struct DirectionData {
    int cell_id;
    double time_spent;
    std::array<double, 3> start_point;
    std::array<double, 3> end_point;
};

class Tetrahedron {
public:
    // Updated constructor to accept CellVectorField
    Tetrahedron(const TetraCell& cell, const std::vector<Node>& nodes, const CellVectorField& field);

    // Use your project's Vector3D class instead of SNSolver::Vector3D
    bool findExit(const std::array<double, 3>& x0, const Vector3D& v, double& t_exit, std::array<double, 3>& x_exit, int& exit_face_id) const;

private:
    std::array<Vector3D, 4> vertices; // Use your project's Vector3D
    Vector3D velocity;
};

#endif // TETRAHEDRON_H