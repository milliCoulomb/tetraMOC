// include/Tetrahedron.hpp

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "MeshHandler.hpp"
#include "Vector3D.hpp" // Ensure this defines Vector3D
#include "Field.hpp"    // Include the Field class
#include <array>

// struct DirectionData {
//     int cell_id;
//     double time_spent;
//     std::array<double, 3> start_point;
//     std::array<double, 3> end_point;
// };

class Tetrahedron {
public:
    // Constructor accepting CellVectorField from Field class
    Tetrahedron(const TetraCell& cell, const std::vector<Node>& nodes, const CellVectorField& field);

    // findExit method using Vector3D
    bool findExit(const Vector3D& x0, const Vector3D& v, double& t_exit, Vector3D& x_exit, int& exit_face_id) const;

private:
    std::array<Vector3D, 4> vertices; // Fully qualified with SNSolver
    Vector3D velocity;                // Fully qualified with SNSolver
};

#endif // TETRAHEDRON_H
