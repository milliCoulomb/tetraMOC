// include/Tetrahedron.hpp

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "MeshHandler.hpp"
#include "Vector3D.hpp" // Ensure this defines SNSolver::Vector3D
#include "Field.hpp"    // Include the Field class
#include <array>

struct DirectionData {
    int cell_id;
    double time_spent;
    std::array<double, 3> start_point;
    std::array<double, 3> end_point;
};

class Tetrahedron {
public:
    // Constructor accepting CellVectorField from Field class
    Tetrahedron(const TetraCell& cell, const std::vector<Node>& nodes, const CellVectorField& field);

    // findExit method using SNSolver::Vector3D
    bool findExit(const std::array<double, 3>& x0, const SNSolver::Vector3D& v, double& t_exit, std::array<double, 3>& x_exit, int& exit_face_id) const;

private:
    std::array<SNSolver::Vector3D, 4> vertices; // Fully qualified with SNSolver
    SNSolver::Vector3D velocity;                // Fully qualified with SNSolver
};

#endif // TETRAHEDRON_H
