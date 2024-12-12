// include/tetrahedron.h
#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "MeshHandler.h"
#include "Vector3D.hpp" // Include your existing Vector3D class
#include <array>

struct DirectionData {
    int cell_id;
    double time_spent;
    std::array<double, 3> start_point;
    std::array<double, 3> end_point;
};

class Tetrahedron {
public:
    Tetrahedron(const TetraCell& cell, const std::vector<Node>& nodes, const CellField& field);

    // Use SNSolver::Vector3D for the velocity vector
    bool findExit(const std::array<double, 3>& x0, const SNSolver::Vector3D& v, double& t_exit, std::array<double, 3>& x_exit, int& exit_face_id) const;

private:
    std::array<SNSolver::Vector3D, 4> vertices;
    SNSolver::Vector3D velocity;
};

#endif // TETRAHEDRON_H