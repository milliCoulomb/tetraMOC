// include/Tetrahedron.hpp

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "MeshHandler.hpp"
#include "Vector3D.hpp" // Ensure this defines Vector3D
#include "Field.hpp"    // Include the Field class
#include <array>

class Tetrahedron {
public:
    // Constructor accepting CellVectorField from Field class
    Tetrahedron(const TetraCell& cell, const std::vector<Vector3D>& nodes, const CellVectorField& field);

    // findExit method using Vector3D
    bool findExit(const Vector3D& x0, const Vector3D& v, double& t_exit, Vector3D& x_exit, int& exit_face_id) const;

    // Getter for vertices
    const std::array<Vector3D, 4>& getVertices() const { return vertices; }

    // Getter for CenterOfMass
    const Vector3D& getCenterOfMass() const { return CenterOfMass; }

    // Getter for velocity
    const Vector3D& getVelocity() const { return velocity; }

private:
    std::array<Vector3D, 4> vertices; // Fully qualified with SNSolver
    Vector3D CenterOfMass;             // Fully qualified with SNSolver
    Vector3D velocity;                // Fully qualified with SNSolver
};

#endif // TETRAHEDRON_H
