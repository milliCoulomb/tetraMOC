// tests/test_tetrahedron.cpp
#include <gtest/gtest.h>
#include "Tetrahedron.hpp"
#include "Vector3D.hpp"

// Helper function to create a simple tetrahedron
TetraCell createSimpleTetraCell(int cell_id, const std::array<int, 4>& node_ids) {
    TetraCell cell;
    cell.cell_id = cell_id;
    cell.node_ids = node_ids;
    return cell;
}

TEST(TetrahedronTest, FindExitBasicTest) {
    // Define nodes
    std::vector<Vector3D> nodes = {
        Vector3D(0.0, 0.0, 0.0), // Node 0
        Vector3D(1.0, 0.0, 0.0), // Node 1
        Vector3D(0.0, 1.0, 0.0), // Node 2
        Vector3D(0.0, 0.0, 1.0)  // Node 3
    };
    
    // Define a single tetrahedron
    TetraCell cell = createSimpleTetraCell(0, {0, 1, 2, 3});
    CellVectorField field = {1.0, 1.0, 1.0}; // Velocity vector
    
    // Create Tetrahedron object
    Tetrahedron tetra(cell, nodes, field);
    
    // Define a starting point inside the tetrahedron
    // std::array<double, 3> start_point = {0.1, 0.1, 0.1};
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Define velocity vector
    Vector3D velocity(1.0, 1.0, 1.0); // Diagonal direction
    
    double t_exit;
    // std::array<double, 3> x_exit;
    Vector3D x_exit;
    int exit_face_id;
    
    // Call findExit
    bool has_exit = tetra.findExit(start_point, velocity, t_exit, x_exit, exit_face_id);
    
    // Assertions
    EXPECT_TRUE(has_exit);
    EXPECT_GT(t_exit, 0.0);
    EXPECT_NE(exit_face_id, -1);
    
    // Check that exit point is within the expected bounds
    EXPECT_GE(x_exit.x, 0.0);
    EXPECT_LE(x_exit.x, 1.0);
    EXPECT_GE(x_exit.y, 0.0);
    EXPECT_LE(x_exit.y, 1.0);
    EXPECT_GE(x_exit.z, 0.0);
    EXPECT_LE(x_exit.z, 1.0);
}

TEST(TetrahedronTest, NoExitTest) {
    // Define nodes
    std::vector<Vector3D> nodes = {
        Vector3D(0.0, 0.0, 0.0), // Node 0
        Vector3D(1.0, 0.0, 0.0), // Node 1
        Vector3D(0.0, 1.0, 0.0), // Node 2
        Vector3D(0.0, 0.0, 1.0)  // Node 3
    };
    
    // Define a single tetrahedron
    TetraCell cell = createSimpleTetraCell(0, {0, 1, 2, 3});
    CellVectorField field = {-1.0, -1.0, -1.0}; // Velocity vector pointing inward
    
    // Create Tetrahedron object
    Tetrahedron tetra(cell, nodes, field);
    
    // Define a starting point inside the tetrahedron
    // std::array<double, 3> start_point = {0.0, 0.0, 0.1};
    Vector3D start_point(0.0, 0.0, 0.1);
    
    // Define velocity vector
    Vector3D velocity(-1.0, -1.0, -1.0); // Inward direction
    
    double t_exit;
    // std::array<double, 3> x_exit;
    Vector3D x_exit;
    int exit_face_id;
    
    // Call findExit
    bool has_exit = tetra.findExit(start_point, velocity, t_exit, x_exit, exit_face_id);
    // Assertions
    EXPECT_FALSE(has_exit);
}

// test Tetrahedron construction and getters for vertices, center of mass, and velocity
TEST(TetrahedronTest, TetrahedronConstructionTest) {
    // Define nodes
    std::vector<Vector3D> nodes = {
        Vector3D(0.0, 0.0, 0.0), // Node 0
        Vector3D(1.0, 0.0, 0.0), // Node 1
        Vector3D(0.0, 1.0, 0.0), // Node 2
        Vector3D(0.0, 0.0, 1.0)  // Node 3
    };
    
    // Define a single tetrahedron
    TetraCell cell = createSimpleTetraCell(0, {0, 1, 2, 3});
    CellVectorField field = {1.0, 1.0, 1.0}; // Velocity vector
    
    // Create Tetrahedron object
    Tetrahedron tetra(cell, nodes, field);
    
    // Get vertices
    const std::array<Vector3D, 4>& vertices = tetra.getVertices();
    EXPECT_EQ(vertices[0], Vector3D(0.0, 0.0, 0.0));
    EXPECT_EQ(vertices[1], Vector3D(1.0, 0.0, 0.0));
    EXPECT_EQ(vertices[2], Vector3D(0.0, 1.0, 0.0));
    EXPECT_EQ(vertices[3], Vector3D(0.0, 0.0, 1.0));
    
    // Get center of mass
    const Vector3D& center_of_mass = tetra.getCenterOfMass();
    EXPECT_EQ(center_of_mass, Vector3D(1.0, 1.0, 1.0));

    // Get velocity
    const Vector3D& velocity = tetra.getVelocity();
    EXPECT_EQ(velocity, Vector3D(1.0, 1.0, 1.0));
}