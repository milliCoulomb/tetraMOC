// tests/test_tetrahedron.cpp
#include <gtest/gtest.h>
#include "Tetrahedron.hpp"
#include "Vector3D.hpp"

namespace SNS = SNSolver;

// Helper function to create a simple tetrahedron
TetraCell createSimpleTetraCell(int cell_id, const std::array<int, 4>& node_ids) {
    TetraCell cell;
    cell.cell_id = cell_id;
    cell.node_ids = node_ids;
    return cell;
}

TEST(TetrahedronTest, FindExitBasicTest) {
    // Define nodes
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0}, // Node 0
        {1.0, 0.0, 0.0}, // Node 1
        {0.0, 1.0, 0.0}, // Node 2
        {0.0, 0.0, 1.0}  // Node 3
    };
    
    // Define a single tetrahedron
    TetraCell cell = createSimpleTetraCell(0, {0, 1, 2, 3});
    CellField field = {1.0, 1.0, 1.0}; // Velocity vector
    
    // Create Tetrahedron object
    Tetrahedron tetra(cell, nodes, field);
    
    // Define a starting point inside the tetrahedron
    std::array<double, 3> start_point = {0.1, 0.1, 0.1};
    
    // Define velocity vector
    SNS::Vector3D velocity(1.0, 1.0, 1.0); // Diagonal direction
    
    double t_exit;
    std::array<double, 3> x_exit;
    int exit_face_id;
    
    // Call findExit
    bool has_exit = tetra.findExit(start_point, velocity, t_exit, x_exit, exit_face_id);
    
    // Assertions
    EXPECT_TRUE(has_exit);
    EXPECT_GT(t_exit, 0.0);
    EXPECT_NE(exit_face_id, -1);
    
    // Check that exit point is within the expected bounds
    EXPECT_GE(x_exit[0], 0.0);
    EXPECT_LE(x_exit[0], 1.0);
    EXPECT_GE(x_exit[1], 0.0);
    EXPECT_LE(x_exit[1], 1.0);
    EXPECT_GE(x_exit[2], 0.0);
    EXPECT_LE(x_exit[2], 1.0);
}

TEST(TetrahedronTest, NoExitTest) {
    // Define nodes
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0}, // Node 0
        {1.0, 0.0, 0.0}, // Node 1
        {0.0, 1.0, 0.0}, // Node 2
        {0.0, 0.0, 1.0}  // Node 3
    };
    
    // Define a single tetrahedron
    TetraCell cell = createSimpleTetraCell(0, {0, 1, 2, 3});
    CellField field = {-1.0, -1.0, -1.0}; // Velocity vector pointing inward
    
    // Create Tetrahedron object
    Tetrahedron tetra(cell, nodes, field);
    
    // Define a starting point inside the tetrahedron
    std::array<double, 3> start_point = {0.1, 0.1, 0.1};
    
    // Define velocity vector
    SNS::Vector3D velocity(-1.0, -1.0, -1.0); // Inward direction
    
    double t_exit;
    std::array<double, 3> x_exit;
    int exit_face_id;
    
    // Call findExit
    bool has_exit = tetra.findExit(start_point, velocity, t_exit, x_exit, exit_face_id);
    
    // Assertions
    EXPECT_FALSE(has_exit);
}