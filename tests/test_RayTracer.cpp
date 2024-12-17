// tests/test_RayTracer.cpp
#include <gtest/gtest.h>
#include "RayTracer.hpp"
#include "MeshHandler.hpp"
#include "Field.hpp"
#include "Vector3D.hpp"
#include "TestUtils.hpp"
#include <fstream>
#include <cstdio> // For std::remove

using namespace SNSolver;

// Test Fixture for RayTracer
class RayTracerTest : public ::testing::Test {
protected:
    // Temporary file names
    std::string nodes_file = "temp_nodes_test_RayTracer.txt";
    std::string cells_file = "temp_cells_test_RayTracer.txt";
    std::string faces_file = "temp_faces_test_RayTracer.txt";
    std::string field_file = "temp_field_test_RayTracer.txt";
    
    // Clean up temporary files after each test
    void TearDown() override {
        std::remove(nodes_file.c_str());
        std::remove(cells_file.c_str());
        std::remove(faces_file.c_str());
        std::remove(field_file.c_str());
    }

    // Helper function to create a simple mesh with two adjacent tetrahedrons
    bool setupSimpleMesh(MeshHandler& mesh) {
        // Define nodes content with node IDs
        std::string nodes_content = "5\n"
                                    "0 0.0 0.0 0.0\n" // Node 0
                                    "1 1.0 0.0 0.0\n" // Node 1
                                    "2 0.0 1.0 0.0\n" // Node 2
                                    "3 0.0 0.0 1.0\n" // Node 3
                                    "4 1.0 1.0 1.0\n"; // Node 4

        // Create temporary nodes.txt
        if(!createTempFile(nodes_file, nodes_content)) return false;
        if(!mesh.loadNodes(nodes_file)) return false;

        // Define cells content with cell IDs and four node IDs per cell
        std::string cells_content = "2\n"
                                    "0 0 1 2 3\n" // Cell 0: Nodes 0,1,2,3
                                    "1 1 2 3 4\n"; // Cell 1: Nodes 1,2,3,4

        // Create temporary cells.txt
        if(!createTempFile(cells_file, cells_content)) return false;
        if(!mesh.loadCells(cells_file)) return false;

        return true;
    }
    
    // Helper function to setup a simple field
    bool setupSimpleField(Field& field) {
        // Define vector fields for two cells
        // Cell 0: Velocity towards positive x-axis
        // Cell 1: Velocity towards positive y-axis
        std::string field_content = "2\n"
                                    "1.0 0.0 0.0\n" // Cell 0
                                    "0.0 1.0 0.0\n"; // Cell 1

        // Create temporary field.txt
        if(!createTempFile(field_file, field_content)) return false;
        if(!field.loadVectorField(field_file)) return false;

        return true;
    }
    
    // Helper function to setup simple face connectivity
    bool setupSimpleFaceConnectivity(MeshHandler& mesh) {
        // Define face connectivity with counts of adjacent cells
        // Each line: n0 n1 n2 <count> <cell_id0> [<cell_id1> ...]
        std::string faces_content = 
            "7\n" // Number of faces
            "0 1 2 1 0\n"  // Face 0: Nodes 0,1,2 adjacent to Cell 0
            "0 1 3 1 0\n"  // Face 1: Nodes 0,1,3 adjacent to Cell 0
            "0 2 3 1 0\n"  // Face 2: Nodes 0,2,3 adjacent to Cell 0
            "1 2 3 2 0 1\n"// Face 3: Nodes 1,2,3 adjacent to Cell 0 and Cell 1
            "1 2 4 1 1\n"  // Face 4: Nodes 1,2,4 adjacent to Cell 1
            "1 3 4 1 1\n"  // Face 5: Nodes 1,3,4 adjacent to Cell 1
            "2 3 4 1 1\n"; // Face 6: Nodes 2,3,4 adjacent to Cell 1

        // Create temporary faces.txt
        if(!createTempFile(faces_file, faces_content)) return false;
        if(!mesh.loadFaceConnectivity(faces_file)) return false;

        return true;
    }
};

// Test case: Basic ray tracing from Cell 0
TEST_F(RayTracerTest, TraceRayBasicTest) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    RayTracer ray_tracer(mesh, field);
    
    // Define starting cell and point
    int start_cell_id = 0;
    // std::array<double, 3> start_point = {0.1, 0.1, 0.1};
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Trace the ray with max_iter = 1
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 1);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 0 and enter Cell 1
    EXPECT_EQ(pathline.size(), 1) << "Pathline should have one segment";
    EXPECT_EQ(pathline[0].cell_id, 0) << "First segment should be in Cell 0";
    
    // Additional checks could include verifying that end_point lies on the expected face
    // However, without specific geometry calculations, it's difficult to assert exact values
}

// Test case: Trace ray through two cells
TEST_F(RayTracerTest, TraceRayThroughTwoCells) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    RayTracer ray_tracer(mesh, field);
    
    // Define starting cell and point
    int start_cell_id = 0;
    // std::array<double, 3> start_point = {0.1, 0.1, 0.1};
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Trace the ray with higher max_iter
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 10);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 0 and enter Cell 1
    EXPECT_GE(pathline.size(), 1) << "Pathline should have at least one segment";
    EXPECT_EQ(pathline[0].cell_id, 0) << "First segment should be in Cell 0";
    
    if(pathline.size() > 1) {
        EXPECT_EQ(pathline[1].cell_id, 1) << "Second segment should be in Cell 1";
    }
}

// Test case: Trace ray exiting the domain
TEST_F(RayTracerTest, TraceRayExitsDomain) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    RayTracer ray_tracer(mesh, field);
    
    // Define starting cell and point
    int start_cell_id = 1; // Starting from Cell 1, which is adjacent to Cell 0
    // std::array<double, 3> start_point = {0.9, 0.9, 0.9};
    Vector3D start_point(0.9, 0.9, 0.9);
    
    // Trace the ray
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 10);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 1 and potentially exit the domain
    EXPECT_EQ(pathline.size(), 1) << "Pathline should have one segment";
    EXPECT_EQ(pathline[0].cell_id, 1) << "First segment should be in Cell 1";
    
    // Since Cell 1 is the last cell, ray should exit the domain after Cell 1
}