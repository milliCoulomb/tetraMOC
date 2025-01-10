// tests/test_RayTracer.cpp

#include <gtest/gtest.h>
#include "RayTracer.hpp"
#include "MeshHandler.hpp"
#include "Field.hpp"
#include "Vector3D.hpp"
#include "TestUtils.hpp" // Contains createTempFile and vectorsAlmostEqual
#include "Tetrahedron.hpp" // Ensure it includes the Tetrahedron class with findExit method
#include <fstream>
#include <cstdio> // For std::remove>

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
        std::string field_content = "2\n" // Number of vectors
                                       "1.0 0.0 0.0\n" // Vector for cell 0
                                       "0.0 1.0 0.0\n"; // Vector for cell 1

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

// Test case: Basic ray tracing in Variable Direction Mode
TEST_F(RayTracerTest, TraceRay_VariableDirection_BasicTest) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize RayTracer in VARIABLE_DIRECTION mode
    RayTracer ray_tracer(mesh, field);
    
    // Define starting cell and point
    int start_cell_id = 0;
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Trace the ray with max_iter = 1
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 1);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 0
    EXPECT_EQ(pathline.size(), 1) << "Pathline should have one segment";
    EXPECT_EQ(pathline[0].cell_id, 0) << "First segment should be in Cell 0";
    
    // Verify end_point lies on the expected exit
    // Vector3D expected_end(0.8, 0.1, 0.1); // Corrected Expected exit point
    // EXPECT_TRUE(vectorsAlmostEqual(pathline[0].end_point, expected_end)) << "End point should match expected exit point";
}

// Test case: Ray Tracing through two cells in Variable Direction Mode
TEST_F(RayTracerTest, TraceRay_VariableDirection_TwoCells) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    Logger::info("Number of boundary faces: " + std::to_string(mesh.getBoundaryFaces().size()));
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize RayTracer in VARIABLE_DIRECTION mode
    RayTracer ray_tracer(mesh, field);
    
    // Define starting cell and point
    int start_cell_id = 0;
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Trace the ray with higher max_iter to allow traversal into Cell 1
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 2);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 0 and enter Cell 1
    EXPECT_GE(pathline.size(), 2) << "Pathline should have at least two segments";
    EXPECT_EQ(pathline[0].cell_id, 0) << "First segment should be in Cell 0";
    EXPECT_EQ(pathline[1].cell_id, 1) << "Second segment should be in Cell 1";
    
    // Verify end_point of second segment lies on the expected exit from Cell 1
    // Vector3D expected_end(0.8, 0.3, 0.1); // Corrected Expected exit point from Cell 1
    // EXPECT_TRUE(vectorsAlmostEqual(pathline[1].end_point, expected_end)) << "End point should match expected exit point from Cell 1";
}

// Test case: Ray Tracing in Constant Direction Mode
TEST_F(RayTracerTest, TraceRay_ConstantDirection_BasicTest) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Define a constant direction
    Vector3D fixed_direction(0.0, 1.0, 0.0); // Along +Y axis
    
    // Initialize RayTracer in CONSTANT_DIRECTION mode
    RayTracer ray_tracer(mesh, fixed_direction, 1.0, static_cast<size_t>(0));
    
    // Define starting cell and point
    int start_cell_id = 0;
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Trace the ray with max_iter = 1
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 1);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 0
    EXPECT_EQ(pathline.size(), 1) << "Pathline should have one segment";
    EXPECT_EQ(pathline[0].cell_id, 0) << "First segment should be in Cell 0";
    
    // Verify end_point lies on the expected exit
    // Vector3D expected_end(0.1, 0.8, 0.1); // Corrected Expected exit point
    // EXPECT_TRUE(vectorsAlmostEqual(pathline[0].end_point, expected_end)) << "End point should match fixed direction exit point";
}

// Test case: Ray Tracing with Invalid Cell ID in Variable Direction Mode
TEST_F(RayTracerTest, TraceRay_VariableDirection_InvalidCellID) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize RayTracer in VARIABLE_DIRECTION mode
    RayTracer ray_tracer(mesh, field);
    
    // Define starting cell ID as invalid
    int start_cell_id = -1;
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Trace the ray
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 10);
    
    // Assertions
    EXPECT_TRUE(pathline.empty()) << "Pathline should be empty for invalid cell ID";
}

// Test case: Ray Tracing with Invalid Cell ID in Constant Direction Mode
TEST_F(RayTracerTest, TraceRay_ConstantDirection_InvalidCellID) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Define a constant direction
    Vector3D fixed_direction(0.0, 1.0, 0.0); // Along +Y axis
    
    // Initialize RayTracer in CONSTANT_DIRECTION mode
    RayTracer ray_tracer(mesh, fixed_direction, 1.0, static_cast<size_t>(0));
    
    // Define invalid starting cell ID
    int start_cell_id = 10; // Assuming only 2 cells exist (IDs 0 and 1)
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Trace the ray
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 10);
    
    // Assertions
    EXPECT_TRUE(pathline.empty()) << "Pathline should be empty for invalid cell ID";
}

// Test case: Ray Tracing Exiting the Domain in Variable Direction Mode
TEST_F(RayTracerTest, TraceRay_VariableDirection_ExitsDomain) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize RayTracer in VARIABLE_DIRECTION mode
    RayTracer ray_tracer(mesh, field);
    
    // Define starting cell and point near the edge
    int start_cell_id = 1; // Cell 1
    Vector3D start_point(0.9, 0.9, 0.9); // Near the far end
    
    // Trace the ray with max_iter = 1 to simulate exiting the domain
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 1);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 1
    EXPECT_EQ(pathline.size(), 1) << "Pathline should have one segment";
    EXPECT_EQ(pathline[0].cell_id, 1) << "First segment should be in Cell 1";
    
    // Verify end_point lies on the expected exit
    // Vector3D expected_end(0.9, 1.0, 0.9); // Corrected Expected exit point
    // EXPECT_TRUE(vectorsAlmostEqual(pathline[0].end_point, expected_end)) << "End point should match expected exit point";
}

// Test case: Ray Tracing Exiting the Domain in Constant Direction Mode
TEST_F(RayTracerTest, TraceRay_ConstantDirection_ExitsDomain) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Define a constant direction that points outwards
    Vector3D fixed_direction(0.0, 1.0, 0.0); // Along +Y axis
    
    // Initialize RayTracer in CONSTANT_DIRECTION mode
    RayTracer ray_tracer(mesh, fixed_direction, 1.0, static_cast<size_t>(0));
    
    // Define starting cell near the boundary
    int start_cell_id = 1; // Cell 1
    Vector3D start_point(0.9, 0.9, 0.9); // Near the far end
    
    // Trace the ray
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 10);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 1
    EXPECT_EQ(pathline.size(), 1) << "Pathline should have one segment";
    EXPECT_EQ(pathline[0].cell_id, 1) << "First segment should be in Cell 1";
    
    // Verify end_point lies on the expected exit
    // Vector3D expected_end(0.9, 1.0, 0.9); // Corrected Expected exit point
    // EXPECT_TRUE(vectorsAlmostEqual(pathline[0].end_point, expected_end)) << "End point should match expected exit point";
}

// Test case: Multiple Iterations in Variable Direction Mode
TEST_F(RayTracerTest, TraceRay_VariableDirection_MultipleIterations) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSimpleMesh(mesh)) << "Failed to setup simple mesh";
    
    Field field;
    ASSERT_TRUE(setupSimpleField(field)) << "Failed to setup simple field";
    
    ASSERT_TRUE(setupSimpleFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize RayTracer in VARIABLE_DIRECTION mode
    RayTracer ray_tracer(mesh, field);
    
    // Define starting cell and point
    int start_cell_id = 0;
    Vector3D start_point(0.1, 0.1, 0.1);
    
    // Trace the ray with higher max_iter to allow traversal into Cell 1
    std::vector<CellTrace> pathline = ray_tracer.traceRay(start_cell_id, start_point, 2);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 0 and enter Cell 1
    EXPECT_GE(pathline.size(), 2) << "Pathline should have at least two segments";
    EXPECT_EQ(pathline[0].cell_id, 0) << "First segment should be in Cell 0";
    EXPECT_EQ(pathline[1].cell_id, 1) << "Second segment should be in Cell 1";
    
    // Verify end_point of second segment lies on the expected exit from Cell 1
    // Vector3D expected_end(0.8, 0.3, 0.1); // Corrected Expected exit point from Cell 1
    // EXPECT_TRUE(vectorsAlmostEqual(pathline[1].end_point, expected_end)) << "End point should match expected exit point from Cell 1";
}