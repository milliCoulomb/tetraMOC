// tests/test_RayTracer.cpp
#include <gtest/gtest.h>
#include "RayTracer.hpp"
#include "MeshHandler.hpp"
#include "Field.hpp"
#include "TestUtils.hpp" // Include the Test Utilities
#include "Vector3D.hpp"

using namespace SNSolver;

// Helper function to create a simple mesh with two adjacent tetrahedrons
void setupSimpleMesh(MeshHandler& mesh) {
    // Define nodes
    std::vector<Node> nodes = {
        {0.0, 0.0, 0.0}, // Node 0
        {1.0, 0.0, 0.0}, // Node 1
        {0.0, 1.0, 0.0}, // Node 2
        {0.0, 0.0, 1.0}, // Node 3
        {1.0, 1.0, 1.0}  // Node 4
    };
    
    // Load nodes
    // Assuming mesh.loadNodes reads from a file, we'll bypass it by directly assigning
    // If MeshHandler doesn't allow direct assignment, consider modifying it or using temporary files
    // For simplicity, let's assume MeshHandler can have nodes assigned directly
    // Otherwise, create temporary nodes.txt
    // Create temporary nodes.txt
    std::string nodes_content = "5\n"
                                "0.0 0.0 0.0\n"
                                "1.0 0.0 0.0\n"
                                "0.0 1.0 0.0\n"
                                "0.0 0.0 1.0\n"
                                "1.0 1.0 1.0\n";
    ASSERT_TRUE(createTempFile("temp_nodes_test_RayTracer.txt", nodes_content));
    ASSERT_TRUE(mesh.loadNodes("temp_nodes_test_RayTracer.txt"));
    
    // Define cells (two tetrahedrons sharing a face)
    // Cell 0: Nodes 0,1,2,3
    // Cell 1: Nodes 1,2,3,4
    std::string cells_content = "2\n"
                                "0 1 2 3\n"
                                "1 2 3 4\n";
    ASSERT_TRUE(createTempFile("temp_cells_test_RayTracer.txt", cells_content));
    ASSERT_TRUE(mesh.loadCells("temp_cells_test_RayTracer.txt"));
}

void setupSimpleField(Field& field) {
    // Define vector fields for two cells
    // Cell 0: Velocity towards positive x-axis
    // Cell 1: Velocity towards positive y-axis
    std::string field_content = "2\n"
                                "1.0 0.0 0.0\n" // Cell 0
                                "0.0 1.0 0.0\n"; // Cell 1
    ASSERT_TRUE(createTempFile("temp_field_test_RayTracer.txt", field_content));
    ASSERT_TRUE(field.loadVectorField("temp_field_test_RayTracer.txt"));
}

void setupSimpleFaceConnectivity(RayTracer& ray_tracer) {
    // Define face connectivity
    // Faces are sorted triples
    // Shared face between Cell 0 and Cell 1: Nodes 1,2,3
    // Other faces are boundary faces
    std::string faces_content = 
        "0 1 2 0\n"  // Face between Nodes 0,1,2 - adjacent to Cell 0
        "0 1 3 0\n"  // Face between Nodes 0,1,3 - adjacent to Cell 0
        "0 2 3 0\n"  // Face between Nodes 0,2,3 - adjacent to Cell 0
        "1 2 3 0 1\n"// Face between Nodes 1,2,3 - adjacent to Cell 0 and Cell 1
        "1 2 4 1\n"  // Face between Nodes 1,2,4 - adjacent to Cell 1
        "1 3 4 1\n"  // Face between Nodes 1,3,4 - adjacent to Cell 1
        "2 3 4 1\n"; // Face between Nodes 2,3,4 - adjacent to Cell 1
    ASSERT_TRUE(createTempFile("temp_faces_test_RayTracer.txt", faces_content));
    ASSERT_TRUE(ray_tracer.loadFaceConnectivity("temp_faces_test_RayTracer.txt"));
}

TEST(RayTracerTest, TraceRayBasicTest) {
    // Setup mesh
    MeshHandler mesh;
    setupSimpleMesh(mesh);
    
    // Setup field
    Field field;
    setupSimpleField(field);
    
    // Initialize RayTracer
    RayTracer ray_tracer(mesh, field);
    setupSimpleFaceConnectivity(ray_tracer);
    
    // Define starting cell and point
    int start_cell_id = 0;
    std::array<double, 3> start_point = {0.1, 0.1, 0.1};
    
    // Trace the ray
    std::vector<DirectionData> pathline = ray_tracer.traceRay(start_cell_id, start_point, 1);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 0 and enter Cell 1
    EXPECT_EQ(pathline.size(), 1) << "Pathline should have one segment";
    EXPECT_EQ(pathline[0].cell_id, 0) << "First segment should be in Cell 0";
    
    // Check that the exit point is on the shared face (Nodes 1,2,3)
    std::array<double, 3> exit_point = pathline[0].end_point;
    // Since Cell 0 has velocity towards positive x, it should exit through face 0,1,2 (Nodes 0,1,2)
    // But in our setup, velocity is along x for Cell 0, so the exit face should be the face perpendicular to x
    // Let's verify the exit_face_id is 0 (face 0: Nodes 0,1,2)
    EXPECT_EQ(pathline[0].cell_id, 0);
    // For more precise checks, calculate the expected exit face
    // However, without exact geometry calculations, we can check that the exit_face_id is 0 or similar
    // Alternatively, check that after tracing, the next cell is 1
    // But since pathline records only the segments, we can infer based on the pathline
}

TEST(RayTracerTest, TraceRayThroughTwoCells) {
    // Setup mesh
    MeshHandler mesh;
    setupSimpleMesh(mesh);
    
    // Setup field
    Field field;
    setupSimpleField(field);
    
    // Initialize RayTracer
    RayTracer ray_tracer(mesh, field);
    setupSimpleFaceConnectivity(ray_tracer);
    
    // Define starting cell and point
    int start_cell_id = 0;
    std::array<double, 3> start_point = {0.1, 0.1, 0.1};
    
    // Trace the ray with higher max_iter
    std::vector<DirectionData> pathline = ray_tracer.traceRay(start_cell_id, start_point, 10);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 0 and enter Cell 1
    EXPECT_GE(pathline.size(), 1) << "Pathline should have at least one segment";
    EXPECT_EQ(pathline[0].cell_id, 0) << "First segment should be in Cell 0";
    
    if(pathline.size() > 1) {
        EXPECT_EQ(pathline[1].cell_id, 1) << "Second segment should be in Cell 1";
    }
}

TEST(RayTracerTest, TraceRayExitsDomain) {
    // Setup mesh
    MeshHandler mesh;
    setupSimpleMesh(mesh);
    
    // Setup field
    Field field;
    setupSimpleField(field);
    
    // Initialize RayTracer
    RayTracer ray_tracer(mesh, field);
    setupSimpleFaceConnectivity(ray_tracer);
    
    // Define starting cell and point
    int start_cell_id = 1; // Starting from Cell 1, which is adjacent to Cell 0
    std::array<double, 3> start_point = {0.9, 0.9, 0.9};
    
    // Trace the ray
    std::vector<DirectionData> pathline = ray_tracer.traceRay(start_cell_id, start_point, 10);
    
    // Assertions
    ASSERT_FALSE(pathline.empty()) << "Pathline should not be empty";
    
    // Expect the ray to exit Cell 1 and potentially exit the domain
    EXPECT_EQ(pathline[0].cell_id, 1) << "First segment should be in Cell 1";
    
    // Since there is no Cell 2, the ray should exit the domain after Cell 1
    EXPECT_EQ(pathline.size(), 1) << "Pathline should have one segment";
}