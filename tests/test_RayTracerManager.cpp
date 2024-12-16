// tests/test_RayTracerManager.cpp
#include <gtest/gtest.h>
#include "MeshHandler.hpp"
#include "Field.hpp"
#include "RayTracer.hpp"
#include "AngularQuadrature.hpp"
#include "RayTracerManager.hpp"
#include "GeometryUtils.hpp"
#include "TestUtils.hpp" // Include TestUtils for temporary file utilities
#include <algorithm> // For std::sort

// Namespace shortcut
namespace SNS = SNSolver;

// Helper function to create a simple mesh with two adjacent tetrahedrons
void setupSimpleMesh(SNS::MeshHandler& mesh) {
    // Define nodes content
    std::string nodes_content = 
        "0.0 0.0 0.0\n" // Node 0
        "1.0 0.0 0.0\n" // Node 1
        "0.0 1.0 0.0\n" // Node 2
        "0.0 0.0 1.0\n" // Node 3
        "1.0 1.0 1.0\n"; // Node 4
    
    // Define cells content
    std::string cells_content = 
        "0 0 1 2 3\n" // Cell 0: Nodes 0,1,2,3
        "1 1 2 3 4\n"; // Cell 1: Nodes 1,2,3,4
    
    // Define face connectivity content
    // Format: face_id n0 n1 n2 cell_id(s)
    std::string faces_content = 
        "0 0 1 2 0\n" // Face (0,1,2) adjacent to Cell 0
        "1 0 1 3 0\n" // Face (0,1,3) adjacent to Cell 0
        "2 0 2 3 0\n" // Face (0,2,3) adjacent to Cell 0
        "3 1 2 3 0 1\n" // Shared Face (1,2,3) adjacent to Cells 0 and 1
        "4 1 2 4 1\n" // Face (1,2,4) adjacent to Cell 1
        "5 1 3 4 1\n" // Face (1,3,4) adjacent to Cell 1
        "6 2 3 4 1\n"; // Face (2,3,4) adjacent to Cell 1
    
    // Generate unique temporary filenames
    std::string nodes_filename = SNS::generateUniqueFilename("test_mesh_nodes_", ".txt");
    std::string cells_filename = SNS::generateUniqueFilename("test_mesh_cells_", ".txt");
    std::string faces_filename = SNS::generateUniqueFilename("test_mesh_faces_", ".txt");
    
    // Create temporary files
    ASSERT_TRUE(SNS::createTempFile(nodes_filename, nodes_content)) << "Failed to create nodes temporary file.";
    ASSERT_TRUE(SNS::createTempFile(cells_filename, cells_content)) << "Failed to create cells temporary file.";
    ASSERT_TRUE(SNS::createTempFile(faces_filename, faces_content)) << "Failed to create faces temporary file.";
    
    // Load mesh data
    ASSERT_TRUE(mesh.loadNodes(nodes_filename)) << "Failed to load nodes from temporary file.";
    ASSERT_TRUE(mesh.loadCells(cells_filename)) << "Failed to load cells from temporary file.";
    ASSERT_TRUE(mesh.loadFaceConnectivity(faces_filename)) << "Failed to load face connectivity from temporary file.";
    
    // Optionally, remove temporary files if not needed
    std::remove(nodes_filename.c_str());
    std::remove(cells_filename.c_str());
    std::remove(faces_filename.c_str());
}

TEST(RayTracerManagerTest, BoundaryFaceIdentification) {
    // Setup mesh
    SNS::MeshHandler mesh;
    setupSimpleMesh(mesh);
    
    // Setup field (dummy, not used in this test)
    SNS::Field field;
    
    // Initialize RayTracer
    SNS::RayTracer tracer(mesh, field);
    
    // Retrieve boundary faces
    std::vector<SNS::Face> boundary_faces = tracer.getBoundaryFaces();
    
    // Expected boundary faces: All except the shared face (1,2,3)
    // Total boundary faces: 6 (since each tetrahedron has 4 faces, one shared)
    ASSERT_EQ(boundary_faces.size(), 6) << "Expected 6 boundary faces.";
    
    // Verify that none of the boundary faces is the shared face (1,2,3)
    for(const auto& face : boundary_faces) {
        // Sort the face node indices for consistent comparison
        std::array<int,3> sorted_face = {face.n0, face.n1, face.n2};
        std::sort(sorted_face.begin(), sorted_face.end());
        
        // Shared face sorted: {1,2,3}
        EXPECT_NE(sorted_face[0], 1) << "Shared face (1,2,3) incorrectly identified as boundary face.";
        EXPECT_NE(sorted_face[1], 2) << "Shared face (1,2,3) incorrectly identified as boundary face.";
        EXPECT_NE(sorted_face[2], 3) << "Shared face (1,2,3) incorrectly identified as boundary face.";
    }
}
TEST(RayTracerManagerTest, PointSamplingOnTriangle) {
    // Define a simple triangle in 3D space
    std::array<SNS::Node, 3> triangle = {
        SNS::Node{0.0, 0.0, 0.0}, // Node A
        SNS::Node{1.0, 0.0, 0.0}, // Node B
        SNS::Node{0.0, 1.0, 0.0}  // Node C
    };
    
    // Number of samples
    int num_samples = 1000;
    
    for(int i = 0; i < num_samples; ++i) {
        std::array<double, 3> point = SNS::samplePointOnTriangle(triangle);
        
        // Barycentric coordinates check
        double u = point[0];
        double v = point[1];
        double w = 1.0 - u - v;
        
        // Point should satisfy: u >= 0, v >= 0, w >= 0
        EXPECT_GE(u, 0.0) << "Sampled u coordinate is negative.";
        EXPECT_GE(v, 0.0) << "Sampled v coordinate is negative.";
        EXPECT_GE(w, 0.0) << "Sampled w coordinate is negative.";
        
        // Additionally, u + v + w should be approximately 1
        EXPECT_NEAR(u + v + w, 1.0, 1e-6) << "Sum of barycentric coordinates deviates from 1.";
    }
}
TEST(RayTracerManagerTest, DirectionValidation) {
    // Define a face normal pointing in +z direction
    std::array<double, 3> face_normal = {0.0, 0.0, 1.0};
    
    // Initialize MeshHandler with simple mesh (not used in this test)
    SNS::MeshHandler mesh;
    setupSimpleMesh(mesh);
    
    // Initialize Field with dummy data
    SNS::Field field;
    
    // Initialize RayTracer
    SNS::RayTracer tracer(mesh, field);
    
    // Initialize AngularQuadrature with dummy orders
    SNS::AngularQuadrature angular_quadrature(1, 1); // Dummy orders
    // Assuming angular_quadrature.getDirections() returns at least one direction
    
    // Initialize RayTracerManager
    SNS::RayTracerManager manager(mesh, field, tracer, angular_quadrature);
    
    // Define directions
    // 0 degrees (parallel, should be rejected)
    SNS::Direction dir_parallel_positive(1.0, 0.0, 1.0); // mu = cos(0°) = 1, phi = 0°
    // 90 degrees (perpendicular, should be accepted)
    SNS::Direction dir_perpendicular(0.0, 0.0, 0.0); // mu = cos(90°) = 0, phi = 0°
    // 45 degrees, should be accepted
    double theta_45 = M_PI / 4; // 45 degrees in radians
    SNS::Direction dir_45_degrees(std::cos(theta_45), M_PI / 4, 1.0); // mu = cos(45°), phi = 45°
    // 179 degrees (parallel in opposite direction, should be rejected)
    double theta_179 = 179.0 * M_PI / 180.0; // 179 degrees in radians
    SNS::Direction dir_parallel_negative(std::cos(theta_179), 0.0, 1.0); // mu = cos(179°), phi = 0°
    
    // Threshold angle
    double threshold = 1.0; // degrees
    
    // Test parallel positive
    EXPECT_FALSE(manager.isValidDirection(face_normal, dir_parallel_positive, threshold)) 
        << "Parallel positive direction should be invalid.";
    
    // Test perpendicular
    EXPECT_TRUE(manager.isValidDirection(face_normal, dir_perpendicular, threshold)) 
        << "Perpendicular direction should be valid.";
    
    // Test 45 degrees
    EXPECT_TRUE(manager.isValidDirection(face_normal, dir_45_degrees, threshold)) 
        << "45-degree direction should be valid.";
    
    // Test parallel negative
    EXPECT_FALSE(manager.isValidDirection(face_normal, dir_parallel_negative, threshold)) 
        << "Parallel negative direction should be invalid.";
}
TEST(RayTracerManagerTest, RayTracingSimpleMesh) {
    // Setup mesh
    SNS::MeshHandler mesh;
    setupSimpleMesh(mesh);
    
    // Setup field (dummy velocities)
    SNS::Field field;
    
    // Initialize RayTracer
    SNS::RayTracer tracer(mesh, field);
    
    // Define a direction that intersects the shared face (1,2,3)
    // For simplicity, use a direction from Cell 0 to Cell 1
    // This corresponds to mu = cos(theta), phi = 45°, theta = 45°
    double theta = M_PI / 4; // 45 degrees in radians
    double mu = std::cos(theta); // cos(45°) ≈ 0.7071
    double phi = M_PI / 4; // 45 degrees in radians
    SNS::Direction dir(mu, phi, 1.0); // weight is arbitrary here
    
    // Initialize AngularQuadrature with a single direction for testing
    SNS::AngularQuadrature angular_quadrature(1, 1); // Dummy orders
    // Assuming angular_quadrature.getDirections() returns the desired direction
    // If not, you may need to mock or adjust the AngularQuadrature class
    
    // Initialize RayTracerManager
    SNS::RayTracerManager manager(mesh, field, tracer, angular_quadrature);
    
    // Generate tracking data with 1 ray per boundary face
    manager.generateTrackingData(1);
    
    // Retrieve tracking data
    const auto& tracking_data = manager.getTrackingData();
    
    // Since there are 6 boundary faces and 1 ray per face with half_directions = angular_quadrature.getDirections().size() / 2 = 0 (assuming 2 directions)
    // Adjust expected_rays based on angular_quadrature's getDirections()
    // For this test, focus on verifying at least one correct traversal
    
    bool found_correct_traversal = false;
    
    for(const auto& data : tracking_data) {
        // Check if any ray traverses from Cell 0 to Cell 1
        if(data.cell_traces.size() >= 2 &&
           data.cell_traces[0].cell_id == 0 &&
           data.cell_traces[1].cell_id == 1) {
            found_correct_traversal = true;
            break;
        }
    }
    
    EXPECT_TRUE(found_correct_traversal) << "Expected to find at least one ray traversing from Cell 0 to Cell 1.";
}
TEST(RayTracerManagerTest, TrackingDataStorage) {
    // Setup mesh
    SNS::MeshHandler mesh;
    setupSimpleMesh(mesh);
    
    // Setup field (dummy velocities)
    SNS::Field field;
    
    // Initialize RayTracer
    SNS::RayTracer tracer(mesh, field);
    
    // Initialize AngularQuadrature with multiple directions
    // For example, orders 2,2 might generate four directions
    // Assuming AngularQuadrature::getDirections() returns four directions
    SNS::AngularQuadrature angular_quadrature(2, 2);
    
    // Initialize RayTracerManager
    SNS::RayTracerManager manager(mesh, field, tracer, angular_quadrature);
    
    // Generate tracking data with 2 rays per boundary face
    manager.generateTrackingData(2);
    
    // Retrieve tracking data
    const auto& tracking_data = manager.getTrackingData();
    
    // Number of boundary faces in simple mesh: 6
    // Rays per face: 2
    // Half_directions: angular_quadrature.getDirections().size() / 2 = 2
    // Total expected rays: 6 * 2 * 2 = 24
    int half_directions = angular_quadrature.getDirections().size() / 2;
    int expected_rays = 6 * 2 * half_directions;
    
    ASSERT_EQ(tracking_data.size(), expected_rays) 
        << "Expected " << expected_rays << " tracking data entries, but got " << tracking_data.size() << ".";
    
    // Verify uniqueness of ray_ids
    std::vector<int> ray_ids;
    for(const auto& data : tracking_data) {
        ray_ids.push_back(data.ray_id);
    }
    std::sort(ray_ids.begin(), ray_ids.end());
    auto last = std::unique(ray_ids.begin(), ray_ids.end());
    ASSERT_EQ(last - ray_ids.begin(), tracking_data.size()) 
        << "ray_id should be unique for each TrackingData.";
    
    // Verify that cell_traces are not empty and contain valid cell IDs
    for(const auto& data : tracking_data) {
        ASSERT_FALSE(data.cell_traces.empty()) 
            << "cell_traces should not be empty for ray_id " << data.ray_id;
        for(const auto& trace : data.cell_traces) {
            EXPECT_GE(trace.cell_id, 0) 
                << "Invalid cell_id " << trace.cell_id << " in ray_id " << data.ray_id;
            EXPECT_LT(trace.cell_id, static_cast<int>(mesh.getCells().size())) 
                << "cell_id " << trace.cell_id << " out of range in ray_id " << data.ray_id;
        }
    }
}
