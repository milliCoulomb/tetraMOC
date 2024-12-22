// tests/test_FluxSolver.cpp

#include <gtest/gtest.h>
#include "FluxSolver.hpp"
#include "MeshHandler.hpp"
#include "Field.hpp"
#include "AngularQuadrature.hpp"
#include "TrackingData.hpp"
#include "Vector3D.hpp"
#include "TestUtils.hpp" // Contains createTempFile and directionVector
#include "Logger.hpp"

#include <fstream>
#include <cstdio> // For std::remove
#include <cmath>
#include <vector>
#include <string>

// Test Fixture for FluxSolver
class FluxSolverTest : public ::testing::Test {
protected:
    // Temporary file names
    std::string nodes_file = "temp_nodes_test_FluxSolver.txt";
    std::string cells_file = "temp_cells_test_FluxSolver.txt";
    std::string faces_file = "temp_faces_test_FluxSolver.txt";
    std::string field_file = "temp_field_test_FluxSolver.txt";

    // Clean up temporary files after each test
    void TearDown() override {
        std::remove(nodes_file.c_str());
        std::remove(cells_file.c_str());
        std::remove(faces_file.c_str());
        std::remove(field_file.c_str());
    }

    // Helper function to create a simple mesh with a single tetrahedron
    bool setupSingleCellMesh(MeshHandler& mesh) {
        // Define nodes content with node IDs
        std::string nodes_content = "4\n"
                                    "0 0.0 0.0 0.0\n" // Node 0
                                    "1 1.0 0.0 0.0\n" // Node 1
                                    "2 0.0 1.0 0.0\n" // Node 2
                                    "3 0.0 0.0 1.0\n"; // Node 3

        // Create temporary nodes.txt
        if(!createTempFile(nodes_file, nodes_content)) return false;
        if(!mesh.loadNodes(nodes_file)) return false;

        // Define cells content with cell IDs and four node IDs per cell
        std::string cells_content = "1\n"
                                    "0 0 1 2 3\n"; // Cell 0: Nodes 0,1,2,3

        // Create temporary cells.txt
        if(!createTempFile(cells_file, cells_content)) return false;
        if(!mesh.loadCells(cells_file)) return false;

        return true;
    }

    // helper function to create two cells mesh
    bool setupTwoCellMesh(MeshHandler& mesh) {
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

    // Helper function to setup a simple field with a single source term
    bool setupSingleCellField(Field& field) {
        // Define scalar fields content
        // Assuming the first number is the number of scalar fields
        // Followed by lines: cell_id Q_k
        std::string field_content = "1\n" // Number of scalar fields
                                    "1.0\n"; // Cell 0: Q_k = 1.0

        // Create temporary field.txt
        if(!createTempFile(field_file, field_content)) return false;
        if(!field.loadScalarField(field_file)) return false;

        return true;
    }

    // Helper function to setup a simple field with two source terms
    bool setupTwoCellField(Field& field) {
        // Define scalar fields content
        // Assuming the first number is the number of scalar fields
        // Followed by lines: cell_id Q_k
        std::string field_content = "2\n" // Number of scalar fields
                                    "1.0\n" // Cell 0: Q_k = 1.0
                                    "2.0\n"; // Cell 1: Q_k = 2.0

        // Create temporary field.txt
        if(!createTempFile(field_file, field_content)) return false;
        if(!field.loadScalarField(field_file)) return false;

        return true;
    }

    // Helper function to setup simple face connectivity
    bool setupSingleCellFaceConnectivity(MeshHandler& mesh) {
        // Define face connectivity with counts of adjacent cells
        // Each line: n0 n1 n2 <count> <cell_id0> [<cell_id1> ...]
        std::string faces_content = 
            "4\n" // Number of faces
            "0 1 2 1 0\n"  // Face 0: Nodes 0,1,2 adjacent to Cell 0
            "0 1 3 1 0\n"  // Face 1: Nodes 0,1,3 adjacent to Cell 0
            "0 2 3 1 0\n"  // Face 2: Nodes 0,2,3 adjacent to Cell 0
            "1 2 3 1 0\n"; // Face 3: Nodes 1,2,3 adjacent to Cell 0

        // Create temporary faces.txt
        if(!createTempFile(faces_file, faces_content)) return false;
        if(!mesh.loadFaceConnectivity(faces_file)) return false;

        return true;
    }

    // Helper function to setup simple face connectivity for two cells
    bool setupTwoCellFaceConnectivity(MeshHandler& mesh) {
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

    // Helper function to create a ray traversing the single cell in a specific direction
    TrackingData createSingleRay(int ray_id, const Vector3D& direction, const int cell_id = 0, const double L_k = 1.0, Vector3D start_point = Vector3D(0.0, 0.0, 0.0), Vector3D end_point = Vector3D(1.0, 0.0, 0.0)) {
        TrackingData ray;
        ray.ray_id = ray_id;
        ray.direction = direction.normalized(); // Ensure direction is unit vector

        // Single CellTrace: traversing Cell 0
        CellTrace trace;
        trace.cell_id = cell_id;
        trace.time_spent = L_k; // Time spent in the cell
        trace.start_point = start_point; // Entry point
        trace.end_point = end_point; // Exit point

        ray.cell_traces.push_back(trace);

        return ray;
    }
};

// Test case: Single Cell, Single Direction, Single Ray
TEST_F(FluxSolverTest, SingleCellSingleDirectionSingleRay) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with one direction
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create TrackingData with a single ray in the x-direction
    Vector3D dir_vector = directionVector(dir.mu, dir.phi);
    TrackingData ray = createSingleRay(0, dir_vector, 0, 1.0, Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    std::vector<TrackingData> tracking_data = { ray };
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();
    
    // Expected calculations:
    // psi_in = 0
    // psi_out = (1.0 / 1.0 - 0) * (1 - e^{-1}) ≈ 0.6321205588
    // line_avg_flux = (1.0 / 1.0) - (0.6321205588 - 0) / (1.0 * 1.0) ≈ 0.3678794412
    // flux = line_avg_flux * L_k = 0.3678794412 * 1.0 = 0.3678794412
    // weight = L_k = 1.0
    // Normalized flux = flux / weight = 0.3678794412 / 1.0 = 0.3678794412
    
    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction";
    
    EXPECT_NEAR(flux_data[0][0].flux, 0.3678794412, 1e-6) << "Flux should match the expected line-averaged flux";
    EXPECT_NEAR(flux_data[0][0].weight, 1.0, 1e-6) << "Weight should be equal to L_k (1.0)";
}

// Test case: Single Cell, Multiple Directions, Multiple Rays
TEST_F(FluxSolverTest, SingleCellMultipleDirectionsMultipleRays) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with two directions: +x and +y
    std::vector<Direction> predefined_directions;
    
    // Direction 0: +x
    Direction dir0;
    dir0.mu = 0.0;
    dir0.phi = 0.0;
    dir0.weight = 1.0;
    predefined_directions.push_back(dir0);
    
    // Direction 1: +y
    Direction dir1;
    dir1.mu = 0.0;
    dir1.phi = M_PI / 2.0;
    dir1.weight = 1.0;
    predefined_directions.push_back(dir1);
    
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create TrackingData with two rays: one in +x, one in +y
    Vector3D dir_vector0 = directionVector(dir0.mu, dir0.phi);
    Vector3D dir_vector1 = directionVector(dir1.mu, dir1.phi);
    
    TrackingData ray1 = createSingleRay(0, dir_vector0);
    TrackingData ray2 = createSingleRay(1, dir_vector1);
    std::vector<TrackingData> tracking_data = { ray1, ray2 };
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();
    
    // Expected calculations for each ray:
    // psi_in = 0
    // psi_out = (1.0 / 1.0 - 0) * (1 - e^{-1}) ≈ 0.6321205588
    // line_avg_flux = (1.0 / 1.0) - (0.6321205588 - 0) / (1.0 * 1.0) ≈ 0.3678794412
    // flux = line_avg_flux * L_k = 0.3678794412 * 1.0 = 0.3678794412
    // weight = L_k = 1.0
    // Normalized flux = flux / weight = 0.3678794412 / 1.0 = 0.3678794412
    
    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 2) << "There should be flux data for 2 directions";
    
    double expected_flux = 0.3678794412;
    
    EXPECT_NEAR(flux_data[0][0].flux, expected_flux, 1e-6) << "Flux for +x direction should match expected value";
    EXPECT_NEAR(flux_data[0][0].weight, 1.0, 1e-6) << "Weight for +x direction should be equal to L_k (1.0)";
    
    EXPECT_NEAR(flux_data[0][1].flux, expected_flux, 1e-6) << "Flux for +y direction should match expected value";
    EXPECT_NEAR(flux_data[0][1].weight, 1.0, 1e-6) << "Weight for +y direction should be equal to L_k (1.0)";
}

// Test case: Multiple Cells, Single Direction, Multiple Rays
TEST_F(FluxSolverTest, MultipleCellsSingleDirectionMultipleRays) {
    MeshHandler mesh;
    ASSERT_TRUE(setupTwoCellMesh(mesh)) << "Failed to setup two cell mesh";

    // Initialize Field with Q_k for both cells
    Field field;
    ASSERT_TRUE(setupTwoCellField(field)) << "Failed to setup two cell field";

    // Setup face connectivity for two cells
    ASSERT_TRUE(setupTwoCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";

    // Initialize AngularQuadrature with one direction: +x
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);

    // Create TrackingData with two rays in +x direction traversing Cell 0 and Cell 1
    Vector3D dir_vector = directionVector(dir.mu, dir.phi);
    
    // Ray 1
    TrackingData ray1 = createSingleRay(0, dir_vector, 0, 1.0, Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    // Add second CellTrace traversing Cell1
    CellTrace trace1_1;
    trace1_1.cell_id = 1;
    trace1_1.time_spent = 1.0;
    trace1_1.start_point = Vector3D(1.0, 0.0, 0.0); // Entry to Cell1
    trace1_1.end_point = Vector3D(2.0, 0.0, 0.0);   // Exit from Cell1
    ray1.cell_traces.push_back(trace1_1);

    // Ray 2
    TrackingData ray2 = createSingleRay(1, dir_vector, 0, 1.0, Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    // Add second CellTrace traversing Cell1
    CellTrace trace2_1;
    trace2_1.cell_id = 1;
    trace2_1.time_spent = 1.0;
    trace2_1.start_point = Vector3D(1.0, 0.0, 0.0); // Entry to Cell1
    trace2_1.end_point = Vector3D(2.0, 0.0, 0.0);   // Exit from Cell1
    ray2.cell_traces.push_back(trace2_1);

    std::vector<TrackingData> tracking_data = { ray1, ray2 };

    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);

    // Compute flux
    flux_solver.computeFlux();

    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();

    // Expected calculations:
    // For each ray traversing two cells in a given direction:
    // Cell 0:
    //   psi_in = 0
    //   psi_out = (1.0 / 1.0 - 0) * (1 - e^{-1}) ≈ 0.6321205588
    //   line_avg_flux = (1.0 / 1.0) - (0.6321205588 - 0) / (1.0 * 1.0) ≈ 0.3678794412
    //   flux = 0.3678794412 * 1.0 = 0.3678794412
    //   weight = 1.0
    //   psi_in for Cell 1 = 0.6321205588

    //   For two rays, total flux for Cell0, Direction0: 0.3678794412 * 2 = 0.7357588824
    //   Total weight for Cell0, Direction0: 1.0 * 2 = 2.0
    //   Normalized flux: 0.7357588824 / 2.0 = 0.3678794412

    // Cell1:
    //   psi_in = 0.6321205588
    //   psi_out = 0.6321205588 * exp(-1.0) + 2.0 * (1.0 - exp(-1.0)) ≈ 1.4967852755919449
    //   line_avg_flux = (2.0 / 1.0) - (1.4967852755919449 - 0.6321205588) / (1.0 * 1.0) ≈ 1.1353352832366128
    //   flux = 1.1353352832366128 * 1.0 = 1.1353352832366128
    //   weight = 1.0

    //   For two rays, total flux for Cell1, Direction0: 1.1353352832366128 * 2 = 2.2706705664732256
    //   Total weight for Cell1, Direction0: 1.0 * 2 = 2.0
    //   Normalized flux: 2.2706705664732256 / 2.0 = 1.1353352832366128

    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 2) << "There should be flux data for 2 cells";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction per cell";
    ASSERT_EQ(flux_data[1].size(), 1) << "There should be flux data for 1 direction per cell";

    double expected_flux_cell0_dir0 = 0.3678794412;
    double expected_flux_cell1_dir0 = 1.1353352832366128;
    EXPECT_NEAR(flux_data[0][0].flux, expected_flux_cell0_dir0, 1e-6) << "Flux for Cell 0, Direction 0 should match expected value";
    EXPECT_NEAR(flux_data[0][0].weight, 2.0, 1e-6) << "Weight for Cell 0, Direction 0 should be 2.0";

    EXPECT_NEAR(flux_data[1][0].flux, expected_flux_cell1_dir0, 1e-6) << "Flux for Cell 1, Direction 0 should match expected value";
    EXPECT_NEAR(flux_data[1][0].weight, 2.0, 1e-6) << "Weight for Cell 1, Direction 0 should be 2.0";
}

// Test case: Multiple Cells, Single Direction, Multiple Rays
TEST_F(FluxSolverTest, MultipleCellsSingleDirectionMultipleRaysDifferentLengths) {
    MeshHandler mesh;
    ASSERT_TRUE(setupTwoCellMesh(mesh)) << "Failed to setup two cell mesh";

    // Initialize Field with Q_k for both cells
    Field field;
    ASSERT_TRUE(setupTwoCellField(field)) << "Failed to setup two cell field";

    // Setup face connectivity for two cells
    ASSERT_TRUE(setupTwoCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";

    // Initialize AngularQuadrature with one direction: +x
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);

    // Create TrackingData with two rays in +x direction traversing Cell 0 and Cell 1
    Vector3D dir_vector = directionVector(dir.mu, dir.phi);
    
    // Ray 1
    TrackingData ray1 = createSingleRay(0, dir_vector, 0, 2.0, Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    // Add second CellTrace traversing Cell1
    CellTrace trace1_1;
    trace1_1.cell_id = 1;
    trace1_1.time_spent = 1.0;
    trace1_1.start_point = Vector3D(1.0, 0.0, 0.0); // Entry to Cell1
    trace1_1.end_point = Vector3D(2.0, 0.0, 0.0);   // Exit from Cell1
    ray1.cell_traces.push_back(trace1_1);

    // Ray 2
    TrackingData ray2 = createSingleRay(1, dir_vector, 0, 1.0, Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    // Add second CellTrace traversing Cell1
    CellTrace trace2_1;
    trace2_1.cell_id = 1;
    trace2_1.time_spent = 3.0;
    trace2_1.start_point = Vector3D(1.0, 0.0, 0.0); // Entry to Cell1
    trace2_1.end_point = Vector3D(2.0, 0.0, 0.0);   // Exit from Cell1
    ray2.cell_traces.push_back(trace2_1);

    std::vector<TrackingData> tracking_data = { ray1, ray2 };

    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);

    // Compute flux
    flux_solver.computeFlux();

    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();


    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 2) << "There should be flux data for 2 cells";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction per cell";
    ASSERT_EQ(flux_data[1].size(), 1) << "There should be flux data for 1 direction per cell";

    double expected_flux_cell0_dir0 = 0.501071574802685;
    double expected_flux_cell1_dir0 = 1.4956386230969625;
    EXPECT_NEAR(flux_data[0][0].flux, expected_flux_cell0_dir0, 1e-6) << "Flux for Cell 0, Direction 0 should match expected value";
    EXPECT_NEAR(flux_data[0][0].weight, 3.0, 1e-6) << "Weight for Cell 0, Direction 0 should be 2.0";

    EXPECT_NEAR(flux_data[1][0].flux, expected_flux_cell1_dir0, 1e-6) << "Flux for Cell 1, Direction 0 should match expected value";
    EXPECT_NEAR(flux_data[1][0].weight, 4.0, 1e-6) << "Weight for Cell 1, Direction 0 should be 2.0";
}

// Test case: Rays with Invalid Direction
TEST_F(FluxSolverTest, RaysWithInvalidDirection) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with one valid direction
    std::vector<Direction> predefined_directions;
    Direction valid_dir;
    valid_dir.mu = 0.0;
    valid_dir.phi = 0.0;
    valid_dir.weight = 1.0;
    predefined_directions.push_back(valid_dir);
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create TrackingData with two rays: one valid, one invalid
    // Valid Ray
    Vector3D dir_vector0 = directionVector(valid_dir.mu, valid_dir.phi);
    TrackingData ray1 = createSingleRay(0, dir_vector0);
    
    // Invalid Ray (direction not in angular quadrature)
    Vector3D invalid_dir_vector = directionVector(0.0, M_PI / 2.0); // +y direction not defined
    TrackingData ray2 = createSingleRay(1, invalid_dir_vector);
    
    std::vector<TrackingData> tracking_data = { ray1, ray2 };
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();
    
    // Expected calculations:
    // Only ray1 is valid and contributes to flux_data_[0][0]
    // flux = 0.3678794412
    // weight = 1.0
    
    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction";
    
    EXPECT_NEAR(flux_data[0][0].flux, 0.3678794412, 1e-6) << "Flux should match the expected value from the valid ray";
    EXPECT_NEAR(flux_data[0][0].weight, 1.0, 1e-6) << "Weight should be equal to L_k (1.0)";
}

// Test case: Rays with Invalid Cell IDs
TEST_F(FluxSolverTest, RaysWithInvalidCellIDs) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with one direction
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create TrackingData with two rays:
    // - One traversing a valid cell
    // - One traversing an invalid cell (e.g., cell_id = 10)
    
    // Valid Ray
    Vector3D valid_dir_vector = directionVector(dir.mu, dir.phi);
    TrackingData ray1 = createSingleRay(0, valid_dir_vector);
    
    // Invalid Ray
    Vector3D invalid_dir_vector = directionVector(dir.mu, dir.phi);
    TrackingData ray2 = createSingleRay(1, invalid_dir_vector);
    // Modify CellTrace to have an invalid cell_id
    ray2.cell_traces[0].cell_id = 10; // Assuming only cell_id = 0 is valid
    
    std::vector<TrackingData> tracking_data = { ray1, ray2 };
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();
    
    // Expected calculations:
    // Only ray1 is valid and contributes to flux_data_[0][0]
    // flux = 0.3678794412
    // weight = 1.0
    
    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction";
    
    EXPECT_NEAR(flux_data[0][0].flux, 0.3678794412, 1e-6) << "Flux should match the expected value from the valid ray";
    EXPECT_NEAR(flux_data[0][0].weight, 1.0, 1e-6) << "Weight should be equal to L_k (1.0)";
}

// Test case: Multiple Cells, Multiple Directions, Multiple Rays
TEST_F(FluxSolverTest, MultipleCellsMultipleDirectionsMultipleRays) {
    MeshHandler mesh;
    ASSERT_TRUE(setupTwoCellMesh(mesh)) << "Failed to setup two cell mesh";
    Field field;
    ASSERT_TRUE(setupTwoCellField(field)) << "Failed to setup two cell field";
    ASSERT_TRUE(setupTwoCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    

    // Initialize AngularQuadrature with two directions: +x and +y
    std::vector<Direction> predefined_directions;
    
    // Direction 0: +x
    Direction dir0;
    dir0.mu = 0.0;
    dir0.phi = 0.0;
    dir0.weight = 1.0;
    predefined_directions.push_back(dir0);
    
    // Direction 1: +y
    Direction dir1;
    dir1.mu = 0.0;
    dir1.phi = M_PI / 2.0;
    dir1.weight = 1.0;
    predefined_directions.push_back(dir1);
    
    AngularQuadrature angular_quadrature(predefined_directions);

    // Create TrackingData with multiple rays:
    // - Two rays in +x direction traversing Cell 0 and Cell 1
    // - Two rays in +y direction traversing Cell 0 and Cell 1
    std::vector<TrackingData> tracking_data;

    // Rays in +x direction
    for(int i = 0; i < 2; ++i) {
        Vector3D dir_vector = directionVector(dir0.mu, dir0.phi);
        TrackingData ray = createSingleRay(i, dir_vector);
        // Add second CellTrace traversing Cell1
        CellTrace trace1;
        trace1.cell_id = 1;
        trace1.time_spent = 1.0;
        trace1.start_point = Vector3D(1.0, 0.0, 0.0); // Entry to Cell1
        trace1.end_point = Vector3D(2.0, 0.0, 0.0);   // Exit from Cell1
        ray.cell_traces.push_back(trace1);
        tracking_data.push_back(ray);
    }

    // Rays in +y direction
    for(int i = 2; i < 4; ++i) {
        Vector3D dir_vector = directionVector(dir1.mu, dir1.phi);
        TrackingData ray = createSingleRay(i, dir_vector);
        // Add second CellTrace traversing Cell1
        CellTrace trace1;
        trace1.cell_id = 1;
        trace1.time_spent = 1.0;
        trace1.start_point = Vector3D(0.0, 1.0, 0.0); // Entry to Cell1
        trace1.end_point = Vector3D(0.0, 2.0, 0.0);   // Exit from Cell1
        ray.cell_traces.push_back(trace1);
        tracking_data.push_back(ray);
    }

    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);

    // Compute flux
    flux_solver.computeFlux();

    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();

    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 2) << "There should be flux data for 2 cells";
    ASSERT_EQ(flux_data[0].size(), 2) << "There should be flux data for 2 directions per cell";
    ASSERT_EQ(flux_data[1].size(), 2) << "There should be flux data for 2 directions per cell";

    double expected_flux_cell0_dir0 = 0.3678794412;
    double expected_flux_cell1_dir0 = 1.1353352832366128;

    double expected_flux_cell0_dir1 = 0.3678794412;
    double expected_flux_cell1_dir1 = 1.1353352832366128;

    EXPECT_NEAR(flux_data[0][0].flux, expected_flux_cell0_dir0, 1e-6) << "Flux for Cell 0, Direction 0 should match expected value";
    EXPECT_NEAR(flux_data[0][0].weight, 2.0, 1e-6) << "Weight for Cell 0, Direction 0 should be 2.0";

    EXPECT_NEAR(flux_data[1][0].flux, expected_flux_cell1_dir0, 1e-6) << "Flux for Cell 1, Direction 0 should match expected value";
    EXPECT_NEAR(flux_data[1][0].weight, 2.0, 1e-6) << "Weight for Cell 1, Direction 0 should be 2.0";

    EXPECT_NEAR(flux_data[0][1].flux, expected_flux_cell0_dir1, 1e-6) << "Flux for Cell 0, Direction 1 should match expected value";
    EXPECT_NEAR(flux_data[0][1].weight, 2.0, 1e-6) << "Weight for Cell 0, Direction 1 should be 2.0";

    EXPECT_NEAR(flux_data[1][1].flux, expected_flux_cell1_dir1, 1e-6) << "Flux for Cell 1, Direction 1 should match expected value";
    EXPECT_NEAR(flux_data[1][1].weight, 2.0, 1e-6) << "Weight for Cell 1, Direction 1 should be 2.0";
}

// Test case: No Rays
TEST_F(FluxSolverTest, NoRays) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with one direction
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create empty TrackingData
    std::vector<TrackingData> tracking_data = {};
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();
    
    // Expect flux_data to be all zeros
    ASSERT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction";
    
    EXPECT_DOUBLE_EQ(flux_data[0][0].flux, 0.0) << "Flux should be zero when no rays are present";
    EXPECT_DOUBLE_EQ(flux_data[0][0].weight, 0.0) << "Weight should be zero when no rays are present";
}

// Test case: Rays with Zero Path Length
TEST_F(FluxSolverTest, RaysWithZeroPathLength) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with one direction
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create TrackingData with a single ray having zero path length
    Vector3D dir_vector = directionVector(dir.mu, dir.phi);
    TrackingData ray = createSingleRay(0, dir_vector);
    // Modify CellTrace to have L_k = 0
    ray.cell_traces[0].time_spent = 0.0;
    
    std::vector<TrackingData> tracking_data = { ray };
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();
    
    // Expected calculations:
    // psi_in = 0
    // psi_out = (1.0 / 1.0 - 0) * (1 - e^{-1 * 0}) = 1.0 * (1 - 1) = 0
    // line_avg_flux = (1.0 / 1.0) - (0 - 0) / (1.0 * 0) = 1.0 - undefined (handle gracefully)
    // Since L_k = 0, line_avg_flux should be set appropriately to avoid division by zero
    // line_avg_flux = 0 because int_{0}^{0} f(x) dx = 0
    
    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction";
    
    EXPECT_NEAR(flux_data[0][0].flux, 0.0, 1e-6) << "Flux should be equal to Q_k / sigma_t when L_k = 0";
    EXPECT_NEAR(flux_data[0][0].weight, 0.0, 1e-6) << "Weight should be zero since L_k = 0";
}

// Test case: Rays Traversing the Same Cell Multiple Times
TEST_F(FluxSolverTest, RaysTraversingSameCellMultipleTimes) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with one direction
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create TrackingData with a single ray traversing the same cell twice
    Vector3D dir_vector = directionVector(dir.mu, dir.phi);
    TrackingData ray = createSingleRay(0, dir_vector, 0, 6.0e-1, Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    
    // Second traversal of Cell 0
    CellTrace trace2;
    trace2.cell_id = 0;
    trace2.time_spent = 1.0e-1;
    trace2.start_point = Vector3D(1.0, 0.0, 0.0);
    trace2.end_point = Vector3D(0.0, 0.0, 0.0);
    
    ray.cell_traces.push_back(trace2);
    
    std::vector<TrackingData> tracking_data = { ray };
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();

    EXPECT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction";
    
    EXPECT_NEAR(flux_data[0][0].flux, 0.28083614827344205, 1e-6) << "Flux should match the expected normalized value";
    EXPECT_NEAR(flux_data[0][0].weight, 7e-1, 1e-6) << "Weight should be equal to total L_k (2.0)";
}

// Test case: Flux Normalization
TEST_F(FluxSolverTest, FluxNormalization) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with one direction
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create TrackingData with three rays traversing Cell 0
    std::vector<TrackingData> tracking_data;
    for(int i = 0; i < 3; ++i) {
        Vector3D dir_vector = directionVector(dir.mu, dir.phi);
        TrackingData ray = createSingleRay(i, dir_vector);
        tracking_data.push_back(ray);
    }
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();
    
    // Expected calculations:
    // Each ray:
    //   psi_in = 0
    //   psi_out = (1.0 / 1.0 - 0) * (1 - e^{-1}) ≈ 0.6321205588
    //   line_avg_flux = (1.0 / 1.0) - (0.6321205588 - 0) / (1.0 * 1.0) ≈ 0.3678794412
    //   flux = 0.3678794412 * 1.0 = 0.3678794412
    //   weight = 1.0
    //   psi_in for next traversal = 0.6321205588

    // Total flux: 0.3678794412 * 3 ≈ 1.1036383236
    // Total weight: 1.0 * 3 = 3.0
    // Normalized flux: 1.1036383236 / 3.0 ≈ 0.3678794412

    // Verify flux_data_[cell][direction]
    ASSERT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction";
    
    EXPECT_NEAR(flux_data[0][0].flux, 0.3678794412, 1e-6) << "Normalized flux should match expected value";
    EXPECT_NEAR(flux_data[0][0].weight, 3.0, 1e-6) << "Weight should be equal to total L_k (3.0)";
}

// Test case: Rays with Zero Direction Vector
TEST_F(FluxSolverTest, RaysWithZeroDirection) {
    MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellMesh(mesh)) << "Failed to setup single cell mesh";
    
    Field field;
    ASSERT_TRUE(setupSingleCellField(field)) << "Failed to setup single cell field";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    // Initialize AngularQuadrature with one direction
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 1.0;
    predefined_directions.push_back(dir);
    AngularQuadrature angular_quadrature(predefined_directions);
    
    // Create TrackingData with a single ray having zero direction
    TrackingData ray;
    ray.ray_id = 0;
    ray.direction = Vector3D(0.0, 0.0, 0.0); // Zero vector
    
    // CellTrace remains the same
    // This should ideally be skipped or handled gracefully
    
    std::vector<TrackingData> tracking_data = { ray };
    
    // Initialize FluxSolver
    double sigma_t = 1.0; // Total cross section
    FluxSolver flux_solver(mesh, field, tracking_data, angular_quadrature, sigma_t);
    
    // Compute flux
    flux_solver.computeFlux();
    
    // Retrieve flux data
    const auto& flux_data = flux_solver.getFluxData();
    
    // Expect flux_data to be all zeros since direction is invalid
    ASSERT_EQ(flux_data.size(), 1) << "There should be flux data for 1 cell";
    ASSERT_EQ(flux_data[0].size(), 1) << "There should be flux data for 1 direction";
    
    EXPECT_DOUBLE_EQ(flux_data[0][0].flux, 0.0) << "Flux should be zero for zero direction";
    EXPECT_DOUBLE_EQ(flux_data[0][0].weight, 0.0) << "Weight should be zero for zero direction";
}
