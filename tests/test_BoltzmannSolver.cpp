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
#include "InputHandler.hpp"
#include "BoltzmannSolver.hpp"
#include "RayTracerManager.hpp"

#include <fstream>
#include <cstdio> // For std::remove
#include <cmath>
#include <vector>
#include <string>

// Test Fixture for FluxSolver
class BoltzmannSolverTest : public ::testing::Test {
protected:
    // Temporary file names
    std::string nodes_file = "temp_nodes_test_FluxSolver.txt";
    std::string cells_file = "temp_cells_test_FluxSolver.txt";
    std::string faces_file = "temp_faces_test_FluxSolver.txt";
    std::string nuc_data_file = "temp_nuc_data_test_FluxSolver.txt";
    std::string field_file = "temp_field_test_FluxSolver.txt";

    // Clean up temporary files after each test
    void TearDown() override {
        std::remove(nodes_file.c_str());
        std::remove(cells_file.c_str());
        std::remove(faces_file.c_str());
        std::remove(nuc_data_file.c_str());
        std::remove(field_file.c_str());
    }

    // Helper function to create a simple mesh with a single tetrahedron
    bool setupSingleCellInfiniteMesh(MeshHandler& mesh) {
        // Define nodes content with node IDs
        std::string nodes_content = "4\n"
                                    "0 0.0 0.0 0.0\n" // Node 0
                                    "1 1.0e10 0.0 0.0\n" // Node 1
                                    "2 0.0 1.0e10 0.0\n" // Node 2
                                    "3 0.0 0.0 1.0e10\n"; // Node 3

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
    bool setupTwoCellsInfiniteMesh(MeshHandler& mesh) {
        // Define nodes content with node IDs
        std::string nodes_content = "5\n"
                                    "0 0.0 0.0 0.0\n" // Node 0
                                    "1 1.0e20 0.0 0.0\n" // Node 1
                                    "2 0.0 1.0e20 0.0\n" // Node 2
                                    "3 0.0 0.0 1.0e20\n" // Node 3
                                    "4 1.0e20 1.0e20 1.0e20\n"; // Node 4

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

    // Helper function to setup a simple field
    bool setupSimpleFieldOneCell(Field& field) {
        // Define vector fields for two cells
        // Cell 0: Velocity towards positive x-axis
        // Cell 1: Velocity towards positive y-axis
        std::string field_content = "1\n" // Number of vectors
                                       "0.0 1.0 0.0\n"; // Vector for cell 0

        // Create temporary field.txt
        if(!createTempFile(field_file, field_content)) return false;
        if(!field.loadVectorField(field_file)) return false;

        return true;
    }

    // Helper function to setup a simple field with two source terms
    bool setupSimpleFieldTwoCells(Field& field) {
        // Define vector fields for two cells
        // Cell 0: Velocity towards positive x-axis
        // Cell 1: Velocity towards positive y-axis
        std::string field_content = "2\n" // Number of vectors
                                       "0.0 1.0 0.0\n" // Vector for cell 0
                                       "1.0 0.0 0.0\n"; // Vector for cell 1

        // Create temporary field.txt
        if(!createTempFile(field_file, field_content)) return false;
        if(!field.loadVectorField(field_file)) return false;

        return true;
    }

    // Helper function to setup a simple field with a single source term
    std::vector<double> setupSingleCellField() {
        // Define scalar fields content
        // returns a std::vector<double> with a single source term, 1.0
        std::vector<double> source = {1.0};
        return source;
    }

    // Helper function to setup a simple field with two source terms
    std::vector<double> setupTwoCellField() {
        // Define scalar fields content
        // returns a std::vector<double> with two source terms, 1.0 and 2.0
        std::vector<double> source = {1.0, 2.0};
        return source;
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

    // Helper function to create input handler
    bool setupInputHandler(InputHandler& input_handler) {
        // Define material data content
        std::string nuc_data_content = 
            "1\n" // Number of materials
            "1.0 0.01 0.1 2.0 0.95 0.05\n"; // Material 0: total_xs, fission_xs, scattering_xs, multiplicity, fission_spectrum, delayed_spectrum
        
        // Create temporary nuc_data.txt
        if(!createTempFile(nuc_data_file, nuc_data_content)) return false;
        if(!input_handler.loadData(nuc_data_file)) return false;

        return true;
    }
};

TEST_F(BoltzmannSolverTest, SingleCellTwoRaysInfiniteOneGroup) {
     MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellInfiniteMesh(mesh)) << "Failed to setup single cell mesh";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";

    InputHandler input_handler;
    ASSERT_TRUE(setupInputHandler(input_handler)) << "Failed to setup input handler";
    
    // Initialize AngularQuadrature with one direction
    std::vector<Direction> predefined_directions;
    Direction dir;
    dir.mu = 0.0;
    dir.phi = 0.0;
    dir.weight = 3.0; // Weight for direction 0
    predefined_directions.push_back(dir);
    Direction dir2;
    dir2.mu = 0.0;
    dir2.phi = M_PI / 2.0;
    dir2.weight = 1.0; // Weight for direction 1
    predefined_directions.push_back(dir2);
    AngularQuadrature angular_quadrature(predefined_directions);
    // test if size of predefined_directions is 2
    ASSERT_EQ(angular_quadrature.getDirections().size(), 2) << "There should be 2 directions";
    // test if sum of weights is 4
    ASSERT_EQ(angular_quadrature.getTotalWeight(), 4.0) << "Sum of weights should be 4.0";

    // Create TrackingData with two rays: one in each direction
    std::vector<TrackingData> tracking_data;
    Vector3D dir_vector = directionVector(0.0, 0.0);
    // we set extremely large path length to simulate infinite path length
    TrackingData ray1 = createSingleRay(0, dir_vector, 0, 1.0e10, Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
    tracking_data.push_back(ray1);
    Vector3D dir_vector2 = directionVector(0.0, M_PI / 2.0);
    TrackingData ray2 = createSingleRay(1, dir_vector2, 0, 2.0e10, Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 1.0, 0.0));
    tracking_data.push_back(ray2);
    // test if size of tracking_data is 2
    ASSERT_EQ(tracking_data.size(), 2) << "There should be 2 rays";

    // Initialize FluxSolver
    double sigma_t = input_handler.getEnergyGroupData(0).total_xs; // Total cross section
    double sigma_s = input_handler.getEnergyGroupData(0).scattering_xs; // Scattering cross section
    FluxSolver flux_solver(mesh, tracking_data, angular_quadrature, sigma_t);
    // test if size of flux_data_ is 1
    ASSERT_EQ(flux_solver.getFluxData().size(), 1) << "There should be flux data for 1 cell";
    // test if size of flux_data_[0] is 2
    ASSERT_EQ(flux_solver.getFluxData()[0].size(), 2) << "There should be flux data for 2 directions";

    // create source term
    std::vector<double> source = setupSingleCellField();
    // Compute flux
    flux_solver.computeFlux(source);
    // get flux data
    const auto& flux_data = flux_solver.getFluxData();
    const double expected_flux_dir0 = source[0] / sigma_t;
    const double expected_flux_dir1 = expected_flux_dir0;
    // test if flux_data_[0][0].flux is expected_flux_dir0
    EXPECT_NEAR(flux_data[0][0].flux, expected_flux_dir0, 1e-6) << "Flux for Direction 0 should match Q / sigma_t";
    // test if flux_data_[0][1].flux is expected_flux_dir1
    EXPECT_NEAR(flux_data[0][1].flux, expected_flux_dir1, 1e-6) << "Flux for Direction 1 should match Q / sigma_t";
    // collapse flux
    std::vector<double> collapsed_flux = flux_solver.collapseFlux();
    // test if size of collapsed_flux is 1
    ASSERT_EQ(collapsed_flux.size(), 1) << "There should be 1 collapsed flux value";
    // test if collapsed_flux[0] is the weighted sum of fluxes
    const double expected_collapsed_flux = (expected_flux_dir0 * dir.weight + expected_flux_dir1 * dir2.weight);
    EXPECT_NEAR(collapsed_flux[0], expected_collapsed_flux, 1e-6) << "Collapsed flux should match the weighted sum of fluxes";
    
    // initialiase struct params from BoltzmannSolver
    BoltzmannSolver::SolverParams params;

    BoltzmannSolver bt_solver(input_handler, mesh, tracking_data, angular_quadrature, params);
    std::vector<double> external_source = setupSingleCellField();
    std::vector<double> scalar_flux = bt_solver.solveOneGroupWithSource(external_source, 0, 1e-7);
    const double expected_scalar_flux_infinite = external_source[0] / (sigma_t - sigma_s);
    EXPECT_NEAR(scalar_flux[0], expected_scalar_flux_infinite, 1e-6) << "Scalar flux should match expected value for infinite medium";
}

// Test case: Multiple Cells, Multiple Directions, Multiple Rays
TEST_F(BoltzmannSolverTest, MultipleCellsMultipleDirectionsMultipleRays) {
    MeshHandler mesh;
    ASSERT_TRUE(setupTwoCellsInfiniteMesh(mesh)) << "Failed to setup two cell mesh";
    
    ASSERT_TRUE(setupTwoCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    InputHandler input_handler;
    ASSERT_TRUE(setupInputHandler(input_handler)) << "Failed to setup input handler";

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

    // check the sum of weights
    ASSERT_EQ(angular_quadrature.getTotalWeight(), 2.0) << "Sum of weights should be 2.0";
    // check the number of directions
    ASSERT_EQ(angular_quadrature.getDirections().size(), 2) << "There should be 2 directions";

    // Create TrackingData with multiple rays:
    // - Two rays in +x direction traversing Cell 0 and Cell 1
    // - Two rays in +y direction traversing Cell 0 and Cell 1
    std::vector<TrackingData> tracking_data;

    // Rays in +x direction go though cell 0 first and then cell 1
    for(int i = 0; i < 2; ++i) {
        Vector3D dir_vector = directionVector(dir0.mu, dir0.phi);
        TrackingData ray = createSingleRay(i, dir_vector, 0, 1.0e10, Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 0.0, 0.0));
        // Add second CellTrace traversing Cell1
        CellTrace trace1;
        trace1.cell_id = 1;
        trace1.time_spent = 1.0e10; // infinite path length
        trace1.start_point = Vector3D(1.0, 0.0, 0.0); // Entry to Cell1
        trace1.end_point = Vector3D(2.0, 0.0, 0.0);   // Exit from Cell1
        ray.cell_traces.push_back(trace1);
        tracking_data.push_back(ray);
    }

    // Rays in +y direction go through cell 1 first and then cell 0
    for(int i = 2; i < 4; ++i) {
        Vector3D dir_vector = directionVector(dir1.mu, dir1.phi);
        TrackingData ray = createSingleRay(i, dir_vector, 1, 1.0e10, Vector3D(1.0, 0.0, 0.0), Vector3D(1.0, 1.0, 0.0));
        // Add second CellTrace traversing Cell0
        CellTrace trace1;
        trace1.cell_id = 0;
        trace1.time_spent = 1.0e10;
        trace1.start_point = Vector3D(0.0, 1.0, 0.0); // Entry to Cell1
        trace1.end_point = Vector3D(0.0, 2.0, 0.0);   // Exit from Cell1
        ray.cell_traces.push_back(trace1);
        tracking_data.push_back(ray);
    }

    // Initialize FluxSolver
    double sigma_t = input_handler.getEnergyGroupData(0).total_xs; // Total cross section
    double sigma_s = input_handler.getEnergyGroupData(0).scattering_xs; // Scattering cross section
    FluxSolver flux_solver(mesh, tracking_data, angular_quadrature, sigma_t);
    // check the size of flux_data_
    ASSERT_EQ(flux_solver.getFluxData().size(), 2) << "There should be flux data for 2 cells";
    // initialiase struct params from BoltzmannSolver
    BoltzmannSolver::SolverParams params;

    BoltzmannSolver bt_solver(input_handler, mesh, tracking_data, angular_quadrature, params);
    std::vector<double> external_source = setupTwoCellField();
    std::vector<double> scalar_flux = bt_solver.solveOneGroupWithSource(external_source, 0, 1e-7);
    const double expected_scalar_flux_infinite_cell0 = external_source[0] / (sigma_t - sigma_s);
    const double expected_scalar_flux_infinite_cell1 = external_source[1] / (sigma_t - sigma_s);
    EXPECT_NEAR(scalar_flux[0], expected_scalar_flux_infinite_cell0, 1e-6) << "Scalar flux for Cell 0 should match expected value for infinite medium";
    EXPECT_NEAR(scalar_flux[1], expected_scalar_flux_infinite_cell1, 1e-6) << "Scalar flux for Cell 1 should match expected value for infinite medium";
}

// Test case: Single Cell, true angular quadrature, true RayTracing, infinite medium, one group
TEST_F(BoltzmannSolverTest, SingleCellTrueAngularQuadratureTrueRayTracingInfiniteOneGroup) {
     MeshHandler mesh;
    ASSERT_TRUE(setupSingleCellInfiniteMesh(mesh)) << "Failed to setup single cell mesh";
    
    ASSERT_TRUE(setupSingleCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";

    InputHandler input_handler;
    ASSERT_TRUE(setupInputHandler(input_handler)) << "Failed to setup input handler";

    Field field;
    ASSERT_TRUE(setupSimpleFieldOneCell(field)) << "Failed to setup simple field";
    
    const int num_azimuthal = 4;
    const int num_polar = 2;
    AngularQuadrature angular_quadrature(num_azimuthal, num_polar);
    // test if size of predefined_directions is 2
    ASSERT_EQ(angular_quadrature.getDirections().size(), 8) << "There should be 8 directions";
    // test if sum of weights is 4
    EXPECT_NEAR(angular_quadrature.getTotalWeight(), 4.0 * M_PI, 1e-6) << "Sum of weights should be 4.0 * M_PI";

    bool use_half_quadrature_for_constant = false;
    bool constant_dir_bool = true;
    RayTracerManager manager(mesh, field, angular_quadrature, constant_dir_bool, use_half_quadrature_for_constant);
    // Generate tracking data
    int rays_per_face = 3;
    manager.generateTrackingData(rays_per_face);
    
    // Retrieve tracking data
    const std::vector<TrackingData>& tracking_data = manager.getTrackingData();

    // Initialize FluxSolver
    double sigma_t = input_handler.getEnergyGroupData(0).total_xs; // Total cross section
    double sigma_s = input_handler.getEnergyGroupData(0).scattering_xs; // Scattering cross section
    FluxSolver flux_solver(mesh, tracking_data, angular_quadrature, sigma_t);
    // test if size of flux_data_ is 1
    ASSERT_EQ(flux_solver.getFluxData().size(), 1) << "There should be flux data for 1 cell";
    // test if size of flux_data_[0] is 2
    ASSERT_EQ(flux_solver.getFluxData()[0].size(), 8) << "There should be flux data for 8 directions";

    // create source term
    std::vector<double> source = setupSingleCellField();
    
    // initialiase struct params from BoltzmannSolver
    BoltzmannSolver::SolverParams params;

    BoltzmannSolver bt_solver(input_handler, mesh, tracking_data, angular_quadrature, params);
    std::vector<double> external_source = setupSingleCellField();
    std::vector<double> scalar_flux = bt_solver.solveOneGroupWithSource(external_source, 0, 1e-7);
    const double expected_scalar_flux_infinite = external_source[0] / (sigma_t - sigma_s);
    EXPECT_NEAR(scalar_flux[0], expected_scalar_flux_infinite, 1e-6) << "Scalar flux should match expected value for infinite medium";
}

// Test case: Multiple Cells, true angular quadrature, true RayTracing, infinite medium, one group
TEST_F(BoltzmannSolverTest, MultipleCellsTrueAngularQuadratureTrueRayTracingInfiniteOneGroup) {
    MeshHandler mesh;
    ASSERT_TRUE(setupTwoCellsInfiniteMesh(mesh)) << "Failed to setup two cell mesh";
    
    ASSERT_TRUE(setupTwoCellFaceConnectivity(mesh)) << "Failed to setup face connectivity";
    
    InputHandler input_handler;
    ASSERT_TRUE(setupInputHandler(input_handler)) << "Failed to setup input handler";

    Field field;
    ASSERT_TRUE(setupSimpleFieldTwoCells(field)) << "Failed to setup simple field";
    
    const int num_azimuthal = 2;
    const int num_polar = 4;
    AngularQuadrature angular_quadrature(num_azimuthal, num_polar);
    // test if size of predefined_directions is 2
    ASSERT_EQ(angular_quadrature.getDirections().size(), num_azimuthal * num_polar) << "There should be " << num_azimuthal * num_polar << " directions";
    // test if sum of weights is 4
    EXPECT_NEAR(angular_quadrature.getTotalWeight(), 4.0 * M_PI, 1e-6) << "Sum of weights should be 4.0 * M_PI";

    bool use_half_quadrature_for_constant = false;
    bool constant_dir_bool = true;
    RayTracerManager manager(mesh, field, angular_quadrature, constant_dir_bool, use_half_quadrature_for_constant);
    // Generate tracking data
    int rays_per_face = 5;
    manager.generateTrackingData(rays_per_face);
    
    // Retrieve tracking data
    const std::vector<TrackingData>& tracking_data = manager.getTrackingData();

    // Initialize FluxSolver
    double sigma_t = input_handler.getEnergyGroupData(0).total_xs; // Total cross section
    double sigma_s = input_handler.getEnergyGroupData(0).scattering_xs; // Scattering cross section
    FluxSolver flux_solver(mesh, tracking_data, angular_quadrature, sigma_t);
    // test if size of flux_data_ is 1
    ASSERT_EQ(flux_solver.getFluxData().size(), 2) << "There should be flux data for 2 cells";
    // test if size of flux_data_[0] is 2
    ASSERT_EQ(flux_solver.getFluxData()[0].size(), num_azimuthal * num_polar) << "There should be flux data for " << num_azimuthal * num_polar << " directions";

    // initialiase struct params from BoltzmannSolver
    BoltzmannSolver::SolverParams params;

    BoltzmannSolver bt_solver(input_handler, mesh, tracking_data, angular_quadrature, params);

    std::vector<double> external_source = setupTwoCellField();

    std::vector<double> scalar_flux = bt_solver.solveOneGroupWithSource(external_source, 0, 1e-7);
    const double expected_scalar_flux_infinite_cell0 = external_source[0] / (sigma_t - sigma_s);
    const double expected_scalar_flux_infinite_cell1 = external_source[1] / (sigma_t - sigma_s);

    EXPECT_NEAR(scalar_flux[0], expected_scalar_flux_infinite_cell0, 1e-5) << "Scalar flux for Cell 0 should match expected value for infinite medium";
    EXPECT_NEAR(scalar_flux[1], expected_scalar_flux_infinite_cell1, 1e-5) << "Scalar flux for Cell 1 should match expected value for infinite medium";

}
