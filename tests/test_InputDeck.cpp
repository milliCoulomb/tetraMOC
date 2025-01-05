// tests/test_InputDeck.cpp

#include "InputDeck.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <cstdio>
#include <yaml-cpp/yaml.h>

// Helper function to create temporary YAML files
bool createTempYAML(const std::string& filename, const std::string& content) {
    std::ofstream file(filename);
    if (!file.is_open()) return false;
    file << content;
    file.close();
    return true;
}

// Test parsing a valid InputDeck
TEST(InputDeckTest, ParseValidInputDeck) {
    std::string filename = "temp_valid_inputdeck.yaml";
    std::string yaml_content = R"(
mesh:
  nodes: "path/to/nodes.dat"
  cells: "path/to/cells.dat"
  faces: "path/to/faces.dat"

cross_sections:
  data_files:
    - "path/to/cs1.dat"
    - "path/to/cs2.dat"

angular_quadrature:
  ntheta: 10
  nphi: 20

solver_parameters:
  multi_group_max_iterations: 1000
  multi_group_tolerance: 1e-7
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  max_power_iterations: 100
  rays_per_face: 8
  max_ray_length: 1000
  use_half_hemisphere: false

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"

logging:
  level: "INFO"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    InputDeck deck;
    EXPECT_NO_THROW(deck = parser.parse(filename));

    // Verify mesh paths
    EXPECT_EQ(deck.mesh.nodes, "path/to/nodes.dat");
    EXPECT_EQ(deck.mesh.cells, "path/to/cells.dat");
    EXPECT_EQ(deck.mesh.faces, "path/to/faces.dat");

    // Verify cross_sections data_files
    ASSERT_EQ(deck.cross_sections.data_files.size(), 2);
    EXPECT_EQ(deck.cross_sections.data_files[0], "path/to/cs1.dat");
    EXPECT_EQ(deck.cross_sections.data_files[1], "path/to/cs2.dat");

    // Verify angular_quadrature parameters
    EXPECT_EQ(deck.angular_quadrature_parameters.ntheta, 10);
    EXPECT_EQ(deck.angular_quadrature_parameters.nphi, 20);

    // Verify solver_parameters
    EXPECT_EQ(deck.solver_parameters.multi_group_max_iterations, 1000);
    EXPECT_DOUBLE_EQ(deck.solver_parameters.multi_group_tolerance, 1e-7);
    EXPECT_EQ(deck.solver_parameters.one_group_max_iterations, 500);
    EXPECT_DOUBLE_EQ(deck.solver_parameters.one_group_tolerance, 1e-7);
    EXPECT_DOUBLE_EQ(deck.solver_parameters.fission_source_tolerance, 1e-7);
    EXPECT_DOUBLE_EQ(deck.solver_parameters.keff_tolerance, 1e-6);

    // Verify output settings
    EXPECT_EQ(deck.output.flux_output_file, "output/flux.dat");
    EXPECT_EQ(deck.output.k_eff_output_file, "output/k_eff.dat");

    // Verify logging settings
    EXPECT_EQ(deck.logging.level, "INFO");

    // Clean up
    std::remove(filename.c_str());
}

// Test parsing InputDeck with missing 'mesh' section
TEST(InputDeckTest, ParseMissingMeshSection) {
    std::string filename = "temp_missing_mesh.yaml";
    std::string yaml_content = R"(
cross_sections:
  data_files:
    - "path/to/cs1.dat"

angular_quadrature:
  ntheta: 10
  nphi: 20

solver_parameters:
  multi_group_max_iterations: 1000
  multi_group_tolerance: 1e-7
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  max_power_iterations: 100
  rays_per_face: 8
  max_ray_length: 1000
  use_half_hemisphere: false

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"

logging:
  level: "INFO"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    EXPECT_THROW(parser.parse(filename), std::runtime_error);

    // Clean up
    std::remove(filename.c_str());
}

// Test parsing InputDeck with missing 'cross_sections' section
TEST(InputDeckTest, ParseMissingCrossSectionsSection) {
    std::string filename = "temp_missing_cross_sections.yaml";
    std::string yaml_content = R"(
mesh:
  nodes: "path/to/nodes.dat"
  cells: "path/to/cells.dat"
  faces: "path/to/faces.dat"

angular_quadrature:
  ntheta: 10
  nphi: 20

solver_parameters:
  multi_group_max_iterations: 1000
  multi_group_tolerance: 1e-7
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  max_power_iterations: 100
  rays_per_face: 8
  max_ray_length: 1000
  use_half_hemisphere: false

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"

logging:
  level: "INFO"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    EXPECT_THROW(parser.parse(filename), std::runtime_error);

    // Clean up
    std::remove(filename.c_str());
}

// Test parsing InputDeck with invalid 'angular_quadrature' data
TEST(InputDeckTest, ParseInvalidAngularQuadratureData) {
    std::string filename = "temp_invalid_angular_quadrature.yaml";
    std::string yaml_content = R"(
mesh:
  nodes: "path/to/nodes.dat"
  cells: "path/to/cells.dat"
  faces: "path/to/faces.dat"

cross_sections:
  data_files:
    - "path/to/cs1.dat"

angular_quadrature:
  ntheta: "ten"  # Invalid type
  nphi: 20

solver_parameters:
  multi_group_max_iterations: 1000
  multi_group_tolerance: 1e-7
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  max_power_iterations: 100
  rays_per_face: 8
  max_ray_length: 1000
  use_half_hemisphere: false

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"

logging:
  level: "INFO"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    EXPECT_THROW(parser.parse(filename), YAML::BadConversion);

    // Clean up
    std::remove(filename.c_str());
}

// Test parsing InputDeck with missing 'solver_parameters' section
TEST(InputDeckTest, ParseMissingSolverParametersSection) {
    std::string filename = "temp_missing_solver_parameters.yaml";
    std::string yaml_content = R"(
mesh:
  nodes: "path/to/nodes.dat"
  cells: "path/to/cells.dat"
  faces: "path/to/faces.dat"

cross_sections:
  data_files:
    - "path/to/cs1.dat"

angular_quadrature:
  ntheta: 10
  nphi: 20

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"

logging:
  level: "INFO"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    EXPECT_THROW(parser.parse(filename), std::runtime_error);

    // Clean up
    std::remove(filename.c_str());
}

// Test parsing InputDeck with invalid 'solver_parameters' values
TEST(InputDeckTest, ParseInvalidSolverParametersValues) {
    std::string filename = "temp_invalid_solver_parameters.yaml";
    std::string yaml_content = R"(
mesh:
  nodes: "path/to/nodes.dat"
  cells: "path/to/cells.dat"
  faces: "path/to/faces.dat"

cross_sections:
  data_files:
    - "path/to/cs1.dat"

angular_quadrature:
  ntheta: 10
  nphi: 20

solver_parameters:
  multi_group_max_iterations: "a thousand"  # Invalid type
  multi_group_tolerance: -1e-7            # Invalid value
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  max_power_iterations: 100
  rays_per_face: 8
  max_ray_length: 1000
  use_half_hemisphere: false

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"

logging:
  level: "INFO"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    EXPECT_THROW(parser.parse(filename), YAML::BadConversion);

    // Clean up
    std::remove(filename.c_str());
}

// Test parsing InputDeck with missing 'output' section
TEST(InputDeckTest, ParseMissingOutputSection) {
    std::string filename = "temp_missing_output.yaml";
    std::string yaml_content = R"(
mesh:
  nodes: "path/to/nodes.dat"
  cells: "path/to/cells.dat"
  faces: "path/to/faces.dat"

cross_sections:
  data_files:
    - "path/to/cs1.dat"

angular_quadrature:
  ntheta: 10
  nphi: 20

solver_parameters:
  multi_group_max_iterations: 1000
  multi_group_tolerance: 1e-7
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  max_power_iterations: 100
  rays_per_face: 8
  max_ray_length: 1000
  use_half_hemisphere: false

logging:
  level: "INFO"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    EXPECT_THROW(parser.parse(filename), std::runtime_error);

    // Clean up
    std::remove(filename.c_str());
}

// Test parsing InputDeck with missing 'logging' section
TEST(InputDeckTest, ParseMissingLoggingSection) {
    std::string filename = "temp_missing_logging.yaml";
    std::string yaml_content = R"(
mesh:
  nodes: "path/to/nodes.dat"
  cells: "path/to/cells.dat"
  faces: "path/to/faces.dat"

cross_sections:
  data_files:
    - "path/to/cs1.dat"

angular_quadrature:
  ntheta: 10
  nphi: 20

solver_parameters:
  multi_group_max_iterations: 1000
  multi_group_tolerance: 1e-7
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  max_power_iterations: 100
  rays_per_face: 8
  max_ray_length: 1000
  use_half_hemisphere: false

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    EXPECT_THROW(parser.parse(filename), std::runtime_error);

    // Clean up
    std::remove(filename.c_str());
}

// // Test parsing InputDeck with invalid 'output' file paths
// TEST(InputDeckTest, ParseInvalidOutputFilePaths) {
//     std::string filename = "temp_invalid_output.yaml";
//     std::string yaml_content = R"(
// mesh:
//   nodes: "path/to/nodes.dat"
//   cells: "path/to/cells.dat"
//   faces: "path/to/faces.dat"

// cross_sections:
//   data_files:
//     - "path/to/cs1.dat"

// angular_quadrature:
//   ntheta: 10
//   nphi: 20

// solver_parameters:
//   multi_group_max_iterations: 1000
//   multi_group_tolerance: 1e-7
//   one_group_max_iterations: 500
//   one_group_tolerance: 1e-7
//   fission_source_tolerance: 1e-7
//   keff_tolerance: 1e-6

// output:
//   flux_output_file: 85886868
//   k_eff_output_file: "output/k_eff.dat"

// logging:
//   level: "INFO"
//   log_file: "logs/solver.log"
// )";

//     ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

//     InputDeckParser parser;
//     EXPECT_THROW(parser.parse(filename), YAML::BadConversion);

//     // Clean up
//     std::remove(filename.c_str());
// }

// Test parsing InputDeck with extra unexpected sections
TEST(InputDeckTest, ParseWithExtraSections) {
    std::string filename = "temp_extra_sections.yaml";
    std::string yaml_content = R"(
mesh:
  nodes: "path/to/nodes.dat"
  cells: "path/to/cells.dat"
  faces: "path/to/faces.dat"

cross_sections:
  data_files:
    - "path/to/cs1.dat"

angular_quadrature:
  ntheta: 10
  nphi: 20

solver_parameters:
  multi_group_max_iterations: 1000
  multi_group_tolerance: 1e-7
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  max_power_iterations: 100
  rays_per_face: 8
  max_ray_length: 1000
  use_half_hemisphere: false

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"

logging:
  level: "INFO"

extra_section:
  key: "value"
)";

    ASSERT_TRUE(createTempYAML(filename, yaml_content)) << "Failed to create temporary YAML file.";

    InputDeckParser parser;
    InputDeck deck;
    EXPECT_NO_THROW(deck = parser.parse(filename));

    // Verify that extra_section is ignored or handled as per implementation
    // Assuming it's ignored
    EXPECT_EQ(deck.mesh.nodes, "path/to/nodes.dat");
    // Additional checks can be added as needed

    // Clean up
    std::remove(filename.c_str());
}