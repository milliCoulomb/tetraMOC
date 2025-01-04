// include/InputDeck.hpp
#ifndef INPUT_DECK_HPP
#define INPUT_DECK_HPP

#include <string>
#include <vector>

// Structure to hold mesh file paths
struct MeshPaths {
    std::string nodes;
    std::string cells;
    std::string faces;
};

// Structure to hold cross-section file paths
struct CrossSectionPaths {
    std::vector<std::string> data_files;
};

// Structure to hold angular quadrature parameters (ntheta, nphi)
struct AngularQuadratureParameters {
    int ntheta;
    int nphi;
};

// Structure to hold solver parameters
struct SolverParameters {
    double multi_group_tolerance;
    int multi_group_max_iterations;
    double one_group_tolerance;
    int one_group_max_iterations;
    double fission_source_tolerance;
    double keff_tolerance;
    int max_power_iterations;
    int rays_per_face;
    int max_ray_length;
    bool use_half_hemisphere; // use only the upper hemisphere for ray tracing and symmetrize afterwards
};

// Structure to hold output settings
struct OutputSettings {
    std::string flux_output_file;
    std::string k_eff_output_file;
};

// Structure to hold logging settings
struct LoggingSettings {
    std::string level;
    std::string log_file;
};

// Main structure to hold the entire input deck
struct InputDeck {
    MeshPaths mesh;
    CrossSectionPaths cross_sections;
    AngularQuadratureParameters angular_quadrature_parameters;
    SolverParameters solver_parameters;
    OutputSettings output;
    LoggingSettings logging;
};

class InputDeckParser {
public:
    /**
     * @brief Parses the YAML input deck file and populates the InputDeck structure.
     * 
     * @param filename Path to the YAML input deck file.
     * @return InputDeck Parsed input deck.
     */
    static InputDeck parse(const std::string& filename);
};

#endif // INPUT_DECK_HPP
