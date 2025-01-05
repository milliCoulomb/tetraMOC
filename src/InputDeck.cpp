// src/InputDeck.cpp

#include "InputDeck.hpp"
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <exception>

InputDeck InputDeckParser::parse(const std::string& filename) {
    InputDeck deck;

    try {
        YAML::Node config = YAML::LoadFile(filename);

        // Parse mesh paths
        YAML::Node mesh_node = config["mesh"];
        if (!mesh_node) throw std::runtime_error("Missing 'mesh' section in YAML.");
        deck.mesh.nodes = mesh_node["nodes"].as<std::string>();
        deck.mesh.cells = mesh_node["cells"].as<std::string>();
        deck.mesh.faces = mesh_node["faces"].as<std::string>();

        // Parse cross sections
        YAML::Node cs_node = config["cross_sections"];
        if (!cs_node) throw std::runtime_error("Missing 'cross_sections' section in YAML.");
        for (const auto& file : cs_node["data_files"]) {
            deck.cross_sections.data_files.push_back(file.as<std::string>());
        }

        // Parse angular quadrature parameters
        YAML::Node aq_node = config["angular_quadrature"];
        if (!aq_node) throw std::runtime_error("Missing 'angular_quadrature' section in YAML.");
        deck.angular_quadrature_parameters.ntheta = aq_node["ntheta"].as<int>();
        deck.angular_quadrature_parameters.nphi = aq_node["nphi"].as<int>();

        // Parse solver parameters
        YAML::Node sp_node = config["solver_parameters"];
        if (!sp_node) throw std::runtime_error("Missing 'solver_parameters' section in YAML.");
        deck.solver_parameters.multi_group_max_iterations = sp_node["multi_group_max_iterations"].as<int>();
        deck.solver_parameters.multi_group_tolerance = sp_node["multi_group_tolerance"].as<double>();
        deck.solver_parameters.one_group_max_iterations = sp_node["one_group_max_iterations"].as<int>();
        deck.solver_parameters.one_group_tolerance = sp_node["one_group_tolerance"].as<double>();
        deck.solver_parameters.fission_source_tolerance = sp_node["fission_source_tolerance"].as<double>();
        deck.solver_parameters.keff_tolerance = sp_node["keff_tolerance"].as<double>();
        deck.solver_parameters.max_power_iterations = sp_node["max_power_iterations"].as<int>();
        deck.solver_parameters.rays_per_face = sp_node["rays_per_face"].as<int>();
        deck.solver_parameters.max_ray_length = sp_node["max_ray_length"].as<int>();
        deck.solver_parameters.use_half_hemisphere = sp_node["use_half_hemisphere"].as<bool>();

        // Parse output settings
        YAML::Node out_node = config["output"];
        if (!out_node) throw std::runtime_error("Missing 'output' section in YAML.");
        if (!out_node["flux_output_file"] || !out_node["flux_output_file"].IsScalar()) {
            throw YAML::BadConversion(out_node["flux_output_file"].Mark());
        }
        deck.output.flux_output_file = out_node["flux_output_file"].as<std::string>();
        deck.output.k_eff_output_file = out_node["k_eff_output_file"].as<std::string>();

        // Parse logging settings
        YAML::Node log_node = config["logging"];
        if (!log_node) throw std::runtime_error("Missing 'logging' section in YAML.");
        deck.logging.level = log_node["level"].as<std::string>();

    } catch (const YAML::Exception& e) {
        std::cerr << "YAML Parsing Error: " << e.what() << std::endl;
        throw;
    } catch (const std::exception& e) {
        std::cerr << "InputDeck Parsing Error: " << e.what() << std::endl;
        throw;
    }

    return deck;
}
