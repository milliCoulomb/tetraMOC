// src/main.cpp
#include <iostream>
#include "InputDeck.hpp"
#include "AngularQuadrature.hpp"
#include "MeshHandler.hpp"
#include "Field.hpp"
#include "RayTracerManager.hpp"
#include "BoltzmannSolver.hpp"
#include "InputHandler.hpp"
#include "Settings.hpp"
#include "OutputHandler.hpp"
#include "Logger.hpp"
#include <fstream>


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_deck.yaml>" << std::endl;
        return 1;
    }

    try {
        // Parse InputDeck
        InputDeck input_deck = InputDeckParser::parse(argv[1]);

        // Set logger level
        if (input_deck.logging.level == "INFO") {
            Logger::setLogLevel(LogLevel::INFO);
        } else if (input_deck.logging.level == "WARNING") {
            Logger::setLogLevel(LogLevel::WARNING);
        } else if (input_deck.logging.level == "ERROR") {
            Logger::setLogLevel(LogLevel::ERROR);
        } else {
            Logger::setLogLevel(LogLevel::RUNNING); // default
        }

        // Initialize AngularQuadrature
        AngularQuadrature angular_quadrature(input_deck.angular_quadrature_parameters.ntheta,
                                             input_deck.angular_quadrature_parameters.nphi);
        
        Logger::info("Number of angles: " + std::to_string(angular_quadrature.getDirections().size()));
        Logger::info("Total weight: " + std::to_string(angular_quadrature.getTotalWeight()));

        Logger::running("Angular quadrature generated.");

        // Initialize MeshHandler and load mesh data
        MeshHandler mesh_handler;

        mesh_handler.loadNodes(input_deck.mesh.nodes);
        mesh_handler.loadCells(input_deck.mesh.cells);
        mesh_handler.loadFaceConnectivity(input_deck.mesh.faces);

        auto mesh_bounds = mesh_handler.getBoundaryFaces();
        // display number of nodes
        const int number_of_nodes = static_cast<int>(mesh_handler.getNodes().size());
        // display number of faces
        const int number_of_faces = static_cast<int>(mesh_handler.getFaces().size());
        // display number of cells
        const int number_of_cells = static_cast<int>(mesh_handler.getCells().size());

        Logger::info("Number of nodes: " + std::to_string(number_of_nodes));
        Logger::info("Number of faces: " + std::to_string(number_of_faces));
        Logger::info("Number of cells: " + std::to_string(number_of_cells));
        Logger::running("Mesh data loaded.");

        // Initialize dummy Field
        Field field;

        bool constant_dir_bool = true;
        bool use_half_quadrature_for_constant = false;

        // Initialize RayTracerManager and generate tracking data
        RayTracerManager ray_tracer_manager(mesh_handler, field, angular_quadrature, constant_dir_bool,
                                            use_half_quadrature_for_constant);
        
        int rays_per_face = input_deck.solver_parameters.rays_per_face;
        ray_tracer_manager.generateTrackingData(rays_per_face);

        Logger::running("Tracking data generated.");

        InputHandler input_handler;
        input_handler.loadData(input_deck.cross_sections.data_files[0]);

        Settings settings;
        settings.setMultiGroupTolerance(input_deck.solver_parameters.multi_group_tolerance);
        settings.setMultiGroupMaxIterations(input_deck.solver_parameters.multi_group_max_iterations);
        settings.setOneGroupTolerance(input_deck.solver_parameters.one_group_tolerance);
        settings.setOneGroupMaxIterations(input_deck.solver_parameters.one_group_max_iterations);
        settings.setFissionSourceTolerance(input_deck.solver_parameters.fission_source_tolerance);
        settings.setKeffTolerance(input_deck.solver_parameters.keff_tolerance);
        settings.setMaxPowerIterations(input_deck.solver_parameters.max_power_iterations);

        // // get the output path
        std::string output_path_flux = input_deck.output.flux_output_file;
        std::string output_path_keff = input_deck.output.k_eff_output_file;

        // Initialize BoltzmannSolver
        BoltzmannSolver solver(input_handler, mesh_handler, ray_tracer_manager.getTrackingData(),
                               angular_quadrature, settings);

        Logger::running("Solving eigenvalue problem...");
        // Solve eigenvalue problem
        bool converged = solver.solveEigenvalueProblem();

        if (!converged) {
            std::cerr << "Eigenvalue problem did not converge." << std::endl;
            return 1;
        }
        // store the flux
        std::vector<std::vector<double>> flux = solver.getScalarFlux();
        
        // OutputHandler instance
        OutputHandler output_handler;
        output_handler.writeScalarFlux(output_path_flux, flux);
        output_handler.writeKEff(output_path_keff, solver.getKEff());
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}