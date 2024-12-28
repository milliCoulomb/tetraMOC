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
#include <fstream>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_deck.yaml>" << std::endl;
        return 1;
    }

    try {
        // Parse InputDeck
        InputDeck input_deck = InputDeckParser::parse(argv[1]);

        // Initialize AngularQuadrature
        AngularQuadrature angular_quadrature(input_deck.angular_quadrature_parameters.ntheta,
                                             input_deck.angular_quadrature_parameters.nphi);
        std::cout << "Number of angles: " << angular_quadrature.getDirections().size() << std::endl;
        std::cout << "Total weight: " << angular_quadrature.getTotalWeight() << std::endl;

        // Initialize MeshHandler and load mesh data
        MeshHandler mesh_handler;

        mesh_handler.loadNodes(input_deck.mesh.nodes);
        mesh_handler.loadCells(input_deck.mesh.cells);
        mesh_handler.loadFaceConnectivity(input_deck.mesh.faces);

        auto mesh_bounds = mesh_handler.getBoundaryFaces();
        std::cout << "Number of boundary faces: " << mesh_bounds.size() << std::endl;
        // display number of nodes
        const int number_of_nodes = static_cast<int>(mesh_handler.getNodes().size());
        // display number of faces
        const int number_of_faces = static_cast<int>(mesh_handler.getFaces().size());
        // display number of cells
        const int number_of_cells = static_cast<int>(mesh_handler.getCells().size());
        std::cout << "Number of nodes: " << number_of_nodes << std::endl;
        std::cout << "Number of faces: " << number_of_faces << std::endl;
        std::cout << "Number of cells: " << number_of_cells << std::endl;

        // Initialize dummy Field
        // Field field;

        // bool constant_dir_bool = true;
        // bool use_half_quadrature_for_constant = false;

        // // Initialize RayTracerManager and generate tracking data
        // RayTracerManager ray_tracer_manager(mesh_handler, field, angular_quadrature, constant_dir_bool,
        //                                     use_half_quadrature_for_constant);
        
        // int rays_per_face = input_deck.solver_parameters.rays_per_face;
        // ray_tracer_manager.generateTrackingData(rays_per_face);

        InputHandler input_handler;
        input_handler.loadData(input_deck.cross_sections.data_files[0]);
        std::cout << "Number of groups: " << input_handler.getNumGroups() << std::endl;

        // Settings settings;
        // settings.setMultiGroupTolerance(input_deck.solver_parameters.multi_group_tolerance);
        // settings.setMultiGroupMaxIterations(input_deck.solver_parameters.multi_group_max_iterations);
        // settings.setOneGroupTolerance(input_deck.solver_parameters.one_group_tolerance);
        // settings.setOneGroupMaxIterations(input_deck.solver_parameters.one_group_max_iterations);
        // settings.setFissionSourceTolerance(input_deck.solver_parameters.fission_source_tolerance);
        // settings.setKeffTolerance(input_deck.solver_parameters.keff_tolerance);

        // // get the output path
        // std::string output_path_flux = input_deck.output.flux_output_file;
        // std::string output_path_keff = input_deck.output.k_eff_output_file;

        // // Initialize BoltzmannSolver
        // BoltzmannSolver solver(input_handler, mesh_handler, ray_tracer_manager.getTrackingData(),
        //                        angular_quadrature, settings);

        // // Solve eigenvalue problem
        // bool converged = solver.solveEigenvalueProblem();

        // if (!converged) {
        //     std::cerr << "Eigenvalue problem did not converge." << std::endl;
        //     return 1;
        // }
        // // Output k_eff
        // std::cout << "Converged k_eff: " << solver.getKEff() << std::endl;
        // // store the flux
        // std::vector<std::vector<double>> flux = solver.getScalarFlux();
        
        // // OutputHandler instance
        // OutputHandler output_handler;
        // output_handler.writeScalarFlux(output_path_flux, flux);
        // output_handler.writeKEff(output_path_keff, solver.getKEff());
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}