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

        // Initialize MeshHandler and load mesh data
        MeshHandler mesh_handler;
        mesh_handler.loadNodes(input_deck.mesh.nodes);
        mesh_handler.loadCells(input_deck.mesh.cells);
        mesh_handler.loadFaceConnectivity(input_deck.mesh.faces);

        // Initialize dummy Field
        Field field;

        bool constant_dir_bool = true;
        bool use_half_quadrature_for_constant = false;

        // Initialize RayTracerManager and generate tracking data
        RayTracerManager ray_tracer_manager(mesh_handler, field, angular_quadrature, constant_dir_bool,
                                            use_half_quadrature_for_constant);
        
        int rays_per_face = input_deck.solver_parameters.rays_per_face;
        ray_tracer_manager.generateTrackingData(rays_per_face);

        InputHandler input_handler;
        input_handler.loadData(input_deck.cross_sections.data_files[0]);

        Settings settings;
        settings.setMultiGroupTolerance(input_deck.solver_parameters.multi_group_tolerance);
        settings.setMultiGroupMaxIterations(input_deck.solver_parameters.multi_group_max_iterations);
        settings.setOneGroupTolerance(input_deck.solver_parameters.one_group_tolerance);
        settings.setOneGroupMaxIterations(input_deck.solver_parameters.one_group_max_iterations);
        settings.setFissionSourceTolerance(input_deck.solver_parameters.fission_source_tolerance);
        settings.setKeffTolerance(input_deck.solver_parameters.keff_tolerance);

        // Initialize BoltzmannSolver
        BoltzmannSolver solver(input_handler, mesh_handler, ray_tracer_manager.getTrackingData(),
                               angular_quadrature, settings);

        // Solve eigenvalue problem
        bool converged = solver.solveEigenvalueProblem();

        if (!converged) {
            std::cerr << "Eigenvalue problem did not converge." << std::endl;
            return 1;
        }
        // Output k_eff
        std::cout << "Converged k_eff: " << solver.getKEff() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}