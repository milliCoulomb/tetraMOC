// src/BoltzmannSolver.cpp
#include "BoltzmannSolver.hpp"
#include "Logger.hpp"

#include <cmath>
#include <algorithm>
#include <omp.h> // For OpenMP parallelization

BoltzmannSolver::BoltzmannSolver(const InputHandler& input_handler,
                                 const MeshHandler& mesh_handler,
                                 const std::vector<TrackingData>& tracking_data,
                                 const AngularQuadrature& angular_quadrature,
                                 const SolverParams& params)
    : input_(input_handler),
      mesh_(mesh_handler),
      tracking_data_(tracking_data),
      angular_quadrature_(angular_quadrature),
      params_(params),
      k_eff_(params.initial_k_eff),
      k_eff_old_(params.initial_k_eff) {
    num_groups_ = input_.getNumGroups();
    num_cells_ = static_cast<int>(mesh_.getCells().size());
    QuadratureTotalWeight_ = angular_quadrature_.getTotalWeight();
    Logger::info("BoltzmannSolver initialized with " + std::to_string(num_groups_) + " groups and " + 
                 std::to_string(num_cells_) + " cells.");
}

std::vector<double> BoltzmannSolver::computeScatteringSource(const std::vector<double>& scalar_flux, const int group) const {
    std::vector<double> scat_source(num_cells_, 0.0);

    // Compute sigma_s * phi for each cell
    double sigma_s = input_.getEnergyGroupData(group).scattering_xs;

    #pragma omp parallel for
    for (int cell = 0; cell < num_cells_; ++cell) {
        scat_source[cell] = sigma_s * scalar_flux[cell] / QuadratureTotalWeight_;
    }

    return scat_source;
}

std::vector<double> BoltzmannSolver::solveOneGroupWithSource(const std::vector<double>& external_source, const int group, double eps) {
    if (external_source.size() != static_cast<size_t>(num_cells_)) {
        Logger::error("External source size does not match number of cells.");
        return {};
    }

    // Initialize scalar flux with initial guess (e.g., ones)
    std::vector<double> old_flux(num_cells_, 1.0);
    std::vector<double> new_flux(num_cells_, 0.0);
    double residual = 1.0;
    int iteration = 0;

    Logger::info("Starting one-group solver with source.");

    while (residual > eps && iteration < params_.max_iterations) {
        Logger::info("One-group Iteration " + std::to_string(iteration + 1));

        // Compute scattering source: sigma_s * phi_old
        std::vector<double> scat_source = computeScatteringSource(old_flux, group);

        // Compute total source: scat_source + external_source
        std::vector<double> total_source(num_cells_, 0.0);

        #pragma omp parallel for
        for (int cell = 0; cell < num_cells_; ++cell) {
            total_source[cell] = scat_source[cell] + external_source[cell] / QuadratureTotalWeight_;
        }

        // Compute flux using FluxSolver
        // Assuming FluxSolver has been initialized with TrackingData corresponding to directions and cells
        FluxSolver flux_solver(mesh_, tracking_data_, angular_quadrature_, input_.getEnergyGroupData(group).total_xs);

        // Compute flux based on the total source
        flux_solver.computeFlux(total_source);

        // Collapse flux to scalar flux
        std::vector<double> collapsed_flux = flux_solver.collapseFlux();

        // Update new_flux
        new_flux = collapsed_flux;

        // Compute residual: ||new_flux - old_flux|| / ||old_flux||
        double norm_diff = 0.0;
        double norm_old = 0.0;

        #pragma omp parallel for reduction(+:norm_diff, norm_old)
        for (int cell = 0; cell < num_cells_; ++cell) {
            double diff = new_flux[cell] - old_flux[cell];
            norm_diff += diff * diff;
            norm_old += old_flux[cell] * old_flux[cell];
        }

        residual = (norm_old > 0.0) ? std::sqrt(norm_diff) / std::sqrt(norm_old) : 0.0;
        Logger::info("Residual: " + std::to_string(residual));

        // Prepare for next iteration
        old_flux = new_flux;
        iteration++;
    }

    if (residual <= eps) {
        Logger::info("One-group solver converged in " + std::to_string(iteration) + " iterations.");
    } else {
        Logger::warning("One-group solver did not converge within the maximum iterations.");
    }

    return new_flux;
}

void BoltzmannSolver::updateKEff(const std::vector<double>& fission_source_new,
                                 const std::vector<double>& fission_source_old) {
    double sum_new = 0.0;
    double sum_old = 0.0;

    // Compute total fission source
    #pragma omp parallel for reduction(+:sum_new, sum_old)
    for (int cell = 0; cell < num_cells_; ++cell) {
        sum_new += fission_source_new[cell];
        sum_old += fission_source_old[cell];
    }

    if (sum_old < EPSILON) {
        Logger::error("Old fission source sum is zero, cannot update k_eff.");
        return;
    }

    double new_k_eff = sum_new / sum_old;

    // Thread-safe update of k_eff
    {
        std::lock_guard<std::mutex> lock(k_eff_mutex_);
        k_eff_old_ = k_eff_;
        k_eff_ = new_k_eff;
    }
}

double BoltzmannSolver::getKEff() const {
    return k_eff_;
}