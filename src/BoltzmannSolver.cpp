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
                                 const Settings& settings)
    : input_(input_handler),
      mesh_(mesh_handler),
      tracking_data_(tracking_data),
      angular_quadrature_(angular_quadrature),
      settings_(settings),
      k_eff_(1.0),
      k_eff_old_(1.0) {
    num_groups_ = input_.getNumGroups();
    num_cells_ = static_cast<int>(mesh_.getCells().size());
    QuadratureTotalWeight_ = angular_quadrature_.getTotalWeight();
    Logger::info("BoltzmannSolver initialized with " + std::to_string(num_groups_) + " groups and " + 
                 std::to_string(num_cells_) + " cells.");
}

std::vector<double> BoltzmannSolver::computeScatteringSource(const std::vector<double>& scalar_flux, const int group) const {
    std::vector<double> scat_source(num_cells_, 0.0);

    // Compute sigma_s * phi for each cell
    double sigma_s = input_.getSelfScatteringXS(group);

    #pragma omp parallel for
    for (int cell = 0; cell < num_cells_; ++cell) {
        scat_source[cell] = sigma_s * scalar_flux[cell] / QuadratureTotalWeight_;
    }

    return scat_source;
}

std::vector<std::vector<double>> BoltzmannSolver::computeMultiGroupScatteringSource(const std::vector<std::vector<double>>& scalar_flux) const {
    std::vector<std::vector<double>> scat_source(num_groups_, std::vector<double>(num_cells_, 0.0));

    // Precompute scattering cross sections to reduce repeated access
    std::vector<std::vector<double>> scattering_xs(num_groups_, std::vector<double>(num_groups_, 0.0));
    for(int g_prime = 0; g_prime < num_groups_; ++g_prime){
        for(int group = 0; group < num_groups_; ++group){
            scattering_xs[g_prime][group] = input_.getEnergyGroupData(g_prime).scattering_xs[group];
        }
    }

    #pragma omp parallel for
    for(int cell = 0; cell < num_cells_; ++cell) {
        for(int group = 0; group < num_groups_; ++group) {
            double sum = 0.0;
            for(int g_prime = 0; g_prime < num_groups_; ++g_prime) {
                if(g_prime != group) {
                    sum += scattering_xs[g_prime][group] * scalar_flux[g_prime][cell];
                }
            }
            scat_source[group][cell] = sum;
        }
    }
    return scat_source;
}

std::vector<std::vector<double>> BoltzmannSolver::computeFissionSource(const std::vector<std::vector<double>>& scalar_flux,
    const double old_keff) const
{
    std::vector<std::vector<double>> fission_source(num_groups_, std::vector<double>(num_cells_, 0.0));

    // Precompute fission cross sections to reduce repeated access
    std::vector<double> fission_xs(num_groups_, 0.0);
    for(int group = 0; group < num_groups_; ++group){
        fission_xs[group] = input_.getEnergyGroupData(group).fission_xs * input_.getEnergyGroupData(group).multiplicity / old_keff;
    }

    #pragma omp parallel for collapse(2)
    for(int cell = 0; cell < num_cells_; ++cell) {
        for(int group = 0; group < num_groups_; ++group) {
            fission_source[group][cell] = fission_xs[group] * scalar_flux[group][cell] * input_.getEnergyGroupData(group).fission_spectrum;
        }
    }
    return fission_source;
}

std::vector<double> BoltzmannSolver::solveOneGroupWithSource(
    const std::vector<double>& external_source, 
    const int group, 
    const std::vector<double>& initial_guess) 
{
    if (external_source.size() != static_cast<size_t>(num_cells_)) {
        Logger::error("External source size does not match number of cells.");
        return {};
    }

    // Initialize scalar flux with initial guess or ones
    std::vector<double> old_flux = initial_guess;
    if (old_flux.empty() || 
        static_cast<int>(old_flux.size()) != num_cells_) 
    {
        old_flux.assign(num_cells_, 1.0);
    }

    // Initialize new_flux once and reuse
    std::vector<double> new_flux(num_cells_, 0.0);
    double residual = 1.0;
    int iteration = 0;

    Logger::info("Starting one-group solver with source.");

    while (residual > settings_.getOneGroupTolerance() && iteration < settings_.getOneGroupMaxIterations()) {
        Logger::info("One-group Iteration " + std::to_string(iteration + 1));

        // Compute scattering source: sigma_s * phi_old
        std::vector<double> scat_source = computeScatteringSource(old_flux, group);

        // Compute total source: scat_source + external_source / QuadratureTotalWeight_
        #pragma omp parallel for
        for (int cell = 0; cell < num_cells_; ++cell) {
            scat_source[cell] += external_source[cell] / QuadratureTotalWeight_;
        }

        // Compute flux using FluxSolver
        FluxSolver flux_solver(mesh_, tracking_data_, angular_quadrature_, input_.getEnergyGroupData(group).total_xs);

        // Compute flux based on the total source
        flux_solver.computeFlux(scat_source);

        // Collapse flux to scalar flux
        std::vector<double> collapsed_flux = flux_solver.collapseFlux();

        // Update new_flux
        new_flux = std::move(collapsed_flux);

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

        // Swap old_flux and new_flux to reuse memory
        std::swap(old_flux, new_flux);
        // Reset new_flux for next iteration
        #pragma omp parallel for
        for (int cell = 0; cell < num_cells_; ++cell) {
            new_flux[cell] = 0.0;
        }

        iteration++;
    }

    if (residual <= settings_.getOneGroupTolerance()) {
        Logger::info("One-group solver converged in " + std::to_string(iteration) + " iterations.");
    } else {
        Logger::warning("One-group solver did not converge within the maximum iterations.");
    }

    return old_flux;
}

std::vector<std::vector<double>> BoltzmannSolver::solveMultiGroupWithSource(
    const std::vector<std::vector<double>>& external_source,
    const std::vector<std::vector<double>>& initial_guess) 
{
    if (external_source.size() != static_cast<size_t>(num_groups_) || 
        external_source[0].size() != static_cast<size_t>(num_cells_)) 
    {
        Logger::error("External source size does not match number of groups or cells.");
        return {};
    }

    // Initialize scalar flux with initial guess or ones
    std::vector<std::vector<double>> old_flux = initial_guess;
    if (old_flux.empty() || 
        static_cast<int>(old_flux.size()) != num_groups_ || 
        static_cast<int>(old_flux[0].size()) != num_cells_) 
    {
        old_flux.assign(num_groups_, std::vector<double>(num_cells_, 1.0));
    }

    // Initialize new_flux once and reuse
    std::vector<std::vector<double>> new_flux(num_groups_, std::vector<double>(num_cells_, 0.0));
    double residual = 1.0;
    int iteration = 0;

    Logger::info("Starting multi-group solver with source.");

    while (residual > settings_.getMultiGroupTolerance() && iteration < settings_.getMultiGroupMaxIterations()) {
        Logger::info("Multi-group Iteration " + std::to_string(iteration + 1));

        // Compute scattering source: sigma_s * phi_old
        std::vector<std::vector<double>> scat_source = computeMultiGroupScatteringSource(old_flux);

        // Compute total source: scat_source + external_source
        #pragma omp parallel for collapse(2)
        for (int group = 0; group < num_groups_; ++group) {
            for (int cell = 0; cell < num_cells_; ++cell) {
                scat_source[group][cell] += external_source[group][cell];
            }
        }

        // Solve one group problem for each group
        #pragma omp parallel for
        for(int group = 0; group < num_groups_; ++group) {
            new_flux[group] = solveOneGroupWithSource(scat_source[group], group, old_flux[group]);
        }

        // Compute residual: ||new_flux - old_flux|| / ||old_flux||
        double norm_diff = 0.0;
        double norm_old = 0.0;

        #pragma omp parallel for reduction(+:norm_diff, norm_old) collapse(2) schedule(dynamic)
        for (int group = 0; group < num_groups_; ++group) {
            for (int cell = 0; cell < num_cells_; ++cell) {
                double diff = new_flux[group][cell] - old_flux[group][cell];
                norm_diff += diff * diff;
                norm_old += old_flux[group][cell] * old_flux[group][cell];
            }
        }

        residual = (norm_old > 0.0) ? std::sqrt(norm_diff) / std::sqrt(norm_old) : 0.0;
        Logger::info("Residual: " + std::to_string(residual));

        // Swap old_flux and new_flux to reuse memory
        std::swap(old_flux, new_flux);
        // Reset new_flux for next iteration
        #pragma omp parallel for collapse(2)
        for (int group = 0; group < num_groups_; ++group) {
            for (int cell = 0; cell < num_cells_; ++cell) {
                new_flux[group][cell] = 0.0;
            }
        }

        iteration++;
    }

    if (residual <= settings_.getMultiGroupTolerance()) {
        Logger::info("Multi-group solver converged in " + std::to_string(iteration) + " iterations.");
    }
    if (iteration >= settings_.getMultiGroupMaxIterations()) {
        Logger::warning("Multi-group solver did not converge within the maximum iterations.");
    }
    return old_flux;
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