// src/FluxSolver.cpp
#include "FluxSolver.hpp"
#include <cmath>

// Constructor
FluxSolver::FluxSolver(const MeshHandler& mesh,
                       const Field& field,
                       const std::vector<TrackingData>& tracking_data,
                       const AngularQuadrature& angular_quadrature,
                       double sigma_t)
    : mesh_(mesh),
      field_(field),
      tracking_data_(tracking_data),
      angular_quadrature_(angular_quadrature),
      sigma_t_(sigma_t) {
    initializeFluxData();
}

// Initialize flux data structure
void FluxSolver::initializeFluxData() {
    size_t num_cells = mesh_.getCells().size();
    size_t num_directions = angular_quadrature_.getDirections().size();

    flux_data_.resize(num_cells, std::vector<CellFlux>(num_directions, CellFlux()));
}

// Compute line-averaged flux for a single CellTrace
std::pair<double, double> FluxSolver::computeLineAveragedFlux(double Q_k, double psi_in, double sigma_t, double L_k) const {
    // Compute psi_out using the provided formula
    double psi_out = psi_in + (Q_k / sigma_t - psi_in) * (1.0 - std::exp(-sigma_t * L_k));

    // Compute line-averaged flux
    double line_avg_flux = (Q_k / sigma_t) - (psi_out - psi_in) / (sigma_t * L_k);

    return {line_avg_flux, psi_out};
}

// Find direction index in angular quadrature
size_t FluxSolver::findDirectionIndex(const Vector3D& direction) const {
    const std::vector<Direction>& directions = angular_quadrature_.getDirections();
    for(size_t i = 0; i < directions.size(); ++i) {
        // Reconstruct the 3D direction vector from mu and phi
        double sqrt_term = std::sqrt(1.0 - directions[i].mu * directions[i].mu);
        double x = sqrt_term * std::cos(directions[i].phi);
        double y = sqrt_term * std::sin(directions[i].phi);
        double z = directions[i].mu;
        Vector3D dir_vec(x, y, z);

        if (dir_vec.isAlmostEqual(direction, direction_tolerance_)) {
            return i;
        }
    }
    // If not found, return an invalid index
    return directions.size();
}

// Method to compute flux
void FluxSolver::computeFlux() {
    int num_threads = omp_get_max_threads();

    // Create a local flux_data buffer for each thread
    std::vector<std::vector<std::vector<CellFlux>>> local_flux_data(num_threads,
        std::vector<std::vector<CellFlux>>(mesh_.getCells().size(),
            std::vector<CellFlux>(angular_quadrature_.getDirections().size(), CellFlux())));

    // Parallelize over tracking_data
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        std::vector<std::vector<CellFlux>>& thread_flux = local_flux_data[thread_id];

        #pragma omp for schedule(dynamic)
        for(size_t ray_idx = 0; ray_idx < tracking_data_.size(); ++ray_idx) {
            const TrackingData& ray = tracking_data_[ray_idx];
            const Vector3D& direction = ray.direction;

            // Find direction index
            size_t dir_index = findDirectionIndex(direction);
            if(dir_index >= angular_quadrature_.getDirections().size()) {
                Logger::warning("Direction not found in angular quadrature for ray ID " + std::to_string(ray.ray_id) + ". Skipping this ray.");
                continue;
            }

            // Initialize psi_in for the ray
            double psi_in = 0.0;

            // Iterate over cell traces in order
            for(const auto& trace : ray.cell_traces) {
                if(trace.cell_id < 0 || trace.cell_id >= static_cast<int>(mesh_.getCells().size())) {
                    Logger::warning("Invalid cell ID " + std::to_string(trace.cell_id) + " in ray ID " + std::to_string(ray.ray_id) + ". Skipping this trace.");
                    continue;
                }

                // Retrieve Q_k for the cell
                if(trace.cell_id >= static_cast<int>(field_.getScalarFields().size())) {
                    // If scalar field not set for this cell, assume Q_k = 0
                    Logger::warning("Scalar field Q_k not set for cell ID " + std::to_string(trace.cell_id) + " in ray ID " + std::to_string(ray.ray_id) + ". Assuming Q_k = 0.");
                }
                double Q_k = (trace.cell_id < static_cast<int>(field_.getScalarFields().size())) ?
                             field_.getScalarFields()[trace.cell_id].value : 0.0;

                // Retrieve L_k (time spent in the cell)
                double L_k = trace.time_spent;

                // Compute line-averaged flux and psi_out
                auto [line_avg_flux, psi_out] = computeLineAveragedFlux(Q_k, psi_in, sigma_t_, L_k);

                // Accumulate flux for the direction
                thread_flux[trace.cell_id][dir_index].flux += line_avg_flux * L_k;
                thread_flux[trace.cell_id][dir_index].weight += L_k;

                // Update psi_in for the next cell
                psi_in = psi_out;
            }
        }
    }

    // Combine local_flux_data into global flux_data_ with parallelization
    #pragma omp parallel for schedule(static) // Parallelize over cells
    for (size_t cell = 0; cell < flux_data_.size(); ++cell) {
        for (size_t dir = 0; dir < flux_data_[cell].size(); ++dir) {
            double total_flux = 0.0;
            double total_weight = 0.0;
            for(int thread_id = 0; thread_id < num_threads; ++thread_id) {
                total_flux += local_flux_data[thread_id][cell][dir].flux;
                total_weight += local_flux_data[thread_id][cell][dir].weight;
            }
            #pragma omp atomic
            flux_data_[cell][dir].flux += total_flux;
            #pragma omp atomic
            flux_data_[cell][dir].weight += total_weight;
        }
    }


    // Normalize flux by weights
    #pragma omp parallel for schedule(dynamic)
    for(size_t cell = 0; cell < flux_data_.size(); ++cell) {
        for(size_t dir = 0; dir < flux_data_[cell].size(); ++dir) {
            if(flux_data_[cell][dir].weight > 0.0) {
                flux_data_[cell][dir].flux /= flux_data_[cell][dir].weight;
            }
        }
    }
}