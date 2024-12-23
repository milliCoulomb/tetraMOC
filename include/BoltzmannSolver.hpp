// include/BoltzmannSolver.hpp
#ifndef BOLTZMANN_SOLVER_HPP
#define BOLTZMANN_SOLVER_HPP

#include "InputHandler.hpp"
#include "MeshHandler.hpp"
#include "FluxSolver.hpp" // Include FluxSolver
#include "AngularQuadrature.hpp" // Include AngularQuadrature
#include "TrackingData.hpp" // Include TrackingData
#include "Vector3D.hpp"

#include <vector>
#include <mutex>

class BoltzmannSolver {
public:
    /**
     * @brief Struct to hold solver parameters.
     */
    struct SolverParams {
        double convergence_threshold;
        int max_iterations;
        double initial_k_eff;
    
        SolverParams()
            : convergence_threshold(1e-5),
              max_iterations(1000),
              initial_k_eff(1.0)
        {}
    };

    /**
     * @brief Constructor for BoltzmannSolver.
     * 
     * @param input_handler Reference to InputHandler containing cross-section data.
     * @param mesh_handler Reference to MeshHandler containing mesh data.
     * @param tracking_data Vector of TrackingData containing ray tracing data.
     * @param angular_quadrature Reference to AngularQuadrature containing angles and weights.
     * @param params Solver parameters.
     */
    BoltzmannSolver(const InputHandler& input_handler,
                   const MeshHandler& mesh_handler,
                   const std::vector<TrackingData>& tracking_data,
                   const AngularQuadrature& angular_quadrature,
                   const SolverParams& params = SolverParams());

    /**
     * @brief Solves the one-group transport equation with an external source.
     * 
     * @param external_source External source vector (size = number of cells).
     * @param eps Convergence threshold.
     * @return std::vector<double> Converged scalar flux vector.
     */
    std::vector<double> solveOneGroupWithSource(const std::vector<double>& external_source, const int group, double eps = 1e-5);

    /**
     * @brief Retrieves the computed effective multiplication factor (k_eff).
     * 
     * @return double k_eff value.
     */
    double getKEff() const;

private:
    const InputHandler& input_;
    const MeshHandler& mesh_;
    const std::vector<TrackingData>& tracking_data_;
    const AngularQuadrature& angular_quadrature_;
    SolverParams params_;

    int num_groups_;
    int num_cells_;

    double k_eff_;
    double k_eff_old_;

    // Mutex for thread-safe k_eff update
    std::mutex k_eff_mutex_;

    /**
     * @brief Updates the effective multiplication factor (k_eff).
     * 
     * @param fission_source_new New fission source vector.
     * @param fission_source_old Old fission source vector.
     */
    void updateKEff(const std::vector<double>& fission_source_new,
                   const std::vector<double>& fission_source_old);

    /**
     * @brief Computes the scattering source.
     * 
     * @param scalar_flux Current scalar flux vector.
     * @return std::vector<double> Scattering source vector.
     */
    std::vector<double> computeScatteringSource(const std::vector<double>& scalar_flux, const int group) const;
};

#endif // BOLTZMANN_SOLVER_HPP