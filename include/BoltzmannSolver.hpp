// include/BoltzmannSolver.hpp
#ifndef BOLTZMANN_SOLVER_HPP
#define BOLTZMANN_SOLVER_HPP

#include "InputHandler.hpp"
#include "MeshHandler.hpp"
#include "FluxSolver.hpp" // Include FluxSolver
#include "AngularQuadrature.hpp" // Include AngularQuadrature
#include "TrackingData.hpp" // Include TrackingData
#include "Vector3D.hpp"
#include "Settings.hpp"

#include <vector>
#include <mutex>

class BoltzmannSolver {
public:
    BoltzmannSolver(const InputHandler& input_handler,
                   const MeshHandler& mesh_handler,
                   const std::vector<TrackingData>& tracking_data,
                   const AngularQuadrature& angular_quadrature,
                   const Settings& settings);

    /**
     * @brief Solves the one-group transport equation with an external source.
     * 
     * @param external_source External source vector (size = number of cells).
     * 
     * @return std::vector<double> Converged scalar flux vector.
     */
    std::vector<double> solveOneGroupWithSource(const std::vector<double>& external_source, const int group, const std::vector<double>& initial_guess = {});
    /**
     * @brief Solves the multi-group transport equation with an external source.
     * 
     * @param external_source External source vector std::vector<std::vector<double>> (size = number of groups, each group size = number of cells).
     * 
     * @return std::vector<std::vector<double>> Converged scalar flux vector.
     */
    std::vector<std::vector<double>> solveMultiGroupWithSource(const std::vector<double>& external_source, const std::vector<std::vector<double>>& initial_guess = {});
    /**
     * @brief Solves the eigenvalue problem using the power iteration method.
     * 
     * 
     * @return bool True if the eigenvalue problem converged, false otherwise.
     */
    // for now return a NotImplemented error
    bool solveEigenvalueProblem(const std::vector<std::vector<double>>& initial_guess = {});

    /**
     * @brief Retrieves the computed effective multiplication factor (k_eff).
     * 
     * @return double k_eff value.
     */
    double getKEff() const;

    /**
     * @brief Retrieves the total weight of the angular quadrature.
     * 
     * @return double Total weight.
     */
    double getQuadratureTotalWeight() const { return QuadratureTotalWeight_; }

    /**
     * @brief Retrieves the scalar flux.
     * 
     * @return std::vector<std::vector<double>> Scalar flux vector.
     */
    std::vector<std::vector<double>> getScalarFlux() const { return scalar_flux_; }

private:
    const InputHandler& input_;
    const MeshHandler& mesh_;
    const std::vector<TrackingData>& tracking_data_;
    const AngularQuadrature& angular_quadrature_;
    const Settings& settings_;

    int num_groups_;
    int num_cells_;
    double QuadratureTotalWeight_;

    double k_eff_;
    double k_eff_old_;

    std::vector<std::vector<double>> scalar_flux_;

    // Mutex for thread-safe k_eff update
    std::mutex k_eff_mutex_;

    /**
     * @brief Updates the effective multiplication factor (k_eff).
     * 
     * @param fission_source_new New fission source vector.
     * @param fission_source_old Old fission source vector.
     */
    void updateKEff(const std::vector<double>& fission_source_new, const std::vector<double>& fission_source_old);

    /**
     * @brief Computes the scattering source.
     * 
     * @param scalar_flux Current scalar flux vector.
     * @return std::vector<double> Scattering source vector.
     */
    std::vector<double> computeScatteringSource(const std::vector<double>& scalar_flux, const int group) const;

    /**
     * @brief Computes the multigroup scattering source.
     * 
     * @param scalar_flux Current scalar flux vector (std::vector<std::vector<double>> size = number of groups, each group size = number of cells).
     * @return std::vector<std::vector<double>> Scattering source vector.
     */
    std::vector<std::vector<double>> computeMultiGroupScatteringSource(const std::vector<std::vector<double>>& scalar_flux) const;

    /**
     * @brief Computes the fission source.
     * 
     * @param scalar_flux Current scalar flux vector.
     * @param old_keff Old effective multiplication factor.
     * @return std::vector<double> Fission source vector.
     */

    std::vector<std::vector<double>> computeFissionSource(const std::vector<std::vector<double>>& scalar_flux, const double old_keff) const;

    /**
     * @brief Computes the nu-fission source (not weighted by fission spectrum).
     * 
     * @param scalar_flux Current scalar flux vector.
     * @param old_keff Old effective multiplication factor.
     * @return std::vector<std::vector<double>> Nu-fission source vector.
     */

    std::vector<double> computeNuFissionSource(const std::vector<std::vector<double>>& scalar_flux, const double old_keff) const;
    
};

#endif // BOLTZMANN_SOLVER_HPP