// include/FluxSolver.hpp
#ifndef FLUX_SOLVER_HPP
#define FLUX_SOLVER_HPP

#include "MeshHandler.hpp"
#include "Field.hpp"
#include "TrackingData.hpp"
#include "AngularQuadrature.hpp"
#include "Vector3D.hpp"
#include <vector>
#include <utility> // For std::pair
#include "Logger.hpp"
#include <omp.h>

// Structure to store flux data for a single cell and direction
struct CellFlux {
    double flux = 0.0;   // Accumulated flux
    double weight = 0.0; // Number of contributions
};

class FluxSolver {
public:
    // Constructor
    FluxSolver(const MeshHandler& mesh,
               const Field& field,
               const std::vector<TrackingData>& tracking_data,
               const AngularQuadrature& angular_quadrature,
               double sigma_t);

    // Method to compute flux
    void computeFlux();

    // Getter for flux data
    const std::vector<std::vector<CellFlux>>& getFluxData() const { return flux_data_; }

    // Method to collapse the flux in all directions to a scalar flux by using the angular quadrature weights
    std::vector<double> collapseFlux() const;

private:
    const MeshHandler& mesh_;
    const Field& field_;
    const std::vector<TrackingData>& tracking_data_;
    const AngularQuadrature& angular_quadrature_;
    double sigma_t_; // Total cross section

    // Data structure to store flux per cell per direction
    // flux_data_[cell_id][direction_id]
    std::vector<std::vector<CellFlux>> flux_data_;

    // Helper method to initialize flux data
    void initializeFluxData();

    // Helper method to compute line-averaged flux for a single CellTrace
    // Returns a pair: {line_avg_flux, psi_out}
    std::pair<double, double> computeLineAveragedFlux(double Q_k, double psi_in, double sigma_t, double L_k) const;

    // Helper method to find direction index in angular quadrature
    // Returns the index if found, else returns angular_quadrature_.getDirections().size()
    size_t findDirectionIndex(const Vector3D& direction) const;

    // Tolerance for direction matching
    const double direction_tolerance_ = 1e-6;
};

#endif // FLUX_SOLVER_HPP
