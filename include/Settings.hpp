// include/Settings.hpp
#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <stdexcept>

/**
 * @brief A class to store configuration settings for the BoltzmannSolver and AngularQuadrature.
 */
class Settings {
public:
    /**
     * @brief Default constructor initializing default settings.
     */
    Settings()
        : multi_group_tolerance_(1e-7),
          multi_group_max_iterations_(1000),
          one_group_tolerance_(1e-7),
          one_group_max_iterations_(500),
          fission_source_tolerance_(1e-7),
          keff_tolerance_(1e-6),
          n_theta_(10),
          n_phi_(20) {}

    /**
     * @brief Parameterized constructor to initialize all settings.
     * 
     * @param multi_group_tolerance       Tolerance for the multi-group solver.
     * @param multi_group_max_iterations Max iterations for the multi-group solver.
     * @param one_group_tolerance         Tolerance for the one-group solver.
     * @param one_group_max_iterations   Max iterations for the one-group solver.
     * @param n_theta                     Number of theta points in AngularQuadrature.
     * @param n_phi                       Number of phi points in AngularQuadrature.
     */
    Settings(double multi_group_tolerance, int multi_group_max_iterations,
             double one_group_tolerance, int one_group_max_iterations,
             double fission_source_tolerance, double keff_tolerance,
             int n_theta, int n_phi)
        : multi_group_tolerance_(multi_group_tolerance),
          multi_group_max_iterations_(multi_group_max_iterations),
          one_group_tolerance_(one_group_tolerance),
          one_group_max_iterations_(one_group_max_iterations),
          fission_source_tolerance_(fission_source_tolerance),
          keff_tolerance_(keff_tolerance),
          n_theta_(n_theta),
          n_phi_(n_phi) 
    {
        validate();
    }

    // Getters
    double getMultiGroupTolerance() const { return multi_group_tolerance_; }
    int getMultiGroupMaxIterations() const { return multi_group_max_iterations_; }
    double getOneGroupTolerance() const { return one_group_tolerance_; }
    int getOneGroupMaxIterations() const { return one_group_max_iterations_; }
    double getFissionSourceTolerance() const { return fission_source_tolerance_; }
    double getKeffTolerance() const { return keff_tolerance_; }
    int getNTheta() const { return n_theta_; }
    int getNPhi() const { return n_phi_; }

    // Setters
    void setMultiGroupTolerance(double tolerance) { 
        if(tolerance <= 0.0) {
            throw std::invalid_argument("Multi-group tolerance must be positive.");
        }
        multi_group_tolerance_ = tolerance; 
    }

    void setMultiGroupMaxIterations(int max_iterations) { 
        if(max_iterations <= 0) {
            throw std::invalid_argument("Multi-group max iterations must be positive.");
        }
        multi_group_max_iterations_ = max_iterations; 
    }

    void setOneGroupTolerance(double tolerance) { 
        if(tolerance <= 0.0) {
            throw std::invalid_argument("One-group tolerance must be positive.");
        }
        one_group_tolerance_ = tolerance; 
    }

    void setOneGroupMaxIterations(int max_iterations) { 
        if(max_iterations <= 0) {
            throw std::invalid_argument("One-group max iterations must be positive.");
        }
        one_group_max_iterations_ = max_iterations; 
    }

    void setFissionSourceTolerance(double tolerance) { 
        if(tolerance <= 0.0) {
            throw std::invalid_argument("Fission source tolerance must be positive.");
        }
        fission_source_tolerance_ = tolerance; 
    }

    void setKeffTolerance(double tolerance) { 
        if(tolerance <= 0.0) {
            throw std::invalid_argument("k-effective tolerance must be positive.");
        }
        keff_tolerance_ = tolerance; 
    }

    void setNTheta(int n_theta) { 
        if(n_theta <= 0) {
            throw std::invalid_argument("n_theta must be positive.");
        }
        if(n_theta % 2 != 0) {
            throw std::invalid_argument("n_theta must be a multiple of 2.");
        }
        n_theta_ = n_theta; 
    }

    void setNPhi(int n_phi) { 
        if(n_phi <= 0) {
            throw std::invalid_argument("n_phi must be positive.");
        }
        if (n_phi % 4 != 0) {
            throw std::invalid_argument("n_phi must be a multiple of 4.");
        }
        n_phi_ = n_phi;
    }

private:
    // Member variables storing settings
    double multi_group_tolerance_;          ///< Tolerance for the multi-group solver
    int multi_group_max_iterations_;       ///< Maximum iterations for the multi-group solver
    double one_group_tolerance_;            ///< Tolerance for the one-group solver
    int one_group_max_iterations_;         ///< Maximum iterations for the one-group solver
    double fission_source_tolerance_;       ///< Tolerance for the fission source convergence
    double keff_tolerance_;                 ///< Tolerance for the k-effective convergence
    int n_theta_;                          ///< Number of theta points in AngularQuadrature
    int n_phi_;                            ///< Number of phi points in AngularQuadrature

    /**
     * @brief Validates the current settings.
     * 
     * Throws an exception if any setting is invalid.
     */
    void validate() {
        if(multi_group_tolerance_ <= 0.0) {
            throw std::invalid_argument("Multi-group tolerance must be positive.");
        }
        if(multi_group_max_iterations_ <= 0) {
            throw std::invalid_argument("Multi-group max iterations must be positive.");
        }
        if(one_group_tolerance_ <= 0.0) {
            throw std::invalid_argument("One-group tolerance must be positive.");
        }
        if(one_group_max_iterations_ <= 0) {
            throw std::invalid_argument("One-group max iterations must be positive.");
        }
        if(n_theta_ <= 0) {
            throw std::invalid_argument("n_theta must be positive.");
        }
        if(n_phi_ <= 0) {
            throw std::invalid_argument("n_phi must be positive.");
        }
        if(n_theta_ % 2 != 0) {
            throw std::invalid_argument("n_theta must be a multiple of 2.");
        }
        if(n_phi_ % 4 != 0) {
            throw std::invalid_argument("n_phi must be a multiple of 4.");
        }
    }
};

#endif // SETTINGS_HPP