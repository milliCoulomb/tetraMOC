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
          max_power_iterations_(100),
          keff_tolerance_(1e-6) {}

    /**
     * @brief Parameterized constructor to initialize all settings.
     * 
     * @param multi_group_tolerance       Tolerance for the multi-group solver.
     * @param multi_group_max_iterations Max iterations for the multi-group solver.
     * @param one_group_tolerance         Tolerance for the one-group solver.
     * @param one_group_max_iterations   Max iterations for the one-group solver.
     */
    Settings(double multi_group_tolerance, int multi_group_max_iterations,
             double one_group_tolerance, int one_group_max_iterations,
             double fission_source_tolerance, int max_power_iterations, 
             double keff_tolerance)
        : multi_group_tolerance_(multi_group_tolerance),
          multi_group_max_iterations_(multi_group_max_iterations),
          one_group_tolerance_(one_group_tolerance),
          one_group_max_iterations_(one_group_max_iterations),
          fission_source_tolerance_(fission_source_tolerance),
          max_power_iterations_(max_power_iterations),
          keff_tolerance_(keff_tolerance)
    {
        validate();
    }

    // Getters
    double getMultiGroupTolerance() const { return multi_group_tolerance_; }
    int getMultiGroupMaxIterations() const { return multi_group_max_iterations_; }
    double getOneGroupTolerance() const { return one_group_tolerance_; }
    int getOneGroupMaxIterations() const { return one_group_max_iterations_; }
    double getFissionSourceTolerance() const { return fission_source_tolerance_; }
    int getMaxPowerIterations() const { return max_power_iterations_; }
    double getKeffTolerance() const { return keff_tolerance_; }

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

    void setMaxPowerIterations(int max_iterations) { 
        if(max_iterations <= 0) {
            throw std::invalid_argument("Max power iterations must be positive.");
        }
        max_power_iterations_ = max_iterations; 
    }

    void setKeffTolerance(double tolerance) { 
        if(tolerance <= 0.0) {
            throw std::invalid_argument("k-effective tolerance must be positive.");
        }
        keff_tolerance_ = tolerance; 
    }

private:
    // Member variables storing settings
    double multi_group_tolerance_;          ///< Tolerance for the multi-group solver
    int multi_group_max_iterations_;       ///< Maximum iterations for the multi-group solver
    double one_group_tolerance_;            ///< Tolerance for the one-group solver
    int one_group_max_iterations_;         ///< Maximum iterations for the one-group solver
    double fission_source_tolerance_;       ///< Tolerance for the fission source convergence
    int max_power_iterations_;             ///< Maximum iterations for the power iteration method
    double keff_tolerance_;                 ///< Tolerance for the k-effective convergence

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
        if(fission_source_tolerance_ <= 0.0) {
            throw std::invalid_argument("Fission source tolerance must be positive.");
        }
        if(max_power_iterations_ <= 0) {
            throw std::invalid_argument("Max power iterations must be positive.");
        }
        if(keff_tolerance_ <= 0.0) {
            throw std::invalid_argument("k-effective tolerance must be positive.");
        }
    }
};

#endif // SETTINGS_HPP