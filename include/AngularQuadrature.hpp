// AngularQuadrature.hpp
#ifndef ANGULARQUADRATURE_HPP
#define ANGULARQUADRATURE_HPP

#include <vector>
#include "Quadrature.hpp"
#include <utility>

/**
 * @struct Direction
 * @brief Structure to store the direction of the quadrature and its weight
 * 
 * @param mu: cosine of the polar angle
 * @param phi: azimuthal angle
 * @param weight: weight of the direction
*/

struct Direction {
    double mu = 0.0; // default value to avoid uninitialized variables
    double phi = 0.0; // default value to avoid uninitialized variables
    double weight = 0.0; // default value to avoid uninitialized variables
};

class AngularQuadrature {
public:
    /**
    @brief Constructor for the AngularQuadrature class

    @param theta_order: number of polar angles
    @param phi_order: number of azimuthal angles
    */
    AngularQuadrature(int theta_order, int phi_order);
    // surcharge constructor to be able to build an angular quadrature with a std vector of directions

    /**
     * @brief Constructor for the AngularQuadrature class
     * 
     * @param directions: vector of directions
     */
    AngularQuadrature(std::vector<Direction> directions);

    
    ~AngularQuadrature() = default;

    /**
     * @brief Get the directions of the quadrature
     * 
     * @return std::vector<Direction> vector of directions
     */
    const std::vector<Direction>& getDirections() const;

    /**
     * @brief add a direction to the quadrature
     * 
     * @param direction: direction to add
     */

    void addDirection(const Direction& direction);

    /**
     * @brief Calculate the total weight of the quadrature
     * 
     * @return double total weight
     */

    double getTotalWeight() const;

private:
    int thetaOrder_;
    int phiOrder_;
    std::vector<Direction> directions_;
    void generateQuadrature();
};

#endif // ANGULARQUADRATURE_HPP
