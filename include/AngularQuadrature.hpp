// AngularQuadrature.hpp
#ifndef ANGULARQUADRATURE_HPP
#define ANGULARQUADRATURE_HPP

#include <vector>
#include "Quadrature.hpp"
#include <utility>

/*
@struct Direction
@brief Structure to store the direction of the quadrature and its weight

@var mu: cosine of the polar angle
@var phi: azimuthal angle
@var weight: weight of the direction
*/

struct Direction {
    double mu = 0.0; // default value to avoid uninitialized variables
    double phi = 0.0; // default value to avoid uninitialized variables
    double weight = 0.0; // default value to avoid uninitialized variables
};

class AngularQuadrature {
public:
    AngularQuadrature(int theta_order, int phi_order);
    // surcharge constructor to be able to build an angular quadrature with a std vector of directions
    AngularQuadrature(std::vector<Direction> directions);
    
    ~AngularQuadrature() = default;

    const std::vector<Direction>& getDirections() const;

    void addDirection(const Direction& direction);

    double getTotalWeight() const;

private:
    int thetaOrder_;
    int phiOrder_;
    std::vector<Direction> directions_;
    void generateQuadrature();
};

#endif // ANGULARQUADRATURE_HPP
