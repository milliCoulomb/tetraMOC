// AngularQuadrature.hpp
#ifndef ANGULARQUADRATURE_HPP
#define ANGULARQUADRATURE_HPP

#include <vector>
#include "Quadrature.hpp"
#include <utility>

struct Direction {
    double mu;
    double phi;
    double weight;
};

class AngularQuadrature {
public:
    AngularQuadrature(int theta_order, int phi_order);
    ~AngularQuadrature() = default;

    const std::vector<Direction>& getDirections() const;

    void addDirection(const Direction& direction);

private:
    int thetaOrder_;
    int phiOrder_;
    std::vector<Direction> directions_;
    void generateQuadrature();
};

#endif // ANGULARQUADRATURE_HPP
