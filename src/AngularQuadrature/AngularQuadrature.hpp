// AngularQuadrature.hpp
#ifndef ANGULARQUADRATURE_HPP
#define ANGULARQUADRATURE_HPP

#include <vector>
#include "Quadrature.hpp"
#include <utility>

namespace SNSolver {

struct Direction {
    double mu;
    double phi;
    double weight;
};

class AngularQuadrature {
public:
    AngularQuadrature(int snOrder);
    ~AngularQuadrature() = default;

    const std::vector<Direction>& getDirections() const;

private:
    int snOrder_;
    std::vector<Direction> directions_;
    void generateQuadrature();
};

} // namespace SNSolver

#endif // ANGULARQUADRATURE_HPP
