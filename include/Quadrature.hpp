// Quadrature.hpp
#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <vector>
#include <cmath>

namespace Quadrature {

// values of the Legendre polynomials of degree n at x
double legendre(int n, double x);

// calculates the derivative of the Legendre polynomial of degree n at x
double dLegendre(int n, double x);

// find the roots of the Legendre polynomial of degree n
std::vector<double> legendreRoots(int n);

// generate Gauss-Legendre quadrature points and weights
std::vector<double> legendreWeights(const std::vector<double>& roots, int n);

// generate Gauss-Chebyshev quadrature points and weights
std::pair<std::vector<double>, std::vector<double>> gaussChebyshev(int n);

} // namespace Quadrature

#endif // QUADRATURE_HPP
