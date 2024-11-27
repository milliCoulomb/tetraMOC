// Quadrature.hpp
#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <vector>
#include <cmath>

// Namespace pour les fonctions de quadrature
namespace Quadrature {

// Fonction pour calculer la valeur du polynôme de Legendre de degré n en x
double legendre(int n, double x);

// Fonction pour calculer la dérivée du polynôme de Legendre de degré n en x
double dLegendre(int n, double x);

// Fonction pour trouver les racines du polynôme de Legendre de degré n utilisant la méthode de Newton-Raphson
std::vector<double> legendreRoots(int n);

// Fonction pour calculer les poids associés aux racines de Gauss-Legendre
std::vector<double> legendreWeights(const std::vector<double>& roots, int n);

// Fonction pour générer les points et poids de Gauss-Chebyshev de première espèce
std::pair<std::vector<double>, std::vector<double>> gaussChebyshev(int n);

} // namespace Quadrature

#endif // QUADRATURE_HPP
