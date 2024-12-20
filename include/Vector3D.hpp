// include/Vector3D.hpp
#ifndef VECTOR3D_HPP
#define VECTOR3D_HPP

#include <cmath>
#include <iostream>

class Vector3D {
public:
    double x, y, z;

    // default constructor
    Vector3D() : x(0.0), y(0.0), z(0.0) {}

    // with 3 scalars
    Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // constructor from a scalar
    Vector3D(const Vector3D& other) = default;

    // assign operator (par défaut)
    Vector3D& operator=(const Vector3D& other) = default;

    // addition of vectors
    Vector3D operator+(const Vector3D& v) const {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }

    // subtraction of vectors
    Vector3D operator-(const Vector3D& v) const {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }

    Vector3D operator-() const { return Vector3D(-x, -y, -z); }

    // right multiplication by scalar
    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    // left multiplication by scalar
    friend Vector3D operator*(double scalar, const Vector3D& v) {
        return v * scalar;
    }

    // division by scalar
    Vector3D operator/(double scalar) const {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }

    // scalar product
    double dot(const Vector3D& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    // cross product
    Vector3D cross(const Vector3D& v) const {
        return Vector3D(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }

    // calculate the norm of the vector
    double norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // normalisation
    Vector3D normalized() const {
        double n = norm();
        if (n == 0.0) {
            // std::cerr << "Warning: Attempt to normalize a zero vector.\n";
            return Vector3D(0.0, 0.0, 0.0);
        }
        return (*this) / n;
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector3D& v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }

    // test of equality
    bool operator==(const Vector3D& v) const {
        const double epsilon = 1e-12;
        return (std::abs(x - v.x) < epsilon) &&
               (std::abs(y - v.y) < epsilon) &&
               (std::abs(z - v.z) < epsilon);
    }

    // Inégalité
    bool operator!=(const Vector3D& v) const {
        return !(*this == v);
    }

    bool isAlmostEqual(const Vector3D& v, double tol = 1e-6) const {
        return ((*this) - v).norm() < tol;
    }
};

#endif // VECTOR3D_HPP
