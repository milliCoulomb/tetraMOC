// include/Field.hpp

#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory> // For std::shared_ptr
#include "Vector3D.hpp"
#include "Logger.hpp"

// Structure to store cell scalar field data (for future extension)
struct CellScalarField {
    double value;
};

// Enumeration to specify field type
enum class FieldType {
    VECTOR,
    SCALAR
};

// Type alias to use Vector3D for vector fields
using CellVectorField = Vector3D;

// Field class to manage field data
class Field {
public:
    Field() = default;

    // Constructor with shared CellVectorField and direction vector
    Field(std::shared_ptr<std::vector<CellVectorField>> shared_vector_fields, const Vector3D& direction)
        : sharedVectorFields(shared_vector_fields), direction_(direction) {}

    // Load a vector field from a file (initializes shared CellVectorField)
    bool loadVectorField(const std::string& filename);

    // Load a scalar field from a file (for future extension)
    bool loadScalarField(const std::string& filename);

    // Getters for field data
    const std::vector<CellVectorField>& getVectorFields() const { return *sharedVectorFields; }
    const std::vector<CellScalarField>& getScalarFields() const { return scalarFields; }

    // Getter and Setter for direction
    const Vector3D& getDirection() const { return direction_; }
    void setDirection(const Vector3D& direction) { direction_ = direction; }

private:
    // Shared pointer to CellVectorField to enable sharing across Field instances
    std::shared_ptr<std::vector<CellVectorField>> sharedVectorFields;

    std::vector<CellScalarField> scalarFields;

    Vector3D direction_; // Separate direction vector
};

#endif // FIELD_H