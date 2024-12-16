// include/Field.hpp

#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

// Structure to store cell vector field data (e.g., velocity)
struct CellVectorField {
    double vx, vy, vz;
};

// Structure to store cell scalar field data (for future extension)
struct CellScalarField {
    double value;
};

// Enumeration to specify field type
enum class FieldType {
    VECTOR,
    SCALAR
};

// Field class to manage field data
class Field {
public:
    Field() = default;

    // Load a vector field from a file
    bool loadVectorField(const std::string& filename);

    // Load a scalar field from a file (for future extension)
    bool loadScalarField(const std::string& filename);

    // Getters for field data
    const std::vector<CellVectorField>& getVectorFields() const { return vectorFields; }
    const std::vector<CellScalarField>& getScalarFields() const { return scalarFields; }

private:
    std::vector<CellVectorField> vectorFields;
    std::vector<CellScalarField> scalarFields;
};

#endif // FIELD_H