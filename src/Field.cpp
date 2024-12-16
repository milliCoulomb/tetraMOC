// src/Field.cpp

#include "Field.hpp"

// Load a vector field (e.g., velocity) from a file
bool Field::loadVectorField(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open vector field file " << filename << std::endl;
        return false;
    }

    int num_cells;
    infile >> num_cells;
    if (infile.fail()) {
        std::cerr << "Error: Unable to read number of vectors from " << filename << std::endl;
        return false;
    }
    vectorFields.resize(num_cells);

    std::string line;
    std::getline(infile, line); // Consume the rest of the first line

    for(int i = 0; i < num_cells; ++i) {
        if (!std::getline(infile, line)) {
            std::cerr << "Error: Unexpected end of vector field file at line " << i + 2 << std::endl;
            return false;
        }
        std::istringstream iss(line);
        iss >> vectorFields[i].vx >> vectorFields[i].vy >> vectorFields[i].vz;
        if(iss.fail()) {
            std::cerr << "Error: Invalid format in vector field file at line " << i + 2 << std::endl;
            return false;
        }
        // Check for extra data on the line
        double extra;
        if(iss >> extra) {
            std::cerr << "Error: Extra data found in vector field file at line " << i + 2 << ": " << extra << std::endl;
            return false;
        }
    }

    infile.close();
    return true;
}

// Load a scalar field from a file (for future extension)
bool Field::loadScalarField(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open scalar field file " << filename << std::endl;
        return false;
    }

    int num_cells;
    infile >> num_cells;
    if (infile.fail()) {
        std::cerr << "Error: Unable to read number of scalar values from " << filename << std::endl;
        return false;
    }
    scalarFields.resize(num_cells);

    std::string line;
    std::getline(infile, line); // Consume the rest of the first line

    for(int i = 0; i < num_cells; ++i) {
        if (!std::getline(infile, line)) {
            std::cerr << "Error: Unexpected end of scalar field file at line " << i + 2 << std::endl;
            return false;
        }
        std::istringstream iss(line);
        iss >> scalarFields[i].value;
        if(iss.fail()) {
            std::cerr << "Error: Invalid format in scalar field file at line " << i + 2 << std::endl;
            return false;
        }
        // Check for extra data on the line
        double extra;
        if(iss >> extra) {
            std::cerr << "Error: Extra data found in scalar field file at line " << i + 2 << ": " << extra << std::endl;
            return false;
        }
    }

    infile.close();
    return true;
}
