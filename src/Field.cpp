// src/field.cpp
#include "Field.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

bool Field::loadField(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open field file " << filename << std::endl;
        return false;
    }
    
    int num_cells;
    infile >> num_cells;
    cellValues.resize(num_cells);
    
    for(int i = 0; i < num_cells; ++i) {
        infile >> cellValues[i].vx >> cellValues[i].vy >> cellValues[i].vz;
        if(infile.fail()) {
            std::cerr << "Error: Invalid format in field file at line " << i+2 << std::endl;
            return false;
        }
    }
    
    infile.close();
    return true;
}
