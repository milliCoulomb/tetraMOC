// src/MeshHandler.cpp

#include "MeshHandler.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

// Load nodes from nodes.txt
bool MeshHandler::loadNodes(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open nodes file " << filename << std::endl;
        return false;
    }
    
    int num_nodes;
    infile >> num_nodes;
    if (infile.fail()) {
        std::cerr << "Error: Unable to read number of nodes from " << filename << std::endl;
        return false;
    }
    nodes.resize(num_nodes);
    
    for(int i = 0; i < num_nodes; ++i) {
        infile >> nodes[i].x >> nodes[i].y >> nodes[i].z;
        if(infile.fail()) {
            std::cerr << "Error: Invalid format in nodes file at line " << i+2 << std::endl;
            return false;
        }
    }
    
    infile.close();
    return true;
}

// Load cells from cells.txt
bool MeshHandler::loadCells(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open cells file " << filename << std::endl;
        return false;
    }
    
    int num_cells;
    infile >> num_cells;
    if (infile.fail()) {
        std::cerr << "Error: Unable to read number of cells from " << filename << std::endl;
        return false;
    }
    cells.resize(num_cells);
    
    std::string line;
    std::getline(infile, line); // Consume the rest of the first line
    
    for(int i = 0; i < num_cells; ++i) {
        if (!std::getline(infile, line)) {
            std::cerr << "Error: Unexpected end of cells file at line " << i + 2 << std::endl;
            return false;
        }
        std::istringstream iss(line);
        iss >> cells[i].node_ids[0] >> cells[i].node_ids[1] >> cells[i].node_ids[2] >> cells[i].node_ids[3];
        cells[i].cell_id = i;
        if(iss.fail()) {
            std::cerr << "Error: Invalid format in cells file at line " << i + 2 << std::endl;
            return false;
        }
        // Check for extra node IDs
        int extra;
        if(iss >> extra) {
            std::cerr << "Error: Extra node ID " << extra << " found in cells file at line " << i + 2 << std::endl;
            return false;
        }
    }
    
    infile.close();
    return true;
}

bool MeshHandler::loadFaceConnectivity(const std::string& filename) {
    std::ifstream infile(filename);
    if(!infile.is_open()) {
        std::cerr << "Error: Unable to open face connectivity file " << filename << std::endl;
        return false;
    }
    
    std::string line;
    // Read the first line which contains the number of faces
    int num_faces;
    if (!(infile >> num_faces)) {
        std::cerr << "Error: Invalid format in faces file (missing number of faces)" << std::endl;
        return false;
    }
    std::getline(infile, line); // consume the rest of the first line
    
    for(int i = 0; i < num_faces; ++i) {
        if (!std::getline(infile, line)) {
            std::cerr << "Error: Unexpected end of faces file at line " << i + 2 << std::endl;
            return false;
        }
        std::istringstream iss(line);
        int n0, n1, n2;
        std::vector<int> cells;
        if(!(iss >> n0 >> n1 >> n2)) {
            std::cerr << "Error: Invalid face format at line " << i + 2 << ": " << line << std::endl;
            return false;
        }
        int cell_id;
        while(iss >> cell_id) {
            cells.push_back(cell_id);
        }
        // Ensure that at least one cell ID is present
        if(cells.empty()) {
            std::cerr << "Error: No cell IDs found for face at line " << i + 2 << ": " << line << std::endl;
            return false;
        }
        // Ensure the face is sorted
        std::array<int, 3> sorted_nodes = {n0, n1, n2};
        std::sort(sorted_nodes.begin(), sorted_nodes.end());
        faces.emplace_back(sorted_nodes[0], sorted_nodes[1], sorted_nodes[2], cells);
    }
    
    infile.close();
    return true;
}

// Load velocity field from field.txt
bool MeshHandler::loadVelocityField(const std::string& filename) {
    std::ifstream infile(filename);
    if(!infile.is_open()) {
        std::cerr << "Error: Unable to open velocity field file " << filename << std::endl;
        return false;
    }
    
    int num_vectors;
    if(!(infile >> num_vectors)) {
        std::cerr << "Error: Unable to read number of velocity vectors from " << filename << std::endl;
        return false;
    }
    
    // Check if the number of vectors matches the number of cells
    if(num_vectors != static_cast<int>(cells.size())) {
        std::cerr << "Error: Number of velocity vectors (" << num_vectors << ") does not match number of cells (" << cells.size() << ")." << std::endl;
        return false;
    }
    
    velocities.resize(num_vectors);
    
    for(int i = 0; i < num_vectors; ++i) {
        infile >> velocities[i].vx >> velocities[i].vy >> velocities[i].vz;
        if(infile.fail()) {
            std::cerr << "Error: Invalid format in velocity field file at line " << i + 2 << std::endl;
            return false;
        }
    }
    
    // Optionally, check for extra data after expected velocity vectors
    double extra;
    if(infile >> extra) {
        std::cerr << "Error: Extra data found in velocity field file after expected velocity vectors." << std::endl;
        return false;
    }
    
    infile.close();
    return true;
}