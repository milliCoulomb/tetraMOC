// src/MeshHandler.cpp

#include "MeshHandler.hpp"
#include "Vector3D.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream> // For debugging purposes

// Implementation of loadNodes
bool MeshHandler::loadNodes(const std::string& filename) {
    std::ifstream infile(filename);
    if(!infile.is_open()) {
        std::cerr << "Error: Unable to open nodes file: " << filename << std::endl;
        return false;
    }
    
    int num_nodes;
    infile >> num_nodes;
    if(infile.fail()) {
        std::cerr << "Error: Failed to read the number of nodes from: " << filename << std::endl;
        return false;
    }
    
    nodes.reserve(num_nodes);
    
    for(int i = 0; i < num_nodes; ++i) {
        int node_id;
        double x, y, z;
        infile >> node_id >> x >> y >> z;
        if(infile.fail()) {
            std::cerr << "Error: Invalid node format in line: " << (i + 2) << std::endl;
            return false;
        }
        nodes.emplace_back(Vector3D(x, y, z));
    }
    
    // Optional: Verify no extra data
    std::string extra;
    if(infile >> extra) {
        std::cerr << "Error: Extra data found in nodes file after expected nodes." << std::endl;
        return false;
    }
    
    return true;
}

bool MeshHandler::loadCells(const std::string& filename) {
    std::ifstream infile(filename);
    if(!infile.is_open()) {
        std::cerr << "Error: Unable to open cells file: " << filename << std::endl;
        return false;
    }
    
    int num_cells;
    infile >> num_cells;
    if(infile.fail()) {
        std::cerr << "Error: Failed to read the number of cells from: " << filename << std::endl;
        return false;
    }
    
    cells.reserve(num_cells);
    
    for(int i = 0; i < num_cells; ++i) {
        int cell_id, n0, n1, n2, n3;
        infile >> cell_id >> n0 >> n1 >> n2 >> n3;
        if(infile.fail()) {
            std::cerr << "Error: Invalid cell format in line: " << (i + 2) << std::endl;
            return false;
        }
        cells.emplace_back(TetraCell{ {n0, n1, n2, n3}, cell_id });
    }
    
    // Optional: Verify no extra data
    std::string extra;
    if(infile >> extra) {
        std::cerr << "Error: Extra data found in cells file after expected cells." << std::endl;
        return false;
    }
    
    return true;
}

// Implementation of loadFaceConnectivity
bool MeshHandler::loadFaceConnectivity(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open face connectivity file " << filename << std::endl;
        return false;
    }

    faces.clear();
    face_to_cells.clear();

    int num_faces;
    infile >> num_faces;
    if(infile.fail()){
        std::cerr << "Error: Unable to read number of faces from " << filename << std::endl;
        return false;
    }

    for(int i = 0; i < num_faces; ++i){
        Face face;
        int n0, n1, n2, count;
        infile >> n0 >> n1 >> n2 >> count;
        if(infile.fail() || count < 1){
            std::cerr << "Error: Invalid face format in line: " << (i + 2) << std::endl;
            return false;
        }

        face.n0 = n0;
        face.n1 = n1;
        face.n2 = n2;

        for(int j = 0; j < count; ++j){
            int cell_id;
            infile >> cell_id;
            if(infile.fail()){
                std::cerr << "Error: Unable to read cell_id for face in line: " << (i + 2) << std::endl;
                return false;
            }
            face.adjacent_cell_ids.push_back(cell_id);
        }

        faces.push_back(face);

        // Sort the node indices to create a unique key
        std::array<int, 3> sorted_nodes = {face.n0, face.n1, face.n2};
        std::sort(sorted_nodes.begin(), sorted_nodes.end());
        std::tuple<int, int, int> face_key = std::make_tuple(sorted_nodes[0], sorted_nodes[1], sorted_nodes[2]);

        face_to_cells[face_key] = face.adjacent_cell_ids;
    }

    // Check for extra data
    std::string extra;
    if(infile >> extra){
        std::cerr << "Error: Extra data found in faces file after expected faces." << std::endl;
        return false;
    }

    infile.close();
    return true;
}

// Implementation of getBoundaryFaces
std::vector<Face> MeshHandler::getBoundaryFaces() const {
    std::vector<Face> boundary_faces;
    for(const auto& face : faces) {
        if(face.adjacent_cell_ids.size() == 1) {
            boundary_faces.push_back(face);
        }
    }
    return boundary_faces;
}

// Implementation of getFaceAdjacentCell
int MeshHandler::getFaceAdjacentCell(const Face& face, bool only_one) const {
    // Sort the node indices to create a unique key
    std::array<int, 3> sorted_nodes = {face.n0, face.n1, face.n2};
    std::sort(sorted_nodes.begin(), sorted_nodes.end());
    std::tuple<int, int, int> face_key = std::make_tuple(sorted_nodes[0], sorted_nodes[1], sorted_nodes[2]);

    auto it = face_to_cells.find(face_key);
    if(it == face_to_cells.end()) {
        return -1; // Face not found
    }

    const std::vector<int>& adj_cells = it->second;

    if(only_one) {
        if(adj_cells.size() == 1) {
            return adj_cells[0];
        } else {
            return -1; // Not a boundary face
        }
    } else {
        if(!adj_cells.empty()) {
            return adj_cells[0]; // Return the first adjacent cell
        } else {
            return -1; // No adjacent cells
        }
    }
}

// Implementation of getNeighborCell
int MeshHandler::getNeighborCell(int current_cell_id, int exit_face_id) const {
    // Validate current_cell_id
    if(current_cell_id < 0 || current_cell_id >= static_cast<int>(cells.size())) {
        std::cerr << "Error: Invalid current cell ID: " << current_cell_id << std::endl;
        return -1;
    }

    const TetraCell& cell = cells[current_cell_id];
    std::array<int, 3> face_nodes;

    // Assuming exit_face_id ranges from 0 to 3 for tetrahedral cells
    switch(exit_face_id) {
        case 0:
            face_nodes = {cell.node_ids[0], cell.node_ids[1], cell.node_ids[2]};
            break;
        case 1:
            face_nodes = {cell.node_ids[0], cell.node_ids[1], cell.node_ids[3]};
            break;
        case 2:
            face_nodes = {cell.node_ids[0], cell.node_ids[2], cell.node_ids[3]};
            break;
        case 3:
            face_nodes = {cell.node_ids[1], cell.node_ids[2], cell.node_ids[3]};
            break;
        default:
            std::cerr << "Error: Invalid exit face ID: " << exit_face_id << std::endl;
            return -1;
    }

    // Sort the node indices to create a unique key
    std::sort(face_nodes.begin(), face_nodes.end());
    std::tuple<int, int, int> face_key = std::make_tuple(face_nodes[0], face_nodes[1], face_nodes[2]);

    auto it = face_to_cells.find(face_key);
    if(it == face_to_cells.end()) {
        return -1; // Face not found
    }

    const std::vector<int>& adj_cells = it->second;

    // Find the cell adjacent to the current cell
    for(auto cell_id : adj_cells) {
        if(cell_id != current_cell_id) {
            return cell_id;
        }
    }

    return -1; // No neighboring cell found
}

Vector3D MeshHandler::getCellCenter(int cell_id) const {
    if(cell_id < 0 || cell_id >= static_cast<int>(cells.size())) {
        throw std::out_of_range("Invalid cell ID");
    }

    const TetraCell& cell = cells[cell_id];
    double x = 0.0, y = 0.0, z = 0.0;
    for(int i = 0; i < 4; ++i) {
        x += nodes[cell.node_ids[i]].x;
        y += nodes[cell.node_ids[i]].y;
        z += nodes[cell.node_ids[i]].z;
    }
    return Vector3D(x / 4.0, y / 4.0, z / 4.0);
}