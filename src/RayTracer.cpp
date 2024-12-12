// src/ray_tracer.cpp
#include "RayTracer.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

RayTracer::RayTracer(const MeshHandler& mesh_handler, const Field& field_handler)
    : mesh(mesh_handler), field(field_handler) {
    // Initialization if needed
}

bool RayTracer::loadFaceConnectivity(const std::string& filename) {
    std::ifstream infile(filename);
    if(!infile.is_open()) {
        std::cerr << "Error: Unable to open face connectivity file " << filename << std::endl;
        return false;
    }

    std::string line;
    while(std::getline(infile, line)) {
        std::istringstream iss(line);
        int n0, n1, n2;
        std::vector<int> cells;
        if(!(iss >> n0 >> n1 >> n2)) {
            std::cerr << "Error: Invalid format in faces file at line: " << line << std::endl;
            continue;
        }
        int cell_id;
        while(iss >> cell_id) {
            cells.push_back(cell_id);
        }
        // Ensure the face is sorted
        std::array<int, 3> sorted_nodes = {n0, n1, n2};
        std::sort(sorted_nodes.begin(), sorted_nodes.end());
        std::tuple<int, int, int> face_key = std::make_tuple(sorted_nodes[0], sorted_nodes[1], sorted_nodes[2]);
        face_to_cells[face_key] = cells;
    }

    infile.close();
    return true;
}

int RayTracer::getNeighborCell(int current_cell_id, int exit_face_id) const {
    // Retrieve the current cell's node IDs
    const TetraCell& cell = mesh.getCells()[current_cell_id];
    std::array<int, 3> face_nodes;
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
            return -1;
    }

    // Sort the face nodes to match the key
    std::sort(face_nodes.begin(), face_nodes.end());
    std::tuple<int, int, int> face_key = std::make_tuple(face_nodes[0], face_nodes[1], face_nodes[2]);

    auto it = face_to_cells.find(face_key);
    if(it == face_to_cells.end()) {
        return -1; // No neighbor found
    }

    const std::vector<int>& adjacent_cells = it->second;
    for(auto cell_id : adjacent_cells) {
        if(cell_id != current_cell_id) {
            return cell_id;
        }
    }

    return -1; // No neighbor found
}

std::vector<DirectionData> RayTracer::traceRay(int start_cell_id, const std::array<double, 3>& start_point, int max_iter) {
    std::vector<DirectionData> pathline;
    int current_cell_id = start_cell_id;
    std::array<double, 3> current_point = start_point;
    double total_time = 0.0;

    for(int iter = 0; iter < max_iter; ++iter) {
        if(current_cell_id < 0 || current_cell_id >= static_cast<int>(mesh.getCells().size())) {
            std::cerr << "Invalid current cell ID: " << current_cell_id << std::endl;
            break;
        }

        const TetraCell& cell = mesh.getCells()[current_cell_id];
        const CellField& field_val = field.getCellValues()[current_cell_id];
        Tetrahedron tetra(cell, mesh.getNodes(), field_val);

        SNSolver::Vector3D v(field_val.vx, field_val.vy, field_val.vz);
        double t_exit;
        std::array<double, 3> x_exit;
        int exit_face_id;

        bool has_exit = tetra.findExit(current_point, v, t_exit, x_exit, exit_face_id);
        if(!has_exit) {
            std::cerr << "No exit found for cell " << current_cell_id << " at iteration " << iter << std::endl;
            break;
        }

        // Record the segment
        DirectionData segment;
        segment.cell_id = current_cell_id;
        segment.time_spent = t_exit;
        segment.start_point = current_point;
        segment.end_point = x_exit;
        pathline.push_back(segment);
        total_time += t_exit;

        // Find the neighboring cell via the exit face
        int neighbor_cell_id = getNeighborCell(current_cell_id, exit_face_id);
        if(neighbor_cell_id == -1) {
            // No neighboring cell found; ray exits the domain
            std::cout << "Ray exited the domain at iteration " << iter << std::endl;
            break;
        }

        current_cell_id = neighbor_cell_id;
        current_point = x_exit;
    }

    std::cout << "Ray tracing completed. Total time elapsed: " << total_time << std::endl;
    return pathline;
}