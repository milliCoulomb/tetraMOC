// include/MeshHandler.hpp

#ifndef MESH_HANDLER_H
#define MESH_HANDLER_H

#include <vector>
#include <string>
#include <array>
#include <tuple>
#include <unordered_map>

// Structure to store node coordinates
struct Node {
    double x, y, z;
};

// Structure to store tetrahedral cell connectivity
struct TetraCell {
    std::array<int, 4> node_ids;
    int cell_id;
};

// Structure to represent a face
struct Face {
    int n0, n1, n2; // Node indices
    std::vector<int> adjacent_cell_ids; // Adjacent cell IDs
};

// Custom hash function for tuple<int, int, int>
struct TupleHash {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        std::size_t seed = 0;
        std::hash<int> hasher;
        seed ^= hasher(std::get<0>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(std::get<1>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(std::get<2>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// MeshHandler class to manage mesh data
class MeshHandler {
public:
    MeshHandler() = default;

    // Methods to load mesh data
    bool loadNodes(const std::string& filename);
    bool loadCells(const std::string& filename);
    bool loadFaceConnectivity(const std::string& filename);
    
    // Getters for accessing mesh data
    const std::vector<Node>& getNodes() const { return nodes; }
    const std::vector<TetraCell>& getCells() const { return cells; }
    const std::vector<Face>& getFaces() const { return faces; }

    // Method to retrieve boundary faces (faces with only one adjacent cell)
    std::vector<Face> getBoundaryFaces() const;

    // Method to retrieve the adjacent cell for a given face
    // If only_one is true, assumes the face is a boundary face with only one adjacent cell
    int getFaceAdjacentCell(const Face& face, bool only_one = false) const;

    // Method to find the neighboring cell given current cell and exit face
    int getNeighborCell(int current_cell_id, int exit_face_id) const;

private:
    std::vector<Node> nodes;
    std::vector<TetraCell> cells;
    std::vector<Face> faces; // Changed from tuple to Face struct

    // For efficient face lookup
    std::unordered_map<std::tuple<int, int, int>, std::vector<int>, TupleHash> face_to_cells;
};

#endif // MESH_HANDLER_H
