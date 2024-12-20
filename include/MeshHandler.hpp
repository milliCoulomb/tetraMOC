// include/MeshHandler.hpp

#ifndef MESH_HANDLER_H
#define MESH_HANDLER_H

#include <vector>
#include <string>
#include <array>
#include <tuple>
#include "Vector3D.hpp"
#include <unordered_map>

// Structure to store node coordinates
// struct Node {
//     double x, y, z;
// };

// Structure to store tetrahedral cell connectivity
struct TetraCell {
    std::array<int, 4> node_ids = {0, 0, 0, 0};
    int cell_id = -1;

    // Constructors
    TetraCell() = default;
    TetraCell(const std::array<int, 4>& nodes, int id)
        : node_ids(nodes), cell_id(id) {}

    // Defaulted copy and move constructors and assignments
    TetraCell(const TetraCell&) = default;
    TetraCell(TetraCell&&) = default;
    TetraCell& operator=(const TetraCell&) = default;
    TetraCell& operator=(TetraCell&&) = default;
};

// Structure to represent a face
struct Face {
    int n0 = 0, n1 = 0, n2 = 0; // Node indices
    std::vector<int> adjacent_cell_ids; // Adjacent cell IDs

    // Constructors
    Face() = default;
    Face(int node0, int node1, int node2, const std::vector<int>& adj_cells)
        : n0(node0), n1(node1), n2(node2), adjacent_cell_ids(adj_cells) {}

    // Defaulted copy and move constructors and assignments
    Face(const Face&) = default;
    Face(Face&&) = default;
    Face& operator=(const Face&) = default;
    Face& operator=(Face&&) = default;
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
    const std::vector<Vector3D>& getNodes() const { return nodes; }
    const std::vector<TetraCell>& getCells() const { return cells; }
    const std::vector<Face>& getFaces() const { return faces; }

    // Method to retrieve boundary faces (faces with only one adjacent cell)
    std::vector<Face> getBoundaryFaces() const;

    // Method to retrieve the adjacent cell for a given face
    // If only_one is true, assumes the face is a boundary face with only one adjacent cell
    int getFaceAdjacentCell(const Face& face, bool only_one = false) const;

    // Method to find the neighboring cell given current cell and exit face
    int getNeighborCell(int current_cell_id, int exit_face_id) const;

    // calculate the centroid of a tetrahedron
    Vector3D getCellCenter(int cell_id) const;

private:
    std::vector<Vector3D> nodes;
    std::vector<TetraCell> cells;
    std::vector<Face> faces; // Changed from tuple to Face struct

    // For efficient face lookup
    std::unordered_map<std::tuple<int, int, int>, std::vector<int>, TupleHash> face_to_cells;
};

#endif // MESH_HANDLER_H
