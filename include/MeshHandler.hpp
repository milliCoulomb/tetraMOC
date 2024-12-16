// include/MeshHandler.hpp

#ifndef MESH_HANDLER_H
#define MESH_HANDLER_H

#include <vector>
#include <string>
#include <array>
#include <tuple>

// Structure to store node coordinates
struct Node {
    double x, y, z;
};

// Structure to store tetrahedral cell connectivity
struct TetraCell {
    std::array<int, 4> node_ids;
    int cell_id;
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
    const std::vector<std::tuple<int, int, int, std::vector<int>>>& getFaces() const { return faces; }

private:
    std::vector<Node> nodes;
    std::vector<TetraCell> cells;
    std::vector<std::tuple<int, int, int, std::vector<int>>> faces; // (n0, n1, n2, [cell_ids])
};

#endif // MESH_HANDLER_H