// include/mesh_handler.h
#ifndef MESH_HANDLER_H
#define MESH_HANDLER_H

#include <vector>
#include <string>
#include <array>
#include <tuple>

struct Node {
    double x, y, z;
};

struct TetraCell {
    std::array<int, 4> node_ids;
    int cell_id;
};

class MeshHandler {
public:
    MeshHandler() = default;
    bool loadNodes(const std::string& filename);
    bool loadCells(const std::string& filename);
    bool loadFaceConnectivity(const std::string& filename);

    const std::vector<Node>& getNodes() const { return nodes; }
    const std::vector<TetraCell>& getCells() const { return cells; }
    const std::vector<std::tuple<int, int, int, std::vector<int>>>& getFaces() const { return faces; }

private:
    std::vector<Node> nodes;
    std::vector<TetraCell> cells;
    std::vector<std::tuple<int, int, int, std::vector<int>>> faces; // (n0, n1, n2, [cell_ids])
};

#endif // MESH_HANDLER_H
