// tests/test_MeshHandler.cpp
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio> // For std::remove
#include "MeshHandler.hpp"
#include "TestUtils.hpp"

using namespace SNSolver;

// Test Fixture for MeshHandler
class MeshHandlerTest : public ::testing::Test {
protected:
    // Temporary file names
    std::string nodes_file = "temp_nodes.txt";
    std::string cells_file = "temp_cells.txt";
    std::string faces_file = "temp_faces.txt";
    std::string field_file = "temp_field.txt";

    // Clean up temporary files after each test
    void TearDown() override {
        std::remove(nodes_file.c_str());
        std::remove(cells_file.c_str());
        std::remove(faces_file.c_str());
        std::remove(field_file.c_str());
    }
};

TEST_F(MeshHandlerTest, LoadNodes_Success) {
    // Define content for nodes.txt with node IDs
    std::string nodes_content = "4\n"
                                 "0 0.0 0.0 0.0\n"
                                 "1 1.0 0.0 0.0\n"
                                 "2 0.0 1.0 0.0\n"
                                 "3 0.0 0.0 1.0\n";
    
    // Create temporary nodes.txt
    ASSERT_TRUE(createTempFile(nodes_file, nodes_content)) << "Failed to create temporary nodes.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadNodes(nodes_file)) << "MeshHandler failed to load nodes.txt";

    // Verify loaded nodes
    const std::vector<Node>& nodes = mesh.getNodes();
    ASSERT_EQ(nodes.size(), 4) << "Number of loaded nodes mismatch";

    EXPECT_DOUBLE_EQ(nodes[0].x, 0.0);
    EXPECT_DOUBLE_EQ(nodes[0].y, 0.0);
    EXPECT_DOUBLE_EQ(nodes[0].z, 0.0);

    EXPECT_DOUBLE_EQ(nodes[1].x, 1.0);
    EXPECT_DOUBLE_EQ(nodes[1].y, 0.0);
    EXPECT_DOUBLE_EQ(nodes[1].z, 0.0);

    EXPECT_DOUBLE_EQ(nodes[2].x, 0.0);
    EXPECT_DOUBLE_EQ(nodes[2].y, 1.0);
    EXPECT_DOUBLE_EQ(nodes[2].z, 0.0);

    EXPECT_DOUBLE_EQ(nodes[3].x, 0.0);
    EXPECT_DOUBLE_EQ(nodes[3].y, 0.0);
    EXPECT_DOUBLE_EQ(nodes[3].z, 1.0);
}

TEST_F(MeshHandlerTest, LoadNodes_MalformedFile) {
    // Define malformed content for nodes.txt (missing a coordinate)
    std::string nodes_content = "3\n"
                                 "0.0 0.0 0.0\n"
                                 "1.0 0.0\n" // Incomplete line
                                 "0.0 1.0 0.0\n";
    
    // Create temporary nodes.txt
    ASSERT_TRUE(createTempFile(nodes_file, nodes_content)) << "Failed to create temporary nodes.txt";

    MeshHandler mesh;
    EXPECT_FALSE(mesh.loadNodes(nodes_file)) << "MeshHandler should fail to load malformed nodes.txt";
}

TEST_F(MeshHandlerTest, LoadNodes_IncorrectNodeCount) {
    // Define content where the number of nodes specified doesn't match the actual data
    std::string nodes_content = "5\n" // Specifying 5 nodes
                                 "0.0 0.0 0.0\n"
                                 "1.0 0.0 0.0\n"
                                 "0.0 1.0 0.0\n"
                                 "0.0 0.0 1.0\n"; // Only 4 nodes provided
    
    // Create temporary nodes.txt
    ASSERT_TRUE(createTempFile(nodes_file, nodes_content)) << "Failed to create temporary nodes.txt";

    MeshHandler mesh;
    EXPECT_FALSE(mesh.loadNodes(nodes_file)) << "MeshHandler should fail due to node count mismatch";
}

TEST_F(MeshHandlerTest, LoadCells_Success) {
    // Define content for cells.txt with cell IDs and four node IDs per cell
    std::string cells_content = "2\n"
                                 "0 0 1 2 3\n"
                                 "1 1 2 3 4\n"; // Assuming node IDs up to 4
    
    // Create temporary cells.txt
    ASSERT_TRUE(createTempFile(cells_file, cells_content)) << "Failed to create temporary cells.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadCells(cells_file)) << "MeshHandler failed to load cells.txt";

    // Verify loaded cells
    const std::vector<TetraCell>& cells = mesh.getCells();
    ASSERT_EQ(cells.size(), 2) << "Number of loaded cells mismatch";

    EXPECT_EQ(cells[0].node_ids[0], 0);
    EXPECT_EQ(cells[0].node_ids[1], 1);
    EXPECT_EQ(cells[0].node_ids[2], 2);
    EXPECT_EQ(cells[0].node_ids[3], 3);
    EXPECT_EQ(cells[0].cell_id, 0);

    EXPECT_EQ(cells[1].node_ids[0], 1);
    EXPECT_EQ(cells[1].node_ids[1], 2);
    EXPECT_EQ(cells[1].node_ids[2], 3);
    EXPECT_EQ(cells[1].node_ids[3], 4);
    EXPECT_EQ(cells[1].cell_id, 1);
}

TEST_F(MeshHandlerTest, LoadCells_MalformedFile) {
    // Define malformed content for cells.txt (missing a node ID)
    std::string cells_content = "1\n"
                                 "0 1 2\n"; // Incomplete cell definition
    
    // Create temporary cells.txt
    ASSERT_TRUE(createTempFile(cells_file, cells_content)) << "Failed to create temporary cells.txt";

    MeshHandler mesh;
    EXPECT_FALSE(mesh.loadCells(cells_file)) << "MeshHandler should fail to load malformed cells.txt";
}

TEST_F(MeshHandlerTest, LoadCells_NonTetrahedral) {
    // Define content for cells.txt with cells not having 4 nodes
    std::string cells_content = "2\n"
                                 "0 1 2 3 4\n" // 5 nodes instead of 4
                                 "1 2 3 4\n";
    
    // Create temporary cells.txt
    ASSERT_TRUE(createTempFile(cells_file, cells_content)) << "Failed to create temporary cells.txt";

    MeshHandler mesh;
    // Assuming MeshHandler is designed to check for exactly 4 nodes per cell
    EXPECT_FALSE(mesh.loadCells(cells_file)) << "MeshHandler should fail for non-tetrahedral cells";
}

TEST_F(MeshHandlerTest, LoadFaces_Success) {
    // Define content for faces.txt with counts of adjacent cells
    std::string faces_content = "3\n" // Number of faces
                                   "0 1 2 2 0 1\n" // Face with nodes 0,1,2 adjacent to cells 0 and 1
                                   "0 1 3 1 0\n"   // Face with nodes 0,1,3 adjacent to cell 0
                                   "1 2 3 1 1\n";  // Face with nodes 1,2,3 adjacent to cell 1
    
    // Create temporary faces.txt
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler failed to load faces.txt";

    // Verify loaded faces
    const auto& faces = mesh.getFaces();
    ASSERT_EQ(faces.size(), 3) << "Number of loaded faces mismatch";

    // Face 0
    EXPECT_EQ(faces[0].n0, 0);
    EXPECT_EQ(faces[0].n1, 1);
    EXPECT_EQ(faces[0].n2, 2);
    ASSERT_EQ(faces[0].adjacent_cell_ids.size(), 2);
    EXPECT_EQ(faces[0].adjacent_cell_ids[0], 0);
    EXPECT_EQ(faces[0].adjacent_cell_ids[1], 1);

    // Face 1
    EXPECT_EQ(faces[1].n0, 0);
    EXPECT_EQ(faces[1].n1, 1);
    EXPECT_EQ(faces[1].n2, 3);
    ASSERT_EQ(faces[1].adjacent_cell_ids.size(), 1);
    EXPECT_EQ(faces[1].adjacent_cell_ids[0], 0);

    // Face 2
    EXPECT_EQ(faces[2].n0, 1);
    EXPECT_EQ(faces[2].n1, 2);
    EXPECT_EQ(faces[2].n2, 3);
    ASSERT_EQ(faces[2].adjacent_cell_ids.size(), 1);
    EXPECT_EQ(faces[2].adjacent_cell_ids[0], 1);
}

TEST_F(MeshHandlerTest, LoadFaces_MalformedFile) {
    // Define malformed content for faces.txt (missing cell IDs)
    std::string faces_content = "2\n"
                                 "0 1 2\n" // Missing cell IDs
                                 "1 2 3 1\n";
    
    // Create temporary faces.txt
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";

    MeshHandler mesh;
    EXPECT_FALSE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler should fail to load malformed faces.txt";
}

TEST_F(MeshHandlerTest, LoadFaces_IncorrectFaceCount) {
    // Define content where the number of faces specified doesn't match the actual data
    std::string faces_content = "4\n" // Specifying 4 faces
                                   "0 1 2 0 1\n"
                                   "0 1 3 0\n"
                                   "1 2 3 1\n"; // Only 3 faces provided
    
    // Create temporary faces.txt
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";

    MeshHandler mesh;
    EXPECT_FALSE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler should fail due to face count mismatch";
}

TEST_F(MeshHandlerTest, LoadFaces_BeforeNodesAndCells) {
    // Define content for faces.txt referencing node and cell IDs that haven't been loaded yet
    std::string faces_content = "1\n"
                                 "0 1 2 0\n";
    
    // Create temporary faces.txt
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler should load faces even if nodes and cells are not loaded yet";
    
    // Access face data, but nodes and cells may not exist
    const auto& faces = mesh.getFaces();
    ASSERT_EQ(faces.size(), 1) << "Number of loaded faces mismatch";
    EXPECT_EQ(faces[0].adjacent_cell_ids.size(), 1);
    EXPECT_EQ(faces[0].adjacent_cell_ids[0], 0);
    
    // Optionally, verify node IDs if nodes are loaded
    // In this test, nodes are not loaded, so node IDs may not correspond to actual nodes
}