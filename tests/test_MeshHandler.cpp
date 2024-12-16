// tests/test_MeshHandler.cpp
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio> // For std::remove
#include "MeshHandler.hpp"

// Utility function to create a temporary file with given content
bool createTempFile(const std::string& filename, const std::string& content) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) return false;
    outfile << content;
    outfile.close();
    return true;
}

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
    // Define content for nodes.txt
    std::string nodes_content = "4\n"
                                 "0.0 0.0 0.0\n"
                                 "1.0 0.0 0.0\n"
                                 "0.0 1.0 0.0\n"
                                 "0.0 0.0 1.0\n";
    
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
    // Define content for cells.txt
    std::string cells_content = "2\n"
                                 "0 1 2 3\n"
                                 "1 2 3 4\n"; // Assuming node IDs up to 4
    
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
    // Define content for faces.txt
    std::string faces_content = "3\n" // Number of faces
                                   "0 1 2 0 1\n" // Face with nodes 0,1,2 adjacent to cells 0 and 1
                                   "0 1 3 0\n"   // Face with nodes 0,1,3 adjacent to cell 0
                                   "1 2 3 1\n";  // Face with nodes 1,2,3 adjacent to cell 1
    
    // Create temporary faces.txt
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler failed to load faces.txt";

    // Verify loaded faces
    const auto& faces = mesh.getFaces();
    ASSERT_EQ(faces.size(), 3) << "Number of loaded faces mismatch";

    // Face 0
    EXPECT_EQ(std::get<0>(faces[0]), 0);
    EXPECT_EQ(std::get<1>(faces[0]), 1);
    EXPECT_EQ(std::get<2>(faces[0]), 2);
    EXPECT_EQ(std::get<3>(faces[0]).size(), 2);
    EXPECT_EQ(std::get<3>(faces[0])[0], 0);
    EXPECT_EQ(std::get<3>(faces[0])[1], 1);

    // Face 1
    EXPECT_EQ(std::get<0>(faces[1]), 0);
    EXPECT_EQ(std::get<1>(faces[1]), 1);
    EXPECT_EQ(std::get<2>(faces[1]), 3);
    EXPECT_EQ(std::get<3>(faces[1]).size(), 1);
    EXPECT_EQ(std::get<3>(faces[1])[0], 0);

    // Face 2
    EXPECT_EQ(std::get<0>(faces[2]), 1);
    EXPECT_EQ(std::get<1>(faces[2]), 2);
    EXPECT_EQ(std::get<2>(faces[2]), 3);
    EXPECT_EQ(std::get<3>(faces[2]).size(), 1);
    EXPECT_EQ(std::get<3>(faces[2])[0], 1);
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

TEST_F(MeshHandlerTest, LoadField_Success) {
    // Define content for nodes.txt
    std::string nodes_content = "8\n" // Number of nodes
                                 "0.0 0.0 0.0\n" // Node 0
                                 "1.0 0.0 0.0\n" // Node 1
                                 "0.0 1.0 0.0\n" // Node 2
                                 "0.0 0.0 1.0\n" // Node 3
                                 "1.0 1.0 1.0\n" // Node 4
                                 "2.0 0.0 0.0\n" // Node 5
                                 "0.0 2.0 0.0\n" // Node 6
                                 "0.0 0.0 2.0\n"; // Node 7
    
    // Define content for cells.txt
    std::string cells_content = "2\n" // Number of cells
                                 "0 1 2 3\n" // Cell 0
                                 "4 5 6 7\n"; // Cell 1
    
    // Define content for field.txt
    std::string field_content = "2\n" // Number of velocity vectors
                                   "1.0 0.0 0.0\n" // Velocity for cell 0
                                   "0.0 1.0 0.0\n"; // Velocity for cell 1
    
    // Create temporary nodes.txt
    ASSERT_TRUE(createTempFile(nodes_file, nodes_content)) << "Failed to create temporary nodes.txt";
    
    // Create temporary cells.txt
    ASSERT_TRUE(createTempFile(cells_file, cells_content)) << "Failed to create temporary cells.txt";
    
    // Create temporary field.txt
    ASSERT_TRUE(createTempFile(field_file, field_content)) << "Failed to create temporary field.txt";
    
    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadNodes(nodes_file)) << "MeshHandler failed to load nodes.txt";
    EXPECT_TRUE(mesh.loadCells(cells_file)) << "MeshHandler failed to load cells.txt";
    // Faces are not required for loading the velocity field, so we can skip loading faces here
    EXPECT_TRUE(mesh.loadVelocityField(field_file)) << "MeshHandler failed to load field.txt";

    // Verify loaded velocities
    const auto& velocities = mesh.getVelocities();
    ASSERT_EQ(velocities.size(), 2) << "Number of loaded velocities mismatch";

    EXPECT_DOUBLE_EQ(velocities[0].vx, 1.0);
    EXPECT_DOUBLE_EQ(velocities[0].vy, 0.0);
    EXPECT_DOUBLE_EQ(velocities[0].vz, 0.0);

    EXPECT_DOUBLE_EQ(velocities[1].vx, 0.0);
    EXPECT_DOUBLE_EQ(velocities[1].vy, 1.0);
    EXPECT_DOUBLE_EQ(velocities[1].vz, 0.0);
}

TEST_F(MeshHandlerTest, LoadField_MalformedFile) {
    // Define malformed content for field.txt (missing a velocity component)
    std::string field_content = "1\n"
                                 "1.0 0.0\n"; // Incomplete velocity vector
    
    // Create temporary field.txt
    ASSERT_TRUE(createTempFile(field_file, field_content)) << "Failed to create temporary field.txt";

    MeshHandler mesh;
    EXPECT_FALSE(mesh.loadVelocityField(field_file)) << "MeshHandler should fail to load malformed field.txt";
}
TEST_F(MeshHandlerTest, LoadField_IncorrectVectorCount) {
    // Define content where the number of velocity vectors specified doesn't match the actual data
    std::string field_content = "3\n" // Specifying 3 velocity vectors
                                   "1.0 0.0 0.0\n"
                                   "0.0 1.0 0.0\n"; // Only 2 vectors provided
    
    // Create temporary field.txt
    ASSERT_TRUE(createTempFile(field_file, field_content)) << "Failed to create temporary field.txt";

    MeshHandler mesh;
    EXPECT_FALSE(mesh.loadVelocityField(field_file)) << "MeshHandler should fail due to velocity vector count mismatch";
}
TEST_F(MeshHandlerTest, CombinedLoading_Success) {
    // Define content for nodes.txt
    std::string nodes_content = "4\n"
                                 "0.0 0.0 0.0\n"
                                 "1.0 0.0 0.0\n"
                                 "0.0 1.0 0.0\n"
                                 "0.0 0.0 1.0\n";
    
    // Define content for cells.txt
    std::string cells_content = "1\n"
                                 "0 1 2 3\n";
    
    // Define content for faces.txt
    std::string faces_content = "4\n"
                                 "0 1 2 0\n"
                                 "0 1 3 0\n"
                                 "0 2 3 0\n"
                                 "1 2 3 0\n";
    
    // Define content for field.txt
    std::string field_content = "1\n"
                                 "1.0 0.0 0.0\n";
    
    // Create temporary files
    ASSERT_TRUE(createTempFile(nodes_file, nodes_content)) << "Failed to create temporary nodes.txt";
    ASSERT_TRUE(createTempFile(cells_file, cells_content)) << "Failed to create temporary cells.txt";
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";
    ASSERT_TRUE(createTempFile(field_file, field_content)) << "Failed to create temporary field.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadNodes(nodes_file)) << "MeshHandler failed to load nodes.txt";
    EXPECT_TRUE(mesh.loadCells(cells_file)) << "MeshHandler failed to load cells.txt";
    EXPECT_TRUE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler failed to load faces.txt";
    EXPECT_TRUE(mesh.loadVelocityField(field_file)) << "MeshHandler failed to load field.txt";

    // Verify nodes
    const auto& nodes = mesh.getNodes();
    ASSERT_EQ(nodes.size(), 4);
    EXPECT_DOUBLE_EQ(nodes[0].x, 0.0);
    EXPECT_DOUBLE_EQ(nodes[3].z, 1.0);

    // Verify cells
    const auto& cells = mesh.getCells();
    ASSERT_EQ(cells.size(), 1);
    EXPECT_EQ(cells[0].node_ids[0], 0);
    EXPECT_EQ(cells[0].cell_id, 0);

    // Verify faces
    const auto& faces = mesh.getFaces();
    ASSERT_EQ(faces.size(), 4);
    EXPECT_EQ(std::get<3>(faces[0]).size(), 1); // Each face has one adjacent cell

    // Verify velocities
    const auto& velocities = mesh.getVelocities();
    ASSERT_EQ(velocities.size(), 1);
    EXPECT_DOUBLE_EQ(velocities[0].vx, 1.0);
}
TEST_F(MeshHandlerTest, Getters_ReturnCorrectData) {
    // Define content for nodes.txt
    std::string nodes_content = "2\n"
                                 "0.0 0.0 0.0\n"
                                 "1.0 1.0 1.0\n";
    
    // Define content for cells.txt
    std::string cells_content = "1\n"
                                 "0 1 2 3\n"; // Assuming node IDs up to 3
    
    // Define content for faces.txt
    std::string faces_content = "2\n"
                                 "0 1 2 0\n"
                                 "1 2 3 0\n";
    
    // Define content for field.txt
    std::string field_content = "1\n"
                                 "0.5 0.5 0.5\n";
    
    // Create temporary files
    ASSERT_TRUE(createTempFile(nodes_file, nodes_content)) << "Failed to create temporary nodes.txt";
    ASSERT_TRUE(createTempFile(cells_file, cells_content)) << "Failed to create temporary cells.txt";
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";
    ASSERT_TRUE(createTempFile(field_file, field_content)) << "Failed to create temporary field.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadNodes(nodes_file));
    EXPECT_TRUE(mesh.loadCells(cells_file));
    EXPECT_TRUE(mesh.loadFaceConnectivity(faces_file));
    EXPECT_TRUE(mesh.loadVelocityField(field_file));

    // Test getNodes
    const auto& nodes = mesh.getNodes();
    ASSERT_EQ(nodes.size(), 2);
    EXPECT_DOUBLE_EQ(nodes[0].x, 0.0);
    EXPECT_DOUBLE_EQ(nodes[1].y, 1.0);

    // Test getCells
    const auto& cells = mesh.getCells();
    ASSERT_EQ(cells.size(), 1);
    EXPECT_EQ(cells[0].cell_id, 0);
    EXPECT_EQ(cells[0].node_ids[3], 3);

    // Test getFaces
    const auto& faces = mesh.getFaces();
    ASSERT_EQ(faces.size(), 2);
    EXPECT_EQ(std::get<0>(faces[1]), 1);
    EXPECT_EQ(std::get<3>(faces[1]).size(), 1);
    EXPECT_EQ(std::get<3>(faces[1])[0], 0);

    // Test getVelocities
    const auto& velocities = mesh.getVelocities();
    ASSERT_EQ(velocities.size(), 1);
    EXPECT_DOUBLE_EQ(velocities[0].vz, 0.5);
}
TEST_F(MeshHandlerTest, LoadFaces_BeforeNodesAndCells) {
    // Define content for faces.txt referencing node and cell IDs that haven't been loaded yet
    std::string faces_content = "1\n"
                                 "0 1 2 0\n";
    
    // Create temporary faces.txt
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler should load faces even if nodes and cells are not loaded yet";

    // However, accessing face data should still work, but might reference non-existent nodes/cells
    const auto& faces = mesh.getFaces();
    ASSERT_EQ(faces.size(), 1);
    EXPECT_EQ(std::get<3>(faces[0]).size(), 1);
    EXPECT_EQ(std::get<3>(faces[0])[0], 0);
}
TEST_F(MeshHandlerTest, LoadLargeMesh) {
    // Define a larger number of nodes, cells, faces, and velocity vectors
    int num_nodes = 1000;
    int num_cells = 2000;
    int num_faces = 4000;
    
    // Generate nodes.txt content
    std::string nodes_content = std::to_string(num_nodes) + "\n";
    for(int i = 0; i < num_nodes; ++i) {
        nodes_content += std::to_string(i*0.1) + " " + std::to_string(i*0.1) + " " + std::to_string(i*0.1) + "\n";
    }
    
    // Generate cells.txt content
    std::string cells_content = std::to_string(num_cells) + "\n";
    for(int i = 0; i < num_cells; ++i) {
        // Simple pattern for node IDs
        cells_content += std::to_string(i % num_nodes) + " " + 
                         std::to_string((i+1) % num_nodes) + " " +
                         std::to_string((i+2) % num_nodes) + " " +
                         std::to_string((i+3) % num_nodes) + "\n";
    }
    
    // Generate faces.txt content
    std::string faces_content = std::to_string(num_faces) + "\n";
    for(int i = 0; i < num_faces; ++i) {
        // Simple pattern for face node IDs and adjacent cells
        faces_content += std::to_string(i % num_nodes) + " " + 
                         std::to_string((i+1) % num_nodes) + " " +
                         std::to_string((i+2) % num_nodes) + " " +
                         std::to_string(i % num_cells) + " " +
                         std::to_string((i+1) % num_cells) + "\n";
    }
    
    // Generate field.txt content
    std::string field_content = std::to_string(num_cells) + "\n";
    for(int i = 0; i < num_cells; ++i) {
        field_content += "1.0 0.0 0.0\n"; // Uniform velocity
    }
    
    // Create temporary files
    ASSERT_TRUE(createTempFile(nodes_file, nodes_content)) << "Failed to create temporary nodes.txt";
    ASSERT_TRUE(createTempFile(cells_file, cells_content)) << "Failed to create temporary cells.txt";
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";
    ASSERT_TRUE(createTempFile(field_file, field_content)) << "Failed to create temporary field.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadNodes(nodes_file)) << "MeshHandler failed to load nodes.txt";
    EXPECT_TRUE(mesh.loadCells(cells_file)) << "MeshHandler failed to load cells.txt";
    EXPECT_TRUE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler failed to load faces.txt";
    EXPECT_TRUE(mesh.loadVelocityField(field_file)) << "MeshHandler failed to load field.txt";

    // Verify counts
    EXPECT_EQ(mesh.getNodes().size(), num_nodes);
    EXPECT_EQ(mesh.getCells().size(), num_cells);
    EXPECT_EQ(mesh.getFaces().size(), num_faces);
    EXPECT_EQ(mesh.getVelocities().size(), num_cells);
}
TEST_F(MeshHandlerTest, LoadField_VelocityCountMismatch) {
    // Define content for nodes.txt
    std::string nodes_content = "4\n"
                                 "0.0 0.0 0.0\n"
                                 "1.0 0.0 0.0\n"
                                 "0.0 1.0 0.0\n"
                                 "0.0 0.0 1.0\n";
    
    // Define content for cells.txt
    std::string cells_content = "1\n"
                                 "0 1 2 3\n";
    
    // Define content for faces.txt
    std::string faces_content = "4\n"
                                 "0 1 2 0\n"
                                 "0 1 3 0\n"
                                 "0 2 3 0\n"
                                 "1 2 3 0\n";
    
    // Define content for field.txt with mismatched velocity count
    std::string field_content = "2\n" // Should be 1
                                 "1.0 0.0 0.0\n"
                                 "0.0 1.0 0.0\n";
    
    // Create temporary files
    ASSERT_TRUE(createTempFile(nodes_file, nodes_content)) << "Failed to create temporary nodes.txt";
    ASSERT_TRUE(createTempFile(cells_file, cells_content)) << "Failed to create temporary cells.txt";
    ASSERT_TRUE(createTempFile(faces_file, faces_content)) << "Failed to create temporary faces.txt";
    ASSERT_TRUE(createTempFile(field_file, field_content)) << "Failed to create temporary field.txt";

    MeshHandler mesh;
    EXPECT_TRUE(mesh.loadNodes(nodes_file)) << "MeshHandler failed to load nodes.txt";
    EXPECT_TRUE(mesh.loadCells(cells_file)) << "MeshHandler failed to load cells.txt";
    EXPECT_TRUE(mesh.loadFaceConnectivity(faces_file)) << "MeshHandler failed to load faces.txt";
    EXPECT_FALSE(mesh.loadVelocityField(field_file)) << "MeshHandler should fail due to velocity vector count mismatch";
}