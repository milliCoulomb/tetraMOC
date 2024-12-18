// tests/test_Field.cpp

#include <gtest/gtest.h>
#include "Field.hpp"
#include <fstream>
#include <cstdio> // For std::remove
#include "TestUtils.hpp" // Assumed to contain createTempFile and vectorsAlmostEqual

// Test Fixture for Field
class FieldTest : public ::testing::Test {
protected:
    // Temporary file names
    std::string vector_field_file = "temp_vector_field.txt";
    std::string scalar_field_file = "temp_scalar_field.txt";

    // Clean up temporary files after each test
    void TearDown() override {
        std::remove(vector_field_file.c_str());
        std::remove(scalar_field_file.c_str());
    }
};

// Test case: Loading a valid vector field
TEST_F(FieldTest, LoadVectorField_Success) {
    // Define content for vector_field.txt
    std::string field_content = "2\n" // Number of vectors
                                   "1.0 0.0 0.0\n" // Vector for cell 0
                                   "0.0 1.0 0.0\n"; // Vector for cell 1

    // Create temporary vector_field.txt
    ASSERT_TRUE(createTempFile(vector_field_file, field_content)) << "Failed to create temporary vector_field.txt";

    Field field;
    EXPECT_TRUE(field.loadVectorField(vector_field_file)) << "Field failed to load vector_field.txt";

    // Verify loaded vectors
    const auto& vectors = field.getVectorFields();
    ASSERT_EQ(vectors.size(), 2) << "Number of loaded vectors mismatch";

    EXPECT_TRUE(vectors[0].isAlmostEqual(Vector3D(1.0, 0.0, 0.0))) << "Vector for cell 0 mismatch";
    EXPECT_TRUE(vectors[1].isAlmostEqual(Vector3D(0.0, 1.0, 0.0))) << "Vector for cell 1 mismatch";
}

// Test case: Loading a malformed vector field (missing a component)
TEST_F(FieldTest, LoadVectorField_MalformedFile) {
    // Define malformed content for vector_field.txt (missing vz for cell 1)
    std::string field_content = "2\n" // Number of vectors
                                   "1.0 0.0 0.0\n" // Vector for cell 0
                                   "0.0 1.0\n"; // Incomplete vector for cell 1

    // Create temporary vector_field.txt
    ASSERT_TRUE(createTempFile(vector_field_file, field_content)) << "Failed to create temporary vector_field.txt";

    Field field;
    EXPECT_FALSE(field.loadVectorField(vector_field_file)) << "Field should fail to load malformed vector_field.txt";
}

// Test case: Loading a vector field with incorrect number of vectors
TEST_F(FieldTest, LoadVectorField_VectorCountMismatch) {
    // Define content where the number of vectors specified doesn't match the actual data
    std::string field_content = "3\n" // Specifying 3 vectors
                                   "1.0 0.0 0.0\n" // Vector for cell 0
                                   "0.0 1.0 0.0\n"; // Vector for cell 1

    // Create temporary vector_field.txt
    ASSERT_TRUE(createTempFile(vector_field_file, field_content)) << "Failed to create temporary vector_field.txt";

    Field field;
    EXPECT_FALSE(field.loadVectorField(vector_field_file)) << "Field should fail due to vector count mismatch";
}

// Test case: Loading a valid scalar field (for future extension)
TEST_F(FieldTest, LoadScalarField_Success) {
    // Define content for scalar_field.txt
    std::string field_content = "2\n" // Number of scalar values
                                   "100.0\n" // Scalar for cell 0
                                   "200.0\n"; // Scalar for cell 1

    // Create temporary scalar_field.txt
    ASSERT_TRUE(createTempFile(scalar_field_file, field_content)) << "Failed to create temporary scalar_field.txt";

    Field field;
    EXPECT_TRUE(field.loadScalarField(scalar_field_file)) << "Field failed to load scalar_field.txt";

    // Verify loaded scalars
    const auto& scalars = field.getScalarFields();
    ASSERT_EQ(scalars.size(), 2) << "Number of loaded scalars mismatch";

    EXPECT_DOUBLE_EQ(scalars[0].value, 100.0) << "Scalar value for cell 0 mismatch";
    EXPECT_DOUBLE_EQ(scalars[1].value, 200.0) << "Scalar value for cell 1 mismatch";
}

// Test case: Loading a malformed scalar field (extra data)
TEST_F(FieldTest, LoadScalarField_MalformedFile) {
    // Define malformed content for scalar_field.txt (extra data for cell 1)
    std::string field_content = "2\n" // Number of scalar values
                                   "100.0\n" // Scalar for cell 0
                                   "200.0 300.0\n"; // Extra data for cell 1

    // Create temporary scalar_field.txt
    ASSERT_TRUE(createTempFile(scalar_field_file, field_content)) << "Failed to create temporary scalar_field.txt";

    Field field;
    EXPECT_FALSE(field.loadScalarField(scalar_field_file)) << "Field should fail to load malformed scalar_field.txt";
}

// Test case: Loading a scalar field with incorrect number of values
TEST_F(FieldTest, LoadScalarField_ScalarCountMismatch) {
    // Define content where the number of scalar values specified doesn't match the actual data
    std::string field_content = "3\n" // Specifying 3 scalars
                                   "100.0\n" // Scalar for cell 0
                                   "200.0\n"; // Scalar for cell 1

    // Create temporary scalar_field.txt
    ASSERT_TRUE(createTempFile(scalar_field_file, field_content)) << "Failed to create temporary scalar_field.txt";

    Field field;
    EXPECT_FALSE(field.loadScalarField(scalar_field_file)) << "Field should fail due to scalar count mismatch";
}

// Additional Test Case: Setting and Getting Direction
TEST_F(FieldTest, SetAndGetDirection) {
    Field field;
    Vector3D direction(1.0, 1.0, 1.0);
    field.setDirection(direction);
    EXPECT_TRUE(field.getDirection().isAlmostEqual(direction)) << "Direction vector mismatch after setting";
}