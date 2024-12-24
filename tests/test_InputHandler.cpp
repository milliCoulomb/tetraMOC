// tests/test_InputHandler.cpp

#include <gtest/gtest.h>
#include "InputHandler.hpp"
#include "TestUtils.hpp"

#include <cstdio> // For std::remove

// Test Fixture for InputHandler
class InputHandlerTest : public ::testing::Test {
protected:
    // Temporary file name
    std::string temp_file = "temp_data.txt";

    // Clean up temporary files after each test
    void TearDown() override {
        std::remove(temp_file.c_str());
    }
};

// Test Case 1: Valid Input with Multiple Energy Groups
TEST_F(InputHandlerTest, ValidInputMultipleGroups) {
    std::string content = "3\n"
                          "1.0 0.1 0.2 0.1 0.3 2.43 0.98 0.02\n"
                          "1.2 0.0 0.6 0.1 0.3 2.50 0.95 0.05\n"
                          "0.8 0.05 0.1 0.1 0.1 2.35 0.97 0.03\n";

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_TRUE(input_handler.loadData(temp_file)) << "Failed to load valid data.";

    ASSERT_EQ(input_handler.getNumGroups(), 3) << "Number of groups should be 3.";

    // Verify Group 0
    InputHandler::EnergyGroupData group0 = input_handler.getEnergyGroupData(0);
    EXPECT_DOUBLE_EQ(group0.total_xs, 1.0);
    EXPECT_DOUBLE_EQ(group0.fission_xs, 0.1);
    EXPECT_DOUBLE_EQ(group0.scattering_xs[0], 0.2);
    EXPECT_DOUBLE_EQ(group0.scattering_xs[1], 0.1);
    EXPECT_DOUBLE_EQ(group0.scattering_xs[2], 0.3);
    EXPECT_DOUBLE_EQ(group0.multiplicity, 2.43);
    EXPECT_DOUBLE_EQ(group0.fission_spectrum, 0.98);
    EXPECT_DOUBLE_EQ(group0.delayed_spectrum, 0.02);

    // Verify Group 1
    InputHandler::EnergyGroupData group1 = input_handler.getEnergyGroupData(1);
    EXPECT_DOUBLE_EQ(group1.total_xs, 1.2);
    EXPECT_DOUBLE_EQ(group1.fission_xs, 0.0);
    EXPECT_DOUBLE_EQ(group1.scattering_xs[0], 0.6);
    EXPECT_DOUBLE_EQ(group1.scattering_xs[1], 0.1);
    EXPECT_DOUBLE_EQ(group1.scattering_xs[2], 0.3);
    EXPECT_DOUBLE_EQ(group1.multiplicity, 2.50);
    EXPECT_DOUBLE_EQ(group1.fission_spectrum, 0.95);
    EXPECT_DOUBLE_EQ(group1.delayed_spectrum, 0.05);

    // Verify Group 2
    InputHandler::EnergyGroupData group2 = input_handler.getEnergyGroupData(2);
    EXPECT_DOUBLE_EQ(group2.total_xs, 0.8);
    EXPECT_DOUBLE_EQ(group2.fission_xs, 0.05);
    EXPECT_DOUBLE_EQ(group2.scattering_xs[0], 0.1);
    EXPECT_DOUBLE_EQ(group2.scattering_xs[1], 0.1);
    EXPECT_DOUBLE_EQ(group2.scattering_xs[2], 0.1);
    EXPECT_DOUBLE_EQ(group2.multiplicity, 2.35);
    EXPECT_DOUBLE_EQ(group2.fission_spectrum, 0.97);
    EXPECT_DOUBLE_EQ(group2.delayed_spectrum, 0.03);
}

// Test Case 2: Input File with Insufficient Data
TEST_F(InputHandlerTest, InsufficientData) {
    std::string content = "2\n"
                          "1.0 0.1 0.9 2.43 0.98 0.02\n"; // Only one group's data provided

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_FALSE(input_handler.loadData(temp_file)) << "Should fail due to insufficient data.";
}

// Test Case 3: Input File with Extra Data
TEST_F(InputHandlerTest, ExtraData) {
    std::string content = "2\n"
                          "1.0 0.1 0.5 0.1 2.43 0.98 0.02\n"
                          "1.2 0.0 0.1 0.9 2.50 0.95 0.05\n"
                          "0.8 0.05 0.0 0.75 2.35 0.97 0.03\n"; // Extra group's data

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    // Depending on implementation, this might pass with a warning
    // Since our implementation logs a warning but does not fail, we expect it to pass
    ASSERT_TRUE(input_handler.loadData(temp_file)) << "Should load data despite extra entries.";

    ASSERT_EQ(input_handler.getNumGroups(), 2) << "Number of groups should be 2.";

    // Verify Group 0
    InputHandler::EnergyGroupData group0 = input_handler.getEnergyGroupData(0);
    EXPECT_DOUBLE_EQ(group0.total_xs, 1.0);
    EXPECT_DOUBLE_EQ(group0.fission_xs, 0.1);
    EXPECT_DOUBLE_EQ(group0.scattering_xs[0], 0.5);
    EXPECT_DOUBLE_EQ(group0.scattering_xs[1], 0.1);
    EXPECT_DOUBLE_EQ(group0.multiplicity, 2.43);
    EXPECT_DOUBLE_EQ(group0.fission_spectrum, 0.98);
    EXPECT_DOUBLE_EQ(group0.delayed_spectrum, 0.02);

    // Verify Group 1
    InputHandler::EnergyGroupData group1 = input_handler.getEnergyGroupData(1);
    EXPECT_DOUBLE_EQ(group1.total_xs, 1.2);
    EXPECT_DOUBLE_EQ(group1.fission_xs, 0.0);
    EXPECT_DOUBLE_EQ(group1.scattering_xs[0], 0.1);
    EXPECT_DOUBLE_EQ(group1.scattering_xs[1], 0.9);
    EXPECT_DOUBLE_EQ(group1.multiplicity, 2.50);
    EXPECT_DOUBLE_EQ(group1.fission_spectrum, 0.95);
    EXPECT_DOUBLE_EQ(group1.delayed_spectrum, 0.05);
}

// Test Case 4: Input File with Non-Numeric Values
TEST_F(InputHandlerTest, NonNumericValues) {
    std::string content = "2\n"
                          "1.0 0.1 0.7 0.1 2.43 0.98 0.02\n"
                          "abc 0.0 1.2 0.1 2.50 0.95 0.05\n"; // Non-numeric total cross section for group 2

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_FALSE(input_handler.loadData(temp_file)) << "Should fail due to non-numeric values.";
}

// Test Case 5: Input File with Zero Energy Groups
TEST_F(InputHandlerTest, ZeroEnergyGroups) {
    std::string content = "0\n"; // Zero groups

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_FALSE(input_handler.loadData(temp_file)) << "Should fail due to zero energy groups.";
}

// Test Case 6: Requesting Data for Invalid Group Indices
TEST_F(InputHandlerTest, InvalidGroupIndices) {
    std::string content = "2\n"
                          "1.0 0.1 0.7 0.1 2.43 0.98 0.02\n"
                          "1.2 0.0 1.2 0.1 2.50 0.95 0.05\n";

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_TRUE(input_handler.loadData(temp_file)) << "Failed to load valid data.";

    // Valid group indices
    EXPECT_NO_THROW({
        InputHandler::EnergyGroupData group0 = input_handler.getEnergyGroupData(0);
    });

    EXPECT_NO_THROW({
        InputHandler::EnergyGroupData group1 = input_handler.getEnergyGroupData(1);
    });

    // Invalid group indices
    EXPECT_THROW({
        input_handler.getEnergyGroupData(-1);
    }, std::out_of_range);

    EXPECT_THROW({
        input_handler.getEnergyGroupData(2);
    }, std::out_of_range);
}

// Test Case 7: Large Number of Energy Groups
TEST_F(InputHandlerTest, LargeNumberOfGroups) {
    int num_groups = 1000;
    std::ostringstream oss;
    oss << num_groups << "\n";
    for(int i = 0; i < num_groups; ++i) {
        oss << (1.0 + i * 0.01) << " "    // total_xs
            << (0.1 + i * 0.001) << " "; // fission_xs
        // scattering_xs for each group
        for(int j = 0; j < num_groups; ++j) {
            oss << (0.9 + i * 0.005 + j * 0.0001) << " ";
        }
        oss << (2.0 + i * 0.01) << " "    // multiplicity
            << (0.95 + i * 0.001) << " "  // fission_spectrum
            << (0.05 + i * 0.0001) << "\n"; // delayed_spectrum
    }
    std::string content = oss.str();

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_TRUE(input_handler.loadData(temp_file)) << "Failed to load large data.";

    ASSERT_EQ(input_handler.getNumGroups(), num_groups) << "Number of groups should match the input.";

    for(int i = 0; i < num_groups; ++i) {
        InputHandler::EnergyGroupData group = input_handler.getEnergyGroupData(i);
        EXPECT_DOUBLE_EQ(group.total_xs, 1.0 + i * 0.01);
        EXPECT_DOUBLE_EQ(group.fission_xs, 0.1 + i * 0.001);
        ASSERT_EQ(group.scattering_xs.size(), static_cast<size_t>(num_groups)) << "Scattering_xs size mismatch.";
        for(int j = 0; j < num_groups; ++j) {
            EXPECT_DOUBLE_EQ(group.scattering_xs[j], 0.9 + i * 0.005 + j * 0.0001);
        }
        EXPECT_DOUBLE_EQ(group.multiplicity, 2.0 + i * 0.01);
        EXPECT_DOUBLE_EQ(group.fission_spectrum, 0.95 + i * 0.001);
        EXPECT_DOUBLE_EQ(group.delayed_spectrum, 0.05 + i * 0.0001);
    }
}

// Test Case 8: Mixed Valid and Invalid Data Lines
TEST_F(InputHandlerTest, MixedValidAndInvalidDataLines) {
    std::string content = "3\n"
                          "1.0 0.1 0.2 0.1 0.3 2.43 0.98 0.02\n"
                          "1.2 0.0 0.6 0.1 0.3 2.50 0.95 caca\n"
                          "0.8 0.05 0.1 0.1 0.1 2.35 0.97 0.03\n";

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_FALSE(input_handler.loadData(temp_file)) << "Should fail due to malformed data.";
}

// Test Case 9: Empty Input File
TEST_F(InputHandlerTest, EmptyInputFile) {
    std::string content = ""; // Empty file

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_FALSE(input_handler.loadData(temp_file)) << "Should fail due to empty file.";
}

// Test Case 10: Input File with Only Header
TEST_F(InputHandlerTest, OnlyHeaderInputFile) {
    std::string content = "2\n"; // Only header, no data lines

    ASSERT_TRUE(createTempFile(temp_file, content)) << "Failed to create temporary input file.";

    InputHandler input_handler;
    ASSERT_FALSE(input_handler.loadData(temp_file)) << "Should fail due to missing data lines.";
}