#include "OutputHandler.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>

// Helper function to read a file's contents into a string
std::string readFile(const std::string& filepath) {
    std::ifstream file(filepath);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

TEST(OutputHandlerTest, WriteScalarFluxSuccessfully) {
    OutputHandler output_handler;
    std::string filepath = "/tmp/test_scalar_flux.txt";
    std::vector<std::vector<double>> flux = {
        {1.1, 2.2, 3.3},
        {4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9}
    };

    // Write scalar flux to file
    EXPECT_NO_THROW(output_handler.writeScalarFlux(filepath, flux));

    // Read the file and verify contents
    std::string expected_content = "1.1 2.2 3.3 \n4.4 5.5 6.6 \n7.7 8.8 9.9 \n";
    std::string actual_content = readFile(filepath);
    EXPECT_EQ(actual_content, expected_content);

    // Clean up
    std::filesystem::remove(filepath);
}

TEST(OutputHandlerTest, WriteKEffSuccessfully) {
    OutputHandler output_handler;
    std::string filepath = "/tmp/test_keff.txt";
    double keff = 1.025;

    // Write KEff to file
    EXPECT_NO_THROW(output_handler.writeKEff(filepath, keff));

    // Read the file and verify contents
    std::string expected_content = "1.025\n";
    std::string actual_content = readFile(filepath);
    EXPECT_EQ(actual_content, expected_content);

    // Clean up
    std::filesystem::remove(filepath);
}

TEST(OutputHandlerTest, WriteScalarFluxThrowsOnInvalidPath) {
    OutputHandler output_handler;
    // Attempt to write to a directory path instead of a file
    std::string invalid_filepath = "/invalid_path/test_scalar_flux.txt";
    std::vector<std::vector<double>> flux = {
        {1.1, 2.2, 3.3},
        {4.4, 5.5, 6.6}
    };

    EXPECT_THROW(output_handler.writeScalarFlux(invalid_filepath, flux), std::runtime_error);
}

TEST(OutputHandlerTest, WriteKEffThrowsOnInvalidPath) {
    OutputHandler output_handler;
    // Attempt to write to an invalid file path
    std::string invalid_filepath = "/invalid_path/test_keff.txt";
    double keff = 1.025;

    EXPECT_THROW(output_handler.writeKEff(invalid_filepath, keff), std::runtime_error);
}

TEST(OutputHandlerTest, WriteScalarFluxHandlesEmptyFlux) {
    OutputHandler output_handler;
    std::string filepath = "/tmp/test_empty_flux.txt";
    std::vector<std::vector<double>> flux = {};

    // Write empty scalar flux to file
    EXPECT_NO_THROW(output_handler.writeScalarFlux(filepath, flux));

    // Read the file and verify it's empty
    std::string actual_content = readFile(filepath);
    EXPECT_TRUE(actual_content.empty());

    // Clean up
    std::filesystem::remove(filepath);
}

TEST(OutputHandlerTest, WriteKEffHandlesZeroValue) {
    OutputHandler output_handler;
    std::string filepath = "/tmp/test_zero_keff.txt";
    double keff = 0.0;

    // Write KEff with zero value to file
    EXPECT_NO_THROW(output_handler.writeKEff(filepath, keff));

    // Read the file and verify contents
    std::string expected_content = "0\n";
    std::string actual_content = readFile(filepath);
    EXPECT_EQ(actual_content, expected_content);

    // Clean up
    std::filesystem::remove(filepath);
}