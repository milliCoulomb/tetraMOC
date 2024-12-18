// tests/TestUtils.hpp
#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "Vector3D.hpp"

// Utility function to create a temporary file with given content
inline bool createTempFile(const std::string& filename, const std::string& content) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) return false;
    outfile << content;
    outfile.close();
    return true;
}

// Utility function to generate a unique temporary filename
inline std::string generateUniqueFilename(const std::string& prefix = "temp_", const std::string& suffix = ".txt") {
    // Generate a unique filename using rand and time
    std::string filename = prefix;
    filename += std::to_string(std::rand());
    filename += suffix;
    return filename;
}

// Helper function to compare two Vector3D instances with a tolerance
inline bool vectorsAlmostEqual(const Vector3D& v1, const Vector3D& v2, double tol = 1e-6) {
    return v1.isAlmostEqual(v2, tol);
}

#endif // TEST_UTILS_HPP