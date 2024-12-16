// tests/TestUtils.hpp
#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>

namespace SNSolver {

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

} // namespace SNSolver

#endif // TEST_UTILS_HPP