// tests/TestUtils.hpp
#ifndef TEST_UTILS_HPP
#define TEST_UTILS_HPP

#include <string>

inline bool createTempFile(const std::string& filename, const std::string& content) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) return false;
    outfile << content;
    outfile.close();
    return true;
}

#endif // TEST_UTILS_HPP