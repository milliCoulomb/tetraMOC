#include "Logger.hpp"
#include <gtest/gtest.h>
#include <sstream>

// Helper function to capture std::cout
std::string captureStdOut(const std::function<void()>& func) {
    std::streambuf* originalCout = std::cout.rdbuf();
    std::ostringstream ss;
    std::cout.rdbuf(ss.rdbuf());
    
    func();
    
    std::cout.rdbuf(originalCout);
    return ss.str();
}

TEST(LoggerTest, LogsInfoWhenLevelIsInfo) {
    Logger::setLogLevel(LogLevel::INFO);
    std::string output = captureStdOut([&]() {
        Logger::info("This is an info message.");
    });
    EXPECT_NE(output.find("[INFO] This is an info message."), std::string::npos);
}

TEST(LoggerTest, DoesNotLogInfoWhenLevelIsWarning) {
    Logger::setLogLevel(LogLevel::WARNING);
    std::string output = captureStdOut([&]() {
        Logger::info("This is an info message.");
    });
    EXPECT_TRUE(output.empty());
}