// Utilities/Logger.cpp
#include "Logger.hpp"

namespace SNSolver {

void Logger::log(const std::string& message, LogLevel level) {
    switch(level) {
        case LogLevel::INFO:
            std::cout << "[INFO] " << message << std::endl;
            break;
        case LogLevel::WARNING:
            std::cout << "[WARNING] " << message << std::endl;
            break;
        case LogLevel::ERROR:
            std::cerr << "[ERROR] " << message << std::endl;
            break;
    }
}

void Logger::info(const std::string& message) {
    log(message, LogLevel::INFO);
}

void Logger::warning(const std::string& message) {
    log(message, LogLevel::WARNING);
}

void Logger::error(const std::string& message) {
    log(message, LogLevel::ERROR);
}

}
