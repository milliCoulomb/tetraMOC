#include "Logger.hpp"

LogLevel Logger::currentLevel = LogLevel::INFO;

void Logger::setLogLevel(LogLevel level) {
    currentLevel = level;
}

#ifndef DISABLE_LOGGING

void Logger::log(const std::string& message, LogLevel level) {
    if (level <= currentLevel) {
        switch(level) {
            case LogLevel::RUNNING:
                std::cout << "[RUNNING] " << message << std::endl;
                break;
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

void Logger::running(const std::string& message) {
    log(message, LogLevel::RUNNING);
}

#else // DISABLE_LOGGING is defined

void Logger::log(const std::string&, LogLevel) {
    // Logging is disabled; do nothing
}

void Logger::info(const std::string&) {
    // Logging is disabled; do nothing
}

void Logger::warning(const std::string&) {
    // Logging is disabled; do nothing
}

void Logger::error(const std::string&) {
    // Logging is disabled; do nothing
}

void Logger::running(const std::string&) {
    // Logging is disabled; do nothing
}

#endif // DISABLE_LOGGING