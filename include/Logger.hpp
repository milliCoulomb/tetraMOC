// Utilities/Logger.hpp
#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>
#include <iostream>

enum class LogLevel { INFO, WARNING, ERROR };

class Logger {
public:
    static void log(const std::string& message, LogLevel level = LogLevel::INFO);
    static void info(const std::string& message);
    static void warning(const std::string& message);
    static void error(const std::string& message);
};

#endif // LOGGER_HPP