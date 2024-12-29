// Utilities/Logger.hpp
#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>
#include <iostream>

enum class LogLevel { RUNNING, ERROR, WARNING, INFO };

class Logger {
public:
    static void setLogLevel(LogLevel level);

    /**
     * @brief Logs a message with the specified log level.
     * 
     * @param message The message to log.
     * @param level The severity level of the log.
     */
    static void log(const std::string& message, LogLevel level = LogLevel::INFO);
    
    /**
     * @brief Logs an informational message.
     * 
     * @param message The message to log.
     */
    static void info(const std::string& message);
    
    /**
     * @brief Logs a warning message.
     * 
     * @param message The message to log.
     */
    static void warning(const std::string& message);
    
    /**
     * @brief Logs an error message.
     * 
     * @param message The message to log.
     */
    static void error(const std::string& message);

    /**
     * @brief Logs a running message.
     * 
     * @param message The message to log.
     */
    static void running(const std::string& message);
    
private:
    static LogLevel currentLevel;
};

#endif // LOGGER_HPP