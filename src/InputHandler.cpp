// Utilities/InputHandler.cpp
#include "InputHandler.hpp"
#include "Logger.hpp" // Assuming Logger is used for error messages
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

bool InputHandler::loadData(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        Logger::error("Failed to open data file: " + filename);
        return false;
    }

    std::string line;
    
    // Read the first line to get the number of groups
    if (!std::getline(infile, line)) {
        Logger::error("Data file is empty: " + filename);
        return false;
    }

    std::istringstream iss_num(line);
    if (!(iss_num >> num_groups_) || num_groups_ <= 0) {
        Logger::error("Invalid number of energy groups in file: " + filename);
        return false;
    }

    // Clear existing data
    energy_groups_.clear();
    energy_groups_.reserve(num_groups_);

    // Read data for each group
    int current_group = 0;
    while (std::getline(infile, line)) {
        if (line.empty()) continue; // Skip empty lines

        std::istringstream iss(line);
        EnergyGroupData group_data;
        if (!(iss >> group_data.total_xs >> group_data.fission_xs >> group_data.scattering_xs
                  >> group_data.multiplicity >> group_data.fission_spectrum >> group_data.delayed_spectrum)) {
            Logger::error("Malformed data at group " + std::to_string(current_group) + " in file: " + filename);
            return false;
        }

        energy_groups_.push_back(group_data);
        current_group++;

        if (current_group > num_groups_) {
            Logger::warning("Extra data found beyond the specified number of groups in file: " + filename);
            break;
        }
    }

    if (current_group < num_groups_) {
        Logger::error("Insufficient data in file: " + filename + ". Expected " + std::to_string(num_groups_) + " groups, found " + std::to_string(current_group));
        return false;
    }

    // Optionally, handle any remaining lines as warnings
    if (current_group < num_groups_) {
        Logger::warning("Number of groups specified (" + std::to_string(num_groups_) + ") does not match the number of data lines (" + std::to_string(current_group) + ") in file: " + filename);
    }

    Logger::info("Successfully loaded data for " + std::to_string(num_groups_) + " energy groups from file: " + filename);
    return true;
}

int InputHandler::getNumGroups() const {
    return num_groups_;
}

InputHandler::EnergyGroupData InputHandler::getEnergyGroupData(int group) const {
    if (group < 0 || group >= num_groups_) {
        Logger::error("Requested data for invalid group index: " + std::to_string(group));
        throw std::out_of_range("Invalid energy group index.");
    }
    return energy_groups_[group];
}