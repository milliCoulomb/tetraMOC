
#include "OutputHandler.hpp"
#include <fstream>

void OutputHandler::writeScalarFlux(const std::string& filepath, const std::vector<std::vector<double>>& flux) {
    std::ofstream flux_file(filepath);
    if (!flux_file.is_open()) {
        throw std::runtime_error("Unable to open flux output file.");
    }

    for (const auto& row : flux) {
        for (const auto& value : row) {
            flux_file << value << " ";
        }
        flux_file << "\n";
    }

    flux_file.close();
}

void OutputHandler::writeKEff(const std::string& filepath, double keff) {
    std::ofstream keff_file(filepath);
    if (!keff_file.is_open()) {
        throw std::runtime_error("Unable to open k_eff output file.");
    }

    keff_file << keff << "\n";
    keff_file.close();
}