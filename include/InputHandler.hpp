// Utilities/InputHandler.hpp
#ifndef INPUT_HANDLER_HPP
#define INPUT_HANDLER_HPP

#include <string>
#include <vector>

class InputHandler {
public:
    /**
     * @brief Struct to hold data for each energy group.
     */
    struct EnergyGroupData {
        double total_xs;
        double fission_xs;
        double scattering_xs;
        double multiplicity;
        double fission_spectrum;
        double delayed_spectrum;
    };

    /**
     * @brief Loads cross-section and nuclear data from the specified input file.
     * 
     * @param filename Path to the input file.
     * @return true if the data was loaded successfully.
     * @return false if there was an error loading the data.
     */
    bool loadData(const std::string& filename);

    /**
     * @brief Retrieves the number of energy groups.
     * 
     * @return int Number of energy groups.
     */
    int getNumGroups() const;

    /**
     * @brief Retrieves the data for a given energy group.
     * 
     * @param group Energy group index (0-based).
     * @return EnergyGroupData Data of the specified energy group.
     * @throws std::out_of_range if the group index is invalid.
     */
    EnergyGroupData getEnergyGroupData(int group) const;

private:
    int num_groups_ = 0;
    std::vector<EnergyGroupData> energy_groups_;
};

#endif // INPUT_HANDLER_HPP