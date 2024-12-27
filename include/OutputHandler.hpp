#ifndef OUTPUTHANDLER_HPP
#define OUTPUTHANDLER_HPP

#include <string>
#include <vector>

class OutputHandler {
public:
    void writeScalarFlux(const std::string& filepath, const std::vector<std::vector<double>>& flux);
    void writeKEff(const std::string& filepath, double keff);
};

#endif // OUTPUTHANDLER_HPP
