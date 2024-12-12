// include/field.h
#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <string>

struct CellField {
    double vx, vy, vz;
};

class Field {
public:
    Field() = default;
    bool loadField(const std::string& filename);

    const std::vector<CellField>& getCellValues() const { return cellValues; }

private:
    std::vector<CellField> cellValues;
};

#endif // FIELD_H