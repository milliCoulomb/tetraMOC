// include/RayTracer.hpp

#ifndef RAY_TRACER_HPP
#define RAY_TRACER_HPP

#include "MeshHandler.hpp"
#include "Field.hpp"
#include "Vector3D.hpp"
#include "TrackingData.hpp"
#include "Logger.hpp"

#include <vector>
#include <memory>
#include <cassert>

// Enumeration to specify RayTracer mode
enum class RayTracerMode {
    VARIABLE_DIRECTION, // Default mode using Field's vector field
    CONSTANT_DIRECTION  // Fixed direction mode
};

class RayTracer {
public:
    // Constructor for variable direction tracing (default)
    RayTracer(const MeshHandler& mesh_handler, const Field& field_handler);
    
    // Constructor for constant direction tracing
    RayTracer(const MeshHandler& mesh_handler, const Vector3D& fixed_direction);
    
    // Trace a ray starting from a cell and point, up to a maximum number of iterations
    std::vector<CellTrace> traceRay(int start_cell_id, const Vector3D& start_point, int max_iter = 100) const;
    
    // Getter for the associated Field (only relevant in VARIABLE_DIRECTION mode)
    // Returns a reference; ensure it's only called in VARIABLE_DIRECTION mode
    const Field& getField() const { 
        assert(mode_ == RayTracerMode::VARIABLE_DIRECTION && "getField() called in CONSTANT_DIRECTION mode");
        assert(field_ptr_ != nullptr && "field_ptr_ is nullptr in VARIABLE_DIRECTION mode");
        return *field_ptr_; 
    }
    
    // Getter for the fixed direction (only relevant in CONSTANT_DIRECTION mode)
    // Returns a reference; ensure it's only called in CONSTANT_DIRECTION mode
    const Vector3D& getFixedDirection() const { 
        assert(mode_ == RayTracerMode::CONSTANT_DIRECTION && "getFixedDirection() called in VARIABLE_DIRECTION mode");
        return fixed_direction_; 
    }
    
    // Getter for the current mode
    RayTracerMode getMode() const { return mode_; }
    
private:
    const MeshHandler& mesh_;
    RayTracerMode mode_;
    
    // For VARIABLE_DIRECTION mode
    const Field* field_ptr_ = nullptr;
    
    // For CONSTANT_DIRECTION mode
    Vector3D fixed_direction_;
};

#endif // RAY_TRACER_HPP