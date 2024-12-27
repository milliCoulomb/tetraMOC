# src/CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(moc_solver)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# List of source files (excluding main.cpp)
set(SOURCE_FILES
    Logger.cpp
    Quadrature.cpp
    AngularQuadrature.cpp
    Tetrahedron.cpp
    MeshHandler.cpp
    Field.cpp
    RayTracer.cpp
    GeometryUtils.cpp
    RayTracerManager.cpp
    FluxSolver.cpp
    InputHandler.cpp
    InputDeck.cpp
    BoltzmannSolver.cpp
)

# Create a static library
add_library(moc_solver_lib STATIC ${SOURCE_FILES})

# Include header directories
target_include_directories(moc_solver_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Link libraries
target_link_libraries(moc_solver_lib PUBLIC 
    OpenMP::OpenMP_CXX 
    yaml-cpp
)

# Add executable for main.cpp
add_executable(tetraMOC main.cpp)

# Link the executable with the library and dependencies
target_link_libraries(tetraMOC PUBLIC 
    moc_solver_lib
    OpenMP::OpenMP_CXX 
    yaml-cpp
)

# Add compiler warning flags based on the compiler
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(moc_solver_lib PRIVATE
        -Wall             # Enable all basic warnings
        -Wextra           # Enable extra warnings
        -pedantic         # Enforce strict ISO compliance
        -Wconversion      # Warn about implicit type conversions
        -Wshadow          # Warn when a local variable shadows another
        -Wnon-virtual-dtor# Warn if a class with virtual functions doesn't have a virtual destructor
        -Wold-style-cast  # Warn about C-style casts
        -Woverloaded-virtual # Warn when a function declaration hides a virtual function in a base class
        -Wfloat-equal     # Warn about comparisons between floating-point numbers for equality
        -Wformat=2        # Warn about printf/scanf format string issues
        -Winit-self       # Warn if a variable is initialized with itself
        -Wmissing-declarations # Warn if a global function is missing a declaration
        -Wmissing-include-dirs # Warn if an #include directive refers to a directory not listed with -I
        -Wredundant-decls # Warn about redundant declarations
        -Wswitch-enum     # Warn if a switch statement doesn't handle all enumeration values
        # -Werror         # Uncomment to treat all warnings as errors
    )
elseif (MSVC)
    target_compile_options(moc_solver_lib PRIVATE
        /W4               # Enable high-level warnings
        # /WX             # Uncomment to treat warnings as errors
    )
endif()