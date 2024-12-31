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
    OutputHandler.cpp
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

# Add compiler warning flags based on the compiler and build type
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    if (CMAKE_BUILD_TYPE STREQUAL "Dev")
        target_compile_options(moc_solver_lib PRIVATE
            -Wall
            -Wextra
            -pedantic
            -Wconversion
            -Wshadow
            -Wnon-virtual-dtor
            -Wold-style-cast
            -Woverloaded-virtual
            -Wfloat-equal
            -Wformat=2
            -Winit-self
            -Wmissing-declarations
            -Wmissing-include-dirs
            -Wredundant-decls
            -Wswitch-enum
        )
    elseif (CMAKE_BUILD_TYPE STREQUAL "perf")
        target_compile_options(moc_solver_lib PRIVATE
            -O2
            -g
            -fno-inline
            -DPROFILE
        )
    endif()
elseif (MSVC)
    if (CMAKE_BUILD_TYPE STREQUAL "Dev")
        target_compile_options(moc_solver_lib PRIVATE
            /W4
        )
    elseif (CMAKE_BUILD_TYPE STREQUAL "perf")
        target_compile_options(moc_solver_lib PRIVATE
            /Zi
            /Od
            /DPROFILE
        )
    endif()
endif()