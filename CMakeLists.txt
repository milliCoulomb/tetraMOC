# CMakeLists.txt (root directory)

cmake_minimum_required(VERSION 3.10)
project(MOC_Solver)

# define C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compiler warning flags based on the compiler
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(ENABLE_LOGGING "Enable logging" ON)

if(NOT ENABLE_LOGGING)
    add_definitions(-DDISABLE_LOGGING)
endif()

# find OpenMP
find_package(OpenMP REQUIRED)
if(NOT OpenMP_FOUND)
    message(FATAL_ERROR "OpenMP not found")
endif()

# Find yaml-cpp
find_package(yaml-cpp REQUIRED)
# if(NOT yaml-cpp_FOUND)
#     message(FATAL_ERROR "yaml-cpp not found")
# endif()

# Add include directories
include_directories(include)

# subdirectories
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(external/yaml-cpp)