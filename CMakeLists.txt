# CMakeLists.txt (root directory)

cmake_minimum_required(VERSION 3.10)
project(MOC_Solver)

# define C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define build types
set(CMAKE_CONFIGURATION_TYPES "Release;Dev" CACHE STRING "" FORCE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Dev)
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

# Add yaml-cpp subdirectory before finding the package
add_subdirectory(external/yaml-cpp)

# # Find yaml-cpp
# find_package(yaml-cpp REQUIRED)

# Add include directories
include_directories(include)

# Set compiler flags based on build type
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL "Dev")
    set(CMAKE_CXX_FLAGS_DEV "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wconversion -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wfloat-equal -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wredundant-decls -Wswitch-enum")
elseif(CMAKE_BUILD_TYPE STREQUAL "perf")
    set(CMAKE_CXX_FLAGS_PERF "${CMAKE_CXX_FLAGS} -O2 -g -fno-inline -DPROFILE")
endif()

# subdirectories
add_subdirectory(src)
add_subdirectory(tests)