# CMakeLists.txt (racine)

cmake_minimum_required(VERSION 3.10)
project(MOC_Solver)

# Définir le standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Activer les optimisations
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(ENABLE_LOGGING "Enable logging" ON)

if(NOT ENABLE_LOGGING)
    add_definitions(-DDISABLE_LOGGING)
endif()

# Trouver OpenMP
find_package(OpenMP REQUIRED)
if(NOT OpenMP_FOUND)
    message(FATAL_ERROR "OpenMP not found")
endif()

# Find yaml-cpp
find_package(yaml-cpp REQUIRED)
if(NOT yaml-cpp_FOUND)
    message(FATAL_ERROR "yaml-cpp not found")
endif()

# Ajouter les répertoires d'en-têtes
include_directories(include)

# Ajouter MEDCoupling
# add_subdirectory(external/medcoupling)

# Ajouter les sous-répertoires
add_subdirectory(src)
add_subdirectory(tests)