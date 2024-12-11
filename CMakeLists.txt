# CMakeLists.txt (racine)

cmake_minimum_required(VERSION 3.10)
project(SN_Solver)

# Définir le standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Activer les optimisations
set(CMAKE_BUILD_TYPE Release)

# Trouver MedCoupling
find_package(MedCoupling REQUIRED)
if(NOT MedCoupling_FOUND)
    message(FATAL_ERROR "MedCoupling not found")
endif()

# Trouver OpenMP
find_package(OpenMP REQUIRED)
if(NOT OpenMP_FOUND)
    message(FATAL_ERROR "OpenMP not found")
endif()

# Ajouter les répertoires d'en-têtes
include_directories(${MedCoupling_INCLUDE_DIRS} include)

# Ajouter les sous-répertoires
add_subdirectory(src)
add_subdirectory(tests)