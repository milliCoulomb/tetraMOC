# src/CMakeLists.txt

# Liste des fichiers sources
set(SOURCE_FILES
    Logger.cpp
    MeshHandler.cpp
    Quadrature.cpp
    AngularQuadrature.cpp
    main.cpp
)

# Ajouter l'exécutable principal
add_executable(MOC_Solver ${SOURCE_FILES})

# Lier les bibliothèques MedCoupling et OpenMP
target_link_libraries(MOC_Solver ${MedCoupling_LIBRARIES} OpenMP::OpenMP_CXX)
