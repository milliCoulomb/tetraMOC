# src/CMakeLists.txt

# Liste des fichiers sources
set(SOURCE_FILES
    Vector3D.cpp
    Quadrature.cpp
    AngularQuadrature.cpp
    main.cpp
)

# Ajouter l'exécutable principal
add_executable(SN_Solver ${SOURCE_FILES})

# Lier les bibliothèques MEDCoupling et OpenMP
target_link_libraries(SN_Solver medcoupling ${OpenMP_LIBRARIES})
