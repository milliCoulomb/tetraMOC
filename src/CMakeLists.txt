# src/CMakeLists.txt

# Liste des fichiers sources
set(SOURCE_FILES
    Logger.cpp
    Quadrature.cpp
    AngularQuadrature.cpp
    main.cpp
)

# Ajouter l'exécutable principal
add_executable(MOC_Solver ${SOURCE_FILES})

# Lier les bibliothèques MEDCoupling et OpenMP
target_link_libraries(MOC_Solver ${OpenMP_LIBRARIES})
# target_link_libraries(MOC_Solver medcoupling ${OpenMP_LIBRARIES})
