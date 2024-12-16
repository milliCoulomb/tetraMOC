# src/CMakeLists.txt

# Liste des fichiers sources (exclure main.cpp)
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
)

# Créer une bibliothèque statique
add_library(moc_solver_lib STATIC ${SOURCE_FILES})

# Inclure les répertoires d'en-têtes
target_include_directories(moc_solver_lib PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Lier MEDCoupling et OpenMP à la bibliothèque
target_link_libraries(moc_solver_lib PUBLIC OpenMP::OpenMP_CXX)
