# tests/CMakeLists.txt

# Activer le test
enable_testing()

# Ajouter Google Test
add_subdirectory(${CMAKE_SOURCE_DIR}/external/googletest ${CMAKE_BINARY_DIR}/googletest-build)

# Inclure les répertoires d'en-têtes de Google Test
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Liste des fichiers de test
set(TEST_SOURCES
    test_Vector3D.cpp
    test_Quadrature.cpp
    test_AngularQuadrature.cpp
    test_Tetrahedron.cpp
    test_MeshHandler.cpp
    test_Field.cpp
    test_RayTracer.cpp
    test_GeometryUtils.cpp
    test_RayTracerManager.cpp
    test_FluxSolver.cpp
    test_InputHandler.cpp
)

# Créer l'exécutable des tests
add_executable(runTests ${TEST_SOURCES})

# Lier les bibliothèques nécessaires (Google Test et ta bibliothèque)
target_link_libraries(runTests gtest gtest_main moc_solver_lib)

# Ajouter les tests à CTest
add_test(NAME RunAllTests COMMAND runTests)
