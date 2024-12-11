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
)

# Créer l'exécutable des tests
add_executable(runTests ${TEST_SOURCES})

# Lier les bibliothèques nécessaires
target_link_libraries(runTests gtest gtest_main ${OpenMP_LIBRARIES})

# Ajouter les tests à CTest
add_test(NAME RunAllTests COMMAND runTests)
