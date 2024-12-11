# tests/CMakeLists.txt

# Activer le test
enable_testing()

# Ajouter Google Test
add_subdirectory(external/googletest)

# Inclure les répertoires d'en-têtes de Google Test
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Liste des fichiers de test
set(TEST_SOURCES
    test_Vector3D.cpp
    test_Quadrature.cpp
    test_AngularQuadrature.cpp
    # Ajoute d'autres fichiers de test si nécessaire
)

# Créer l'exécutable des tests
add_executable(runTests ${TEST_SOURCES})

# Lier les bibliothèques nécessaires
target_link_libraries(runTests gtest gtest_main ${MedCoupling_LIBRARIES} OpenMP::OpenMP_CXX)

# Ajouter les tests à CTest
add_test(NAME RunAllTests COMMAND runTests)
