# tests/CMakeLists.txt

# Activate testing
enable_testing()

# Add Google Test as a subdirectory
add_subdirectory(${CMAKE_SOURCE_DIR}/external/googletest ${CMAKE_BINARY_DIR}/googletest-build)

# include header files
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# test sources
set(TEST_SOURCES
    test_Vector3D.cpp
    test_Settings.cpp
    test_Logger.cpp
    test_InputDeck.cpp
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
    test_OutputHandler.cpp
    test_BoltzmannSolver.cpp
)

# create the test executable
add_executable(runTests ${TEST_SOURCES})

# link libraries
target_link_libraries(runTests gtest gtest_main moc_solver_lib)

# add tests
add_test(NAME RunAllTests COMMAND runTests)
