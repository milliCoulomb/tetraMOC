// tests/test_Settings.cpp
#include "Settings.hpp"
#include <gtest/gtest.h>

// Test default constructor
TEST(SettingsTest, DefaultConstructor) {
    Settings settings;
    EXPECT_DOUBLE_EQ(settings.getMultiGroupTolerance(), 1e-7);
    EXPECT_EQ(settings.getMultiGroupMaxIterations(), 1000);
    EXPECT_DOUBLE_EQ(settings.getOneGroupTolerance(), 1e-7);
    EXPECT_EQ(settings.getOneGroupMaxIterations(), 500);
    EXPECT_DOUBLE_EQ(settings.getFissionSourceTolerance(), 1e-7);
    EXPECT_DOUBLE_EQ(settings.getKeffTolerance(), 1e-6);
    EXPECT_EQ(settings.getMaxPowerIterations(), 100);
}

// Test parameterized constructor with valid parameters
TEST(SettingsTest, ParameterizedConstructorValid) {
    Settings settings(1e-4, 2000, 1e-5, 1000, 2e-6, 200, 2e-6);
    EXPECT_DOUBLE_EQ(settings.getMultiGroupTolerance(), 1e-4);
    EXPECT_EQ(settings.getMultiGroupMaxIterations(), 2000);
    EXPECT_DOUBLE_EQ(settings.getOneGroupTolerance(), 1e-5);
    EXPECT_EQ(settings.getOneGroupMaxIterations(), 1000);
    EXPECT_DOUBLE_EQ(settings.getFissionSourceTolerance(), 2e-6);
    EXPECT_EQ(settings.getMaxPowerIterations(), 200);
    EXPECT_DOUBLE_EQ(settings.getKeffTolerance(), 2e-6);
}

// Test parameterized constructor with invalid multi_group_tolerance
TEST(SettingsTest, ParameterizedConstructorInvalidMultiGroupTolerance) {
    EXPECT_THROW(Settings(0.0, 2000, 1e-5, 1000, 2e-6, 200, 2e-6), std::invalid_argument);
    EXPECT_THROW(Settings(-1e-5, 2000, 1e-5, 1000, 2e-6, 200, 2e-6), std::invalid_argument);
}

// Test parameterized constructor with invalid multi_group_max_iterations
TEST(SettingsTest, ParameterizedConstructorInvalidMultiGroupMaxIterations) {
    EXPECT_THROW(Settings(1e-4, 0, 1e-5, 1000, 2e-6, 200, 2e-6), std::invalid_argument);
    EXPECT_THROW(Settings(1e-4, -100, 1e-5, 1000, 2e-6, 200, 2e-6), std::invalid_argument);
}

// Test parameterized constructor with invalid one_group_tolerance
TEST(SettingsTest, ParameterizedConstructorInvalidOneGroupTolerance) {
    EXPECT_THROW(Settings(1e-4, 2000, 0.0, 1000, 2e-6, 200, 2e-6), std::invalid_argument);
    EXPECT_THROW(Settings(1e-4, 2000, -1e-5, 1000, 2e-6, 200, 2e-6), std::invalid_argument);
}

// Test parameterized constructor with invalid one_group_max_iterations
TEST(SettingsTest, ParameterizedConstructorInvalidOneGroupMaxIterations) {
    EXPECT_THROW(Settings(1e-4, 2000, 1e-5, 0, 2e-6, 200, 2e-6), std::invalid_argument);
    EXPECT_THROW(Settings(1e-4, 2000, 1e-5, -1000, 2e-6, 200, 2e-6), std::invalid_argument);
}

// Test parameterized constructor with invalid power_iterations
TEST(SettingsTest, ParameterizedConstructorInvalidPowerIterations) {
    EXPECT_THROW(Settings(1e-4, 2000, 1e-5, 1000, 2e-6, 0, 2e-6), std::invalid_argument);
    EXPECT_THROW(Settings(1e-4, 2000, 1e-5, 1000, 2e-6, -200, 2e-6), std::invalid_argument);
}

// Test getters after parameterized constructor
TEST(SettingsTest, GettersAfterParameterizedConstructor) {
    Settings settings(2e-4, 3000, 2e-5, 1500, 3e-6, 200, 3e-6);
    EXPECT_DOUBLE_EQ(settings.getMultiGroupTolerance(), 2e-4);
    EXPECT_EQ(settings.getMultiGroupMaxIterations(), 3000);
    EXPECT_DOUBLE_EQ(settings.getOneGroupTolerance(), 2e-5);
    EXPECT_EQ(settings.getOneGroupMaxIterations(), 1500);
    EXPECT_DOUBLE_EQ(settings.getFissionSourceTolerance(), 3e-6);
    EXPECT_DOUBLE_EQ(settings.getKeffTolerance(), 3e-6);
    EXPECT_EQ(settings.getMaxPowerIterations(), 200);
}

// Test setters with valid values
TEST(SettingsTest, SettersValid) {
    Settings settings;
    
    // Multi-group tolerance
    settings.setMultiGroupTolerance(5e-5);
    EXPECT_DOUBLE_EQ(settings.getMultiGroupTolerance(), 5e-5);
    
    // Multi-group max iterations
    settings.setMultiGroupMaxIterations(1500);
    EXPECT_EQ(settings.getMultiGroupMaxIterations(), 1500);
    
    // One-group tolerance
    settings.setOneGroupTolerance(5e-7);
    EXPECT_DOUBLE_EQ(settings.getOneGroupTolerance(), 5e-7);
    
    // One-group max iterations
    settings.setOneGroupMaxIterations(750);
    EXPECT_EQ(settings.getOneGroupMaxIterations(), 750);
    
    // Fission source tolerance
    settings.setFissionSourceTolerance(5e-7);
    EXPECT_DOUBLE_EQ(settings.getFissionSourceTolerance(), 5e-7);
    
    // k-effective tolerance
    settings.setKeffTolerance(5e-7);
    EXPECT_DOUBLE_EQ(settings.getKeffTolerance(), 5e-7);
    
    // Max power iterations
    settings.setMaxPowerIterations(150);
    EXPECT_EQ(settings.getMaxPowerIterations(), 150);
}

// Test setters with invalid multi_group_tolerance
TEST(SettingsTest, SetMultiGroupToleranceInvalid) {
    Settings settings;
    EXPECT_THROW(settings.setMultiGroupTolerance(0.0), std::invalid_argument);
    EXPECT_THROW(settings.setMultiGroupTolerance(-1e-5), std::invalid_argument);
}

// Test setters with invalid multi_group_max_iterations
TEST(SettingsTest, SetMultiGroupMaxIterationsInvalid) {
    Settings settings;
    EXPECT_THROW(settings.setMultiGroupMaxIterations(0), std::invalid_argument);
    EXPECT_THROW(settings.setMultiGroupMaxIterations(-100), std::invalid_argument);
}

// Test setters with invalid one_group_tolerance
TEST(SettingsTest, SetOneGroupToleranceInvalid) {
    Settings settings;
    EXPECT_THROW(settings.setOneGroupTolerance(0.0), std::invalid_argument);
    EXPECT_THROW(settings.setOneGroupTolerance(-1e-5), std::invalid_argument);
}

// Test setters with invalid one_group_max_iterations
TEST(SettingsTest, SetOneGroupMaxIterationsInvalid) {
    Settings settings;
    EXPECT_THROW(settings.setOneGroupMaxIterations(0), std::invalid_argument);
    EXPECT_THROW(settings.setOneGroupMaxIterations(-500), std::invalid_argument);
}

// Test setters with invalid fission_source_tolerance
TEST(SettingsTest, SetFissionSourceToleranceInvalid) {
    Settings settings;
    EXPECT_THROW(settings.setFissionSourceTolerance(0.0), std::invalid_argument);
    EXPECT_THROW(settings.setFissionSourceTolerance(-2e-6), std::invalid_argument);
}

// Test setters with invalid keff_tolerance
TEST(SettingsTest, SetKeffToleranceInvalid) {
    Settings settings;
    EXPECT_THROW(settings.setKeffTolerance(0.0), std::invalid_argument);
    EXPECT_THROW(settings.setKeffTolerance(-2e-6), std::invalid_argument);
}

// Test setters with invalid max_power_iterations
TEST(SettingsTest, SetInvalidMaxPowerIterations) {
    Settings settings;
    EXPECT_THROW(settings.setMaxPowerIterations(0), std::invalid_argument);
    EXPECT_THROW(settings.setMaxPowerIterations(-100), std::invalid_argument);
}