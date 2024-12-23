# TetraMOC

TetraMOC is a C++ library designed for solving the Boltzmann transport equation using tetrahedral meshes. It includes modules for mesh handling, vector operations, angular quadrature, ray tracing, flux solving, and logging utilities.

## Table of Contents

- Features
- Installation
- Usage
- Modules
  - MeshHandler.hpp
  - Vector3D.hpp
  - Field.hpp
  - AngularQuadrature.hpp
  - TrackingData.hpp
  - RayTracer.hpp
  - RayTracerManager.hpp
  - FluxSolver.hpp
  - BoltzmannSolver.hpp
  - InputHandler.hpp
  - Logger.hpp
  - Quadrature.hpp
  - GeometryUtils.hpp
  - Tetrahedron.hpp
- Contributing
- License

## Features

- **Mesh Handling**: Load and manage tetrahedral meshes.
- **Vector Operations**: 3D vector mathematics.
- **Angular Quadrature**: Generate and manage angular quadrature sets.
- **Ray Tracing**: Trace rays through the mesh for flux calculations.
- **Flux Solving**: Compute and collapse flux data.
- **Boltzmann Solver**: Solve the Boltzmann transport equation.
- **Logging**: Integrated logging utilities for debugging and information.
- **Input Handling**: Load and manage input data for simulations.

## Installation

1. **Clone the Repository**:
    ```shell
    git clone https://github.com/yourusername/TetraMOC.git
    cd TetraMOC
    ```

2. **Build the Project**:
    ```shell
    mkdir build
    cd build
    cmake ..
    make
    ```

## Usage

Refer to the examples directory for sample applications using TetraMOC modules.

## Modules

### 

MeshHandler.hpp



Manages mesh data including nodes, cells, and face connectivity.

**Key Structures**:
- 

TetraCell

: Represents tetrahedral cells.
- 

MeshFace

: Represents mesh faces.

**Key Methods**:
- `loadNodes()`, `loadCells()`, `loadFaceConnectivity()`: Load mesh data from files.
- `getCellCenter()`: Compute the center of a cell.

### 

Vector3D.hpp



Provides a 3D vector class for mathematical operations.

**Features**:
- Vector addition, subtraction, scalar multiplication/division.
- Dot and cross products.
- Norm calculation and normalization.
- Equality checks with tolerance.

### 

Field.hpp



Handles field data, supporting both vector and scalar fields.

**Features**:
- Load fields from files.
- Access and modify field data.
- Manage direction vectors.

### 

AngularQuadrature.hpp



Generates and manages angular quadrature sets.

**Key Structures**:
- 

Direction

: Represents a direction with angles and weight.

**Key Methods**:
- `generateQuadrature()`: Generate quadrature points and weights.

### 

TrackingData.hpp



Stores tracking information for rays traced through the mesh.

**Key Structures**:
- 

CellTrace

: Information about a single cell traversal.
- 

TrackingData

: Aggregates traces for a single ray.

### 

RayTracer.hpp



Performs ray tracing through the mesh using either variable or constant directions.

**Key Features**:
- Supports both variable and constant direction modes.
- Traces rays and records cell traversals.

### 

RayTracerManager.hpp



Manages multiple 

RayTracer

 instances for comprehensive ray tracing.

**Key Features**:
- Initializes ray tracers based on angular quadrature.
- Aggregates tracking data.
- Supports parallel processing with OpenMP.

### 

FluxSolver.hpp



Calculates flux within each cell based on tracking data.

**Key Features**:
- Compute flux contributions from traced rays.
- Collapse directional flux to scalar flux.

### 

BoltzmannSolver.hpp



Solves the Boltzmann transport equation using the computed flux.

**Key Features**:
- Iterative solver with convergence criteria.
- Computes effective multiplication factor (k_eff).

### 

InputHandler.hpp



Loads and manages input data including cross-sections and nuclear data.

**Key Features**:
- Load data from files.
- Access energy group data.

### 

Logger.hpp



Provides logging utilities for the application.

**Features**:
- Log messages with different severity levels: INFO, WARNING, ERROR.

### 

Quadrature.hpp



Implements quadrature methods for numerical integration.

**Key Methods**:
- Compute Legendre polynomials and their derivatives.
- Generate Gauss-Legendre and Gauss-Chebyshev quadrature points and weights.

### 

GeometryUtils.hpp



Provides geometric utilities for mesh operations.

**Key Functions**:
- 

computeFaceNormal()

: Compute the normal vector of a triangle face.
- 

samplePointOnTriangle()

: Sample a point uniformly within a triangle.

### 

Tetrahedron.hpp



Represents a tetrahedral element and provides methods for ray exiting.

**Key Features**:
- Store vertex positions.
- Find exit points of rays intersecting the tetrahedron.

## Contributing

1. **Fork the Repository**
2. **Create a Feature Branch**: `git checkout -b feature/YourFeature`
3. **Commit Your Changes**: `git commit -m 'Add some feature'`
4. **Push to the Branch**: `git push origin feature/YourFeature`
5. **Open a Pull Request**

Refer to CONTRIBUTING.md for more details.

## License

This project is licensed under the MIT License. See the LICENSE file for details.