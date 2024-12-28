# TetraMOC

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![C++](https://img.shields.io/badge/language-C++-blue.svg)
![CMake](https://img.shields.io/badge/build-CMake-green.svg)

**TetraMOC** is a small C++ project for solving the Boltzmann transport equation using tetrahedral meshes. Designed with modularity and efficiency in mind, TetraMOC provides comprehensive tools for mesh handling, vector operations, angular quadrature, ray tracing, flux solving, and logging utilities. Because the MEDCoupling library was not available without the SALOME platform when this code was written, a preprocessing of the mesh .MED files is needed. The python script will convert the mesh information to .txt files which are read by the C++ code.

TetraMOC reads a YAML file as an input deck, which contains paths to cross-sections, mesh topology preprocessed with Python and solver parameters (number of directions in the angular quadrature, number of rays per boundary face, convergence threshold, etc).

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Modules](#modules)
  - [MeshHandler.hpp](#meshhandlerhpp)
  - [Vector3D.hpp](#vector3dhpp)
  - [Field.hpp](#fieldhpp)
  - [AngularQuadrature.hpp](#angularquadraturehpp)
  - [TrackingData.hpp](#trackingdatahpp)
  - [RayTracer.hpp](#raytracerhpp)
  - [RayTracerManager.hpp](#raytracermanagerhpp)
  - [FluxSolver.hpp](#fluxsolverhpp)
  - [BoltzmannSolver.hpp](#boltzmannsolverhpp)
  - [InputHandler.hpp](#inputhandlerhpp)
  - [Logger.hpp](#loggerhpp)
  - [Quadrature.hpp](#quadraturehpp)
  - [GeometryUtils.hpp](#geometryutilshpp)
  - [Tetrahedron.hpp](#tetrahedronhpp)
- [Contributing](#contributing)
- [License](#license)

## Features

- **Mesh Handling**: Efficiently load and manage tetrahedral meshes.
- **Vector Operations**: Comprehensive 3D vector mathematics.
- **Angular Quadrature**: Generate and manage angular quadrature sets for accurate integrations.
- **Ray Tracing**: Perform ray tracing through meshes for precise flux calculations.
- **Flux Solving**: Compute and collapse flux data seamlessly.
- **Boltzmann Solver**: Robust solver for the Boltzmann transport equation.
- **Logging**: Integrated logging utilities for debugging and informational purposes.
- **Input Handling**: Flexible loading and management of simulation input data.

## Installation

### Prerequisites

- **C++ Compiler**: Supporting C++17 or higher.
- **CMake**: Version 3.10 or higher.
- **Boost Libraries**: Required for certain modules.

### Steps

1. **Clone the Repository**

    ```bash
    git clone https://github.com/yourusername/TetraMOC.git
    cd TetraMOC
    ```

2. **Build the Project**

    ```bash
    mkdir build
    cd build
    cmake ..
    make
    ```

By default, the unit tests are compiled and logging is enabled.
<!-- 3. **Run Tests** (Optional)

    ```bash
    make test
    ``` -->
2. **Install Python venv with MEDCoupling**

    ```bash
    cd TetraMOC
    python -m venv venv
    pip install -r requirements.txt
    ```


## Usage

Refer to the [examples](examples/) directory for sample applications demonstrating the usage of TetraMOC modules. Below is a basic example to get you started:

```cpp
#include <TetraMOC/MeshHandler.hpp>
#include <TetraMOC/BoltzmannSolver.hpp>
#include <TetraMOC/Logger.hpp>

int main() {
    TetraMOC::Logger::init();
    TetraMOC::MeshHandler mesh;
    mesh.loadNodes("nodes.dat");
    mesh.loadCells("cells.dat");
    mesh.loadFaceConnectivity("faces.dat");

    TetraMOC::BoltzmannSolver solver(mesh);
    solver.solve();

    return 0;
}
```

## Modules

### MeshHandler.hpp

Manages mesh data, including nodes, cells, and face connectivity.

- **Key Structures**:
  - `TetraCell`: Represents a tetrahedral cell.
  - `MeshFace`: Represents a mesh face.

- **Key Methods**:
  - `loadNodes()`, `loadCells()`, `loadFaceConnectivity()`: Load mesh data from files.
  - `getCellCenter()`: Compute the center of a cell.

### Vector3D.hpp

Provides a 3D vector class for mathematical operations.

- **Features**:
  - Vector addition, subtraction, scalar multiplication/division.
  - Dot and cross products.
  - Norm calculation and normalization.
  - Equality checks with tolerance.

### Field.hpp

Handles field data, supporting both vector and scalar fields.

- **Features**:
  - Load fields from files.
  - Access and modify field data.
  - Manage direction vectors.

### AngularQuadrature.hpp

Generates and manages angular quadrature sets.

- **Key Structures**:
  - `Direction`: Represents a direction with angles and weight.

- **Key Methods**:
  - `generateQuadrature()`: Generate quadrature points and weights.

### TrackingData.hpp

Stores tracking information for rays traced through the mesh.

- **Key Structures**:
  - `CellTrace`: Information about a single cell traversal.
  - `TrackingData`: Aggregates traces for a single ray.

### RayTracer.hpp

Performs ray tracing through the mesh using either variable or constant directions.

- **Key Features**:
  - Supports both variable and constant direction modes.
  - Traces rays and records cell traversals.

### RayTracerManager.hpp

Manages multiple `RayTracer` instances for comprehensive ray tracing.

- **Key Features**:
  - Initializes ray tracers based on angular quadrature.
  - Aggregates tracking data.
  - Supports parallel processing with OpenMP.

### FluxSolver.hpp

Calculates flux within each cell based on tracking data.

- **Key Features**:
  - Compute flux contributions from traced rays.
  - Collapse directional flux to scalar flux.

### BoltzmannSolver.hpp

Solves the Boltzmann transport equation using the computed flux.

- **Key Features**:
  - Iterative solver with convergence criteria.
  - Computes effective multiplication factor ($k_{eff}$).

### InputHandler.hpp

Loads and manages input data, including cross-sections and nuclear data.

- **Key Features**:
  - Load data from files.
  - Access energy group data.

### Logger.hpp

Provides logging utilities for the application.

- **Features**:
  - Log messages with different severity levels: INFO, WARNING, ERROR.

### Quadrature.hpp

Implements quadrature methods for numerical integration.

- **Key Methods**:
  - Compute Legendre polynomials and their derivatives.
  - Generate Gauss-Legendre and Gauss-Chebyshev quadrature points and weights.

### GeometryUtils.hpp

Provides geometric utilities for mesh operations.

- **Key Functions**:
  - `computeFaceNormal()`: Compute the normal vector of a triangular face.
  - `samplePointOnTriangle()`: Sample a point uniformly within a triangle.

### Tetrahedron.hpp

Represents a tetrahedral element and provides methods for ray exiting.

- **Key Features**:
  - Store vertex positions.
  - Find exit points of rays intersecting the tetrahedron.

## Contributing

We welcome contributions! To contribute, please follow these steps:

1. **Fork the Repository**
2. **Create a Feature Branch**

    ```bash
    git checkout -b feature/YourFeature
    ```

3. **Commit Your Changes**

    ```bash
    git commit -m "Add some feature"
    ```

4. **Push to the Branch**

    ```bash
    git push origin feature/YourFeature
    ```

5. **Open a Pull Request**

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for more details.

## License

This project is licensed under the [MIT License](LICENSE). See the [LICENSE](LICENSE) file for details.