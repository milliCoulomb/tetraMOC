# TetraMOC

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![C++](https://img.shields.io/badge/language-C++-blue.svg)
![CMake](https://img.shields.io/badge/build-CMake-green.svg)

**TetraMOC** is a small C++ project for solving the Boltzmann transport equation using tetrahedral meshes, with **constant macroscopic cross-section** first and **void boundary conditions**. Designed with modularity and efficiency in mind, TetraMOC provides comprehensive tools for mesh handling, vector operations, angular quadrature, ray tracing, flux solving, and logging utilities. Because the MEDCoupling library was not available without the SALOME platform when this code was written, a preprocessing of the mesh .MED files is needed. The python script will convert the mesh information to .txt files which are read by the C++ code.

TetraMOC reads a YAML file as an input deck, which contains paths to cross-sections, mesh topology preprocessed with Python and solver parameters (number of directions in the angular quadrature, number of rays per boundary face, convergence threshold, etc). Loops over characteristics and directions are parallelized with OpenMP, as well as the ray tracing part.

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
3. **Install Python venv with MEDCoupling**

    ```bash
    cd TetraMOC
    python -m venv venv
    pip install -r requirements.txt
    ```


## Usage

You can use the Python script in the preprocess folder to convert a MED mesh into three text files, one containing nodes IDs with their coordinates in space (nodes.txt), one containing the cell-node connectivity (which node belong to which cell), cells.txt and finally the face to cell connectivity (which face connects which node), faces.txt. To generate these files, use:
```bash
python preprocess_mesh.py --med_file=../examples/cube/cube.med --field_file=../examples/cube/cube.med --output_dir=../examples/cube
```
with *med_file* the path to the .MED file and *output_dir* the directory where the Python code will write the .txt files. Then, create the YAML input deck with the wanted solver parameters and the correct path to .txt files, for example:

```yaml
mesh:
  nodes: "../examples/cube/nodes.txt"
  cells: "../examples/cube/cells.txt"
  faces: "../examples/cube/faces.txt"

cross_sections:
  data_files:
    - "../examples/cube/xs.txt"

angular_quadrature:
  ntheta: 2
  nphi: 4

solver_parameters:
  multi_group_max_iterations: 1000
  multi_group_tolerance: 1e-7
  one_group_max_iterations: 500
  one_group_tolerance: 1e-7
  fission_source_tolerance: 1e-7
  keff_tolerance: 1e-6
  rays_per_face: 8

output:
  flux_output_file: "output/flux.dat"
  k_eff_output_file: "output/k_eff.dat"

logging:
  level: "INFO"
  log_file: "logs/solver.log"
```
and run tetraMOC after compiling it with:
```bash
./src/tetraMOC ../examples/cube/cube.yaml 
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