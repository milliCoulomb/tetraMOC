# preprocess_vtk.py

import vtk
import numpy as np
import argparse
import logging
import os
from collections import defaultdict
from typing import Dict, List, Tuple

def setup_logging():
    """
    Configures the logging settings.
    """
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def build_faces_cell_connectivity(unstructured_grid) -> Tuple[Dict[Tuple[int, ...], List[int]], np.ndarray, np.ndarray]:
    """
    Builds the face-to-cell connectivity for the unstructured tetrahedral mesh.

    Args:
        unstructured_grid (vtk.vtkUnstructuredGrid): The VTK unstructured grid object.

    Returns:
        Tuple containing:
            - Dictionary mapping sorted face node IDs to cell IDs.
            - NumPy array of face IDs.
            - NumPy array of neighboring cell IDs.
    """
    logger = logging.getLogger(__name__)
    logger.info("Building face-to-cell connectivity.")

    face_dict = defaultdict(list)
    num_cells = unstructured_grid.GetNumberOfCells()

    for cell_id in range(num_cells):
        cell = unstructured_grid.GetCell(cell_id)
        if cell.GetCellType() != vtk.VTK_TETRA:
            logger.error("Non-tetrahedral cell found.")
            raise ValueError("Mesh contains non-tetrahedral cells.")
        points = cell.GetPointIds()
        faces = [
            (points.GetId(0), points.GetId(1), points.GetId(2)),
            (points.GetId(0), points.GetId(1), points.GetId(3)),
            (points.GetId(0), points.GetId(2), points.GetId(3)),
            (points.GetId(1), points.GetId(2), points.GetId(3)),
        ]
        for face in faces:
            sorted_face = tuple(sorted(face))
            face_dict[sorted_face].append(cell_id)

    logger.info("Completed building face-to-cell connectivity.")
    face_ids = np.arange(len(face_dict))
    cell_neighbors = np.array([neighbors for neighbors in face_dict.values()])
    return face_dict, face_ids, cell_neighbors

def export_nodes(unstructured_grid, output_dir: str):
    """
    Exports node coordinates to nodes.txt.

    Args:
        unstructured_grid (vtk.vtkUnstructuredGrid): The VTK unstructured grid object.
        output_dir (str): Directory where the output file will be saved.
    """
    logger = logging.getLogger(__name__)
    logger.info("Exporting nodes to nodes.txt")

    nodes_file_path = os.path.join(output_dir, "nodes.txt")
    points = unstructured_grid.GetPoints()
    num_nodes = points.GetNumberOfPoints()

    with open(nodes_file_path, "w") as f:
        f.write(f"{num_nodes}\n")
        for node_id in range(num_nodes):
            x, y, z = points.GetPoint(node_id)
            f.write(f"{node_id} {x} {y} {z}\n")

    logger.info("Nodes exported successfully.")

def export_cells(unstructured_grid, output_dir: str):
    """
    Exports cell connectivity to cells.txt.

    Args:
        unstructured_grid (vtk.vtkUnstructuredGrid): The VTK unstructured grid object.
        output_dir (str): Directory where the output file will be saved.
    """
    logger = logging.getLogger(__name__)
    logger.info("Exporting cells to cells.txt")

    cells_file_path = os.path.join(output_dir, "cells.txt")
    num_cells = unstructured_grid.GetNumberOfCells()

    with open(cells_file_path, "w") as f:
        f.write(f"{num_cells}\n")
        for cell_id in range(num_cells):
            cell = unstructured_grid.GetCell(cell_id)
            point_ids = cell.GetPointIds()
            nodes = [point_ids.GetId(i) for i in range(point_ids.GetNumberOfIds())]
            nodes_str = " ".join(map(str, nodes))
            f.write(f"{cell_id} {nodes_str}\n")

    logger.info("Cells exported successfully.")

def export_face_to_cell(face_dict: Dict[Tuple[int, ...], List[int]], output_dir: str):
    """
    Exports face connectivity to face_to_cell.txt.

    Args:
        face_dict (Dict[Tuple[int, ...], List[int]]): Mapping from face node tuples to adjacent cell IDs.
        output_dir (str): Directory where the output file will be saved.
    """
    logger = logging.getLogger(__name__)
    logger.info("Exporting face connectivity to face_to_cell.txt")

    face_file_path = os.path.join(output_dir, "face_to_cell.txt")
    with open(face_file_path, "w") as f:
        for face_id, (face, cells) in enumerate(face_dict.items()):
            num_neighbors = len(cells)
            cells_str = " ".join(map(str, cells))
            f.write(f"{face_id} {num_neighbors} {cells_str}\n")

    logger.info("Face connectivity exported successfully.")

def preprocess_mesh(vtk_file: str, output_dir: str):
    """
    Preprocesses the VTK unstructured tetrahedral mesh, extracting necessary information and exporting them
    to text files.

    Args:
        vtk_file (str): Path to the VTK mesh file.
        output_dir (str): Directory where the output files will be saved.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Loading VTK mesh file: {vtk_file}")

    # Read the VTK file
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_file)
    reader.Update()
    unstructured_grid = reader.GetOutput()

    if unstructured_grid is None:
        logger.error("Failed to read the VTK file.")
        raise IOError("Unable to read the VTK file.")

    num_nodes = unstructured_grid.GetNumberOfPoints()
    num_cells = unstructured_grid.GetNumberOfCells()
    logger.info(f"Number of nodes: {num_nodes}")
    logger.info(f"Number of cells: {num_cells}")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Export nodes and cells
    export_nodes(unstructured_grid, output_dir)
    export_cells(unstructured_grid, output_dir)

    # Build and export face connectivity
    face_dict, face_ids, cell_neighbors = build_faces_cell_connectivity(unstructured_grid)
    export_face_to_cell(face_dict, output_dir)

    logger.info("Preprocessing completed successfully.")

if __name__ == "__main__":
    setup_logging()

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Preprocess VTK unstructured tetrahedral mesh and export to text files.")
    parser.add_argument("--vtk_file", type=str, required=True, help="Path to the VTK mesh file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for exported files.")

    args = parser.parse_args()

    # Run the preprocessing function with provided arguments
    preprocess_mesh(args.vtk_file, args.output_dir)