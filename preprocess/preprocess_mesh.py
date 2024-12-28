# preprocess_mesh.py

import numpy as np
import medcoupling as mc
from typing import Dict, List, Tuple
import argparse
import logging
import os
from collections import defaultdict

def setup_logging():
    """
    Configures the logging settings.
    """
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def build_faces_cell_connectivity(mesh: mc.MEDCouplingUMesh) -> Tuple[Dict[Tuple[int, ...], List[int]], np.ndarray, np.ndarray]:
    """
    Builds the face-to-cell connectivity for the mesh.

    Args:
        mesh (mc.MEDCouplingUMesh): The mesh object.

    Returns:
        Tuple containing:
            - Dictionary mapping sorted face node IDs to cell IDs.
            - NumPy array of reversed connectivity.
            - NumPy array of reversed connectivity indices.
    """
    logger = logging.getLogger(__name__)
    logger.info("Building face-to-cell connectivity.")
    mesh_faces, desc, descIndex, revDesc, revDescIndex = mesh.buildDescendingConnectivity()

    revDesc_np = revDesc.toNumPyArray()
    revDescIndex_np = revDescIndex.toNumPyArray()
    mesh_faces_nodes = mesh_faces.getNodalConnectivity().toNumPyArray().reshape(-1, 4)
    face_node_ids = mesh_faces_nodes[:, 1:]

    faces_cell_connectivity = defaultdict(list)
    num_faces = mesh_faces.getNumberOfCells()

    for face_id in range(num_faces):
        sorted_face = tuple(sorted(face_node_ids[face_id]))
        # print(f"Face {face_id}: {sorted_face}")
        start_idx = revDescIndex_np[face_id]
        end_idx = revDescIndex_np[face_id + 1]
        cell_ids = revDesc_np[start_idx:end_idx]
        faces_cell_connectivity[sorted_face].extend(cell_ids)
        # logger.debug(f"Face {sorted_face} connected to cells {cell_ids}.")

    logger.info("Completed building face-to-cell connectivity.")
    return faces_cell_connectivity, revDesc_np, revDescIndex_np

def export_faces_connectivity(face_dict: Dict[Tuple[int, ...], List[int]], output_dir: str):
    """
    Exports the face connectivity to a text file named 'faces.txt'.

    Each line in the file contains the node indices of the face followed by the adjacent cell IDs.
    If a face is on the boundary, only one cell ID is listed.

    Args:
        face_dict (Dict[Tuple[int, ...], List[int]]): Mapping from face node tuples to adjacent cell IDs.
        output_dir (str): Directory where the output file will be saved.
    """
    logger = logging.getLogger(__name__)
    logger.info("Exporting face connectivity to faces.txt")

    faces_file_path = os.path.join(output_dir, "faces.txt")
    with open(faces_file_path, "w") as f:
        f.write(f"{len(face_dict)}\n")
        for face, cells in face_dict.items():
            # Convert node indices and cell IDs to strings
            face_str = " ".join(map(str, face))
            cells_str = " ".join(map(str, cells))
            # Write the face nodes followed by adjacent cell IDs
            f.write(f"{face_str} {cells_str}\n")

    logger.info("Face connectivity exported successfully.")

def preprocess_mesh(med_file: str, field_file: str, output_dir: str):
    """
    Preprocesses the MED mesh and field files, extracting necessary information and exporting them
    to text files for use in C++ ray tracing code.

    Args:
        med_file (str): Path to the MED mesh file.
        field_file (str): Path to the MED field file.
        output_dir (str): Directory where the output files will be saved.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Loading MED mesh file: {med_file}")
    # field_name = 'VITESSE_ELEM_dom'
    # try:
    #     # Load the field from the MED file
    #     field_obj = mc.ReadField(field_file, field_name, 1, -1)
    #     mesh = field_obj.getMesh()
    # except Exception as e:
    #     logger.error(f"Failed to load field file: {e}")
    #     raise
    try:
        # Load the mesh from the MED file
        mesh = mc.ReadMeshFromFile(med_file)
    except Exception as e:
        logger.error(f"Failed to load mesh file: {e}")
        raise

    # Extract node coordinates
    coords = mesh.getCoords().toNumPyArray().reshape(-1, 3)
    num_nodes = coords.shape[0]
    logger.info(f"Number of nodes: {num_nodes}")
    mesh_faces, desc, descIndex, revDesc, revDescIndex = mesh.buildDescendingConnectivity()
    cell_conn = mesh_faces.getNodalConnectivity().toNumPyArray().reshape(-1, 4)
    # Extract cell connectivity (assuming tetrahedral cells)
    if cell_conn.shape[1] != 4:
        logger.error("Mesh is not tetrahedral. Each cell must have 4 nodes.")
        raise ValueError("Non-tetrahedral mesh.")
    num_cells = cell_conn.shape[0]
    logger.info(f"Number of cells: {num_cells}")
    logger.info(f"Number of faces: {mesh.getNumberOfCells()}")

    # # Extract the velocity field per cell
    # try:
    #     v_field = field_obj.getArray().toNumPyArray().reshape(-1, 3)
    # except Exception as e:
    #     logger.error(f"Failed to extract velocity field: {e}")
    #     raise

    # if v_field.shape[0] != num_cells:
    #     logger.error("Number of velocity vectors does not match the number of cells.")
    #     raise ValueError("Mismatch between velocity field and mesh.")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Export node coordinates to nodes.txt
    logger.info("Exporting nodes to nodes.txt")
    nodes_file_path = os.path.join(output_dir, "nodes.txt")
    # with open(nodes_file_path, "w") as f:
    #     f.write(f"{num_nodes}\n")
    #     for coord in coords:
    #         f.write(f"{coord[0]} {coord[1]} {coord[2]}\n")
    number_of_nodes = coords.shape[0]
    logger.info(f'Number of nodes: {number_of_nodes}')
    with open(nodes_file_path, "w") as f:
        for i in range(number_of_nodes):
            f.write(f"{i} {coords[i][0]} {coords[i][1]} {coords[i][2]}\n")
    

    # Export cell connectivity to cells.txt
    logger.info("Exporting cells to cells.txt")
    cells_file_path = os.path.join(output_dir, "cells.txt")
    # with open(cells_file_path, "w") as f:
    #     f.write(f"{num_cells}\n")
    #     for cell in cell_conn:
    #         # also write the cell idx, and then the nodes idx
    #         f.write(f"{cell} " + " ".join(map(str, cell)) + "\n")
    #         # f.write(" ".join(map(str, cell)) + "\n")
    num_cells = mesh.getNumberOfCells()
    logger.info(f'Number of cells: {num_cells}')
    with open(cells_file_path, "w") as f:
        for i in range(num_cells):
            nodes = mesh.getNodeIdsOfCell(i)
            # print(f"Cell {i}: {nodes}")
            f.write(f"{i} " + " ".join(map(str, nodes)) + "\n")

    # Export velocity field to field.txt
    # logger.info("Exporting velocity field to field.txt")
    # field_file_path = os.path.join(output_dir, "field.txt")
    # with open(field_file_path, "w") as f:
    #     f.write(f"{v_field.shape[0]}\n")
    #     for vec in v_field:
    #         f.write(f"{vec[0]} {vec[1]} {vec[2]}\n")

    # Build and export face connectivity to faces.txt
    logger.info("Building face connectivity.")
    face_dict = build_faces_cell_connectivity(mesh)[0]
    export_faces_connectivity(face_dict, output_dir)

    logger.info("Preprocessing completed successfully.")

if __name__ == "__main__":
    setup_logging()

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Preprocess MED mesh and field, exporting to text files.")
    parser.add_argument("--med_file", type=str, required=True, help="Path to the MED mesh file.")
    parser.add_argument("--field_file", type=str, required=True, help="Path to the MED field file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for exported files.")

    args = parser.parse_args()

    # Run the preprocessing function with provided arguments
    preprocess_mesh(args.med_file, args.field_file, args.output_dir)