# preprocess_mesh.py

import numpy as np
import medcoupling as mc
import argparse
import logging
import os
from collections import defaultdict

def setup_logging():
    """
    Configures the logging settings.
    """
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def build_faces_connectivity(cell_conn: np.ndarray) -> Dict[Tuple[int, ...], List[int]]:
    """
    Constructs a dictionary mapping each unique face (sorted tuple of node IDs) to the list of cell IDs
    that share that face.

    Args:
        cell_conn (np.ndarray): Array of cell connectivity, where each row contains node indices of a cell.

    Returns:
        Dict[Tuple[int, ...], List[int]]: Mapping from face node tuples to adjacent cell IDs.
    """
    logger = logging.getLogger(__name__)
    logger.info("Building face connectivity.")

    face_dict = defaultdict(list)

    for cell_id, cell in enumerate(cell_conn):
        # Define the four faces of a tetrahedron
        faces = [
            tuple(sorted([cell[0], cell[1], cell[2]])),
            tuple(sorted([cell[0], cell[1], cell[3]])),
            tuple(sorted([cell[0], cell[2], cell[3]])),
            tuple(sorted([cell[1], cell[2], cell[3]])),
        ]
        for face in faces:
            face_dict[face].append(cell_id)

    logger.info(f"Total unique faces: {len(face_dict)}")
    return face_dict

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

    try:
        # Load the field from the MED file
        field_obj = mc.ReadField(field_file, 1, -1)
        mesh = field_obj.getMesh()
    except Exception as e:
        logger.error(f"Failed to load field file: {e}")
        raise

    # Extract node coordinates
    coords = mesh.getCoords().toNumPyArray().reshape(-1, 3)
    num_nodes = coords.shape[0]
    logger.info(f"Number of nodes: {num_nodes}")

    # Extract cell connectivity (assuming tetrahedral cells)
    cell_conn = mesh.getNodalConnectivity().toNumPyArray()
    if cell_conn.shape[1] != 4:
        logger.error("Mesh is not tetrahedral. Each cell must have 4 nodes.")
        raise ValueError("Non-tetrahedral mesh.")
    num_cells = cell_conn.shape[0]
    logger.info(f"Number of cells: {num_cells}")

    # Extract the velocity field per cell
    try:
        v_field = field_obj.getArray().toNumPyArray().reshape(-1, 3)
    except Exception as e:
        logger.error(f"Failed to extract velocity field: {e}")
        raise

    if v_field.shape[0] != num_cells:
        logger.error("Number of velocity vectors does not match the number of cells.")
        raise ValueError("Mismatch between velocity field and mesh.")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Export node coordinates to nodes.txt
    logger.info("Exporting nodes to nodes.txt")
    nodes_file_path = os.path.join(output_dir, "nodes.txt")
    with open(nodes_file_path, "w") as f:
        f.write(f"{num_nodes}\n")
        for coord in coords:
            f.write(f"{coord[0]} {coord[1]} {coord[2]}\n")

    # Export cell connectivity to cells.txt
    logger.info("Exporting cells to cells.txt")
    cells_file_path = os.path.join(output_dir, "cells.txt")
    with open(cells_file_path, "w") as f:
        f.write(f"{num_cells}\n")
        for cell in cell_conn:
            f.write(" ".join(map(str, cell)) + "\n")

    # Export velocity field to field.txt
    logger.info("Exporting velocity field to field.txt")
    field_file_path = os.path.join(output_dir, "field.txt")
    with open(field_file_path, "w") as f:
        f.write(f"{v_field.shape[0]}\n")
        for vec in v_field:
            f.write(f"{vec[0]} {vec[1]} {vec[2]}\n")

    # Build and export face connectivity to faces.txt
    logger.info("Building face connectivity.")
    face_dict = build_faces_connectivity(cell_conn)
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