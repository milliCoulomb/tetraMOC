# preprocess_mesh.py
import numpy as np
import medcoupling as mc
import argparse
import logging
import os

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def preprocess_mesh(med_file: str, field_file: str, output_dir: str):
    logger = logging.getLogger(__name__)
    logger.info(f"Loading MED file: {med_file}")
    
    # Load the mesh
    med_file_obj = mc.MEDFileUMesh.New(med_file)
    mesh = med_file_obj.getMeshAtLevel(0)
    
    # Load the field
    field_obj = mc.MEDFileField1D.New(field_file)
    field = field_obj.getFieldAtLevel(0, 0)
    
    # Extract node coordinates
    coords = mesh.getCoords().toNumPyArray().reshape(-1, 3)
    num_nodes = coords.shape[0]
    logger.info(f"Number of nodes: {num_nodes}")
    
    # Extract cell connectivity (tetrahedrons)
    cell_conn = mesh.getNodalConnectivityArray().toNumPyArray()
    if cell_conn.shape[1] != 4:
        logger.error("Mesh is not tetrahedral. Each cell must have 4 nodes.")
        raise ValueError("Non-tetrahedral mesh.")
    num_cells = cell_conn.shape[0]
    logger.info(f"Number of cells: {num_cells}")
    
    # Extract the vector field per cell
    v_field = field.getArray().toNumPyArray().reshape(-1, 3)
    if v_field.shape[0] != num_cells:
        logger.error("Number of vectors in the field does not match the number of cells.")
        raise ValueError("Mismatch between field and mesh.")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Export nodes
    logger.info("Exporting nodes to nodes.txt")
    with open(os.path.join(output_dir, "nodes.txt"), "w") as f:
        f.write(f"{num_nodes}\n")
        for coord in coords:
            f.write(f"{coord[0]} {coord[1]} {coord[2]}\n")
    
    # Export cells
    logger.info("Exporting cells to cells.txt")
    with open(os.path.join(output_dir, "cells.txt"), "w") as f:
        f.write(f"{num_cells}\n")
        for cell in cell_conn:
            f.write(" ".join(map(str, cell)) + "\n")
    
    # Export field
    logger.info("Exporting field to field.txt")
    with open(os.path.join(output_dir, "field.txt"), "w") as f:
        f.write(f"{v_field.shape[0]}\n")
        for vec in v_field:
            f.write(f"{vec[0]} {vec[1]} {vec[2]}\n")
    
    logger.info("Preprocessing completed successfully.")

if __name__ == "__main__":
    setup_logging()
    parser = argparse.ArgumentParser(description="Preprocess MED mesh and field, exporting to text files.")
    parser.add_argument("--med_file", type=str, required=True, help="Path to the MED mesh file.")
    parser.add_argument("--field_file", type=str, required=True, help="Path to the MED field file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for exported files.")
    
    args = parser.parse_args()
    preprocess_mesh(args.med_file, args.field_file, args.output_dir)