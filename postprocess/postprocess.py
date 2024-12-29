import argparse
import logging
import os
import medcoupling as mc
import numpy as np

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def load_flux(file_name: str) -> list:
    logging.info(f"Loading flux values from {file_name}")
    with open(file_name, 'r') as f:
        flux = [float(value) for value in f.read().split()]
    logging.info(f"Loaded {len(flux)} flux values.")
    return flux

def load_mesh(mesh_file: str) -> mc.MEDCouplingUMesh:
    logging.info(f"Loading MED mesh from {mesh_file}")
    try:
        mesh = mc.ReadMeshFromFile(mesh_file)
        logging.info("Mesh loaded successfully.")
        return mesh
    except Exception as e:
        logging.error(f"Failed to load mesh: {e}")
        raise FileNotFoundError(f"Mesh file {mesh_file} not found or invalid.")

def add_flux_to_mesh(mesh: mc.MEDCouplingUMesh, flux: list) -> mc.MEDCouplingUMesh:
    logging.info("Adding flux values to the mesh as a scalar field.")
    num_cells = mesh.getNumberOfCells()
    flux_to_numpy = np.array(flux)
    if flux_to_numpy.shape[0] != num_cells:
        raise ValueError("Concentration array does not match the number of cells in the mesh.")
    field_on_cells = mc.MEDCouplingFieldDouble(mc.ON_CELLS)
    field_on_cells.setMesh(mesh)
    data_array = mc.DataArrayDouble()
    data_array.alloc(num_cells, 1) # one group flux
    for cell_id in range(num_cells):
        data_array[int(cell_id)] = flux[cell_id]
    
    field_on_cells.setArray(data_array)
    field_on_cells.setTime(0.0, 0, -1)
    field_on_cells.setName('FLUX')
    return field_on_cells

def write_mesh(mesh: mc.MEDCouplingUMesh, output_file: str):
    logging.info(f"Writing updated mesh to {output_file}")
    try:
        mc.Write(mesh, output_file, 0)  # 0 for ASCII, use 1 for binary if needed
        logging.info("Mesh written successfully.")
    except Exception as e:
        logging.error(f"Failed to write mesh: {e}")
        raise IOError(f"Could not write mesh to {output_file}.")

def main():
    setup_logging()
    
    parser = argparse.ArgumentParser(description="Postprocess flux data and integrate with MED mesh.")
    parser.add_argument("--file_name", type=str, required=True, help="Path to the flux.txt file.")
    parser.add_argument("--mesh", type=str, required=True, help="Path to the mesh.med file.")
    parser.add_argument("--output_file", type=str, required=False, help="Path to the output MED file.")
    
    args = parser.parse_args()
    
    flux = load_flux(args.file_name)
    mesh = load_mesh(args.mesh)
    
    if args.output_file:
        output_file = args.output_file
    else:
        output_file = os.path.splitext(args.mesh)[0] + "_with_flux.med"
    
    updated_mesh = add_flux_to_mesh(mesh, flux)
    updated_mesh.writeVTK("mesh_with_flux.vtk")

if __name__ == "__main__":
    main()