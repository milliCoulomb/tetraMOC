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

def add_flux_to_mesh(mesh: mc.MEDCouplingUMesh, flux: list, number_of_groups: int = 6) -> mc.MEDCouplingUMesh:
    logging.info("Adding flux values to the mesh as a scalar field.")
    num_cells = mesh.getNumberOfCells()
    flux_to_numpy = np.array(flux)
    # reshape the flux array so that it has shape (num_cells, number_of_groups)
    flux_to_numpy_reshape = flux_to_numpy
    logging.info(f"Flux array shape: {flux_to_numpy_reshape.shape}")
    # if there is only one group, we need to add an extra dimension
    if number_of_groups == 1:
        flux_to_numpy_reshape = flux_to_numpy_reshape[np.newaxis, :]
    field_on_cells = mc.MEDCouplingFieldDouble(mc.ON_CELLS)
    field_on_cells.setMesh(mesh)
    data_array = mc.DataArrayDouble()
    logging.info(f"Allocating data array with shape ({num_cells}, {number_of_groups})")
    data_array.alloc(num_cells, number_of_groups)
    for cell_id in range(num_cells):
            data_array[int(cell_id)] = tuple(flux_to_numpy_reshape[:, cell_id])
    field_on_cells.setArray(data_array)
    print(field_on_cells)
    print(data_array)
    field_on_cells.setTime(0.0, 0, -1)
    field_on_cells.setName('FLUX')
    return field_on_cells

def main():
    setup_logging()
    
    parser = argparse.ArgumentParser(description="Postprocess flux data and integrate with MED mesh.")
    parser.add_argument("--file_name", type=str, required=True, help="Path to the flux.txt file.")
    parser.add_argument("--mesh", type=str, required=True, help="Path to the mesh.med file.")
    parser.add_argument("--number_of_groups", type=int, required=False, default=1, help="Number of groups in the flux data.")
    parser.add_argument("--output_file", type=str, required=False, help="Path to the output MED file.")
    
    args = parser.parse_args()
    
    flux = load_flux(args.file_name)
    flux = np.genfromtxt(args.file_name, delimiter=' ')
    mesh = load_mesh(args.mesh)
    
    if args.output_file:
        output_file = args.output_file
    else:
        output_file = os.path.splitext(args.mesh)[0] + "_with_flux"
    
    updated_mesh = add_flux_to_mesh(mesh, flux , args.number_of_groups)
    updated_mesh.writeVTK(output_file)
    updated_mesh.write(output_file + ".med")

if __name__ == "__main__":
    main()