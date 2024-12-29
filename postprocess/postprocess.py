import argparse
import logging
import os
import vtk

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def load_flux(file_name: str) -> list:
    logging.info(f"Loading flux values from {file_name}")
    with open(file_name, 'r') as f:
        flux = [float(value) for value in f.read().split()]
    logging.info(f"Loaded {len(flux)} flux values.")
    return flux

def load_mesh(mesh_file: str) -> vtk.vtkPolyData:
    logging.info(f"Loading mesh from {mesh_file}")
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(mesh_file)
    reader.Update()
    mesh = reader.GetOutput()
    if not mesh:
        logging.error("Failed to load mesh.")
        raise FileNotFoundError(f"Mesh file {mesh_file} not found or invalid.")
    logging.info("Mesh loaded successfully.")
    return mesh

def add_flux_to_mesh(mesh: vtk.vtkPolyData, flux: list) -> vtk.vtkPolyData:
    logging.info("Adding flux values to the mesh.")
    flux_array = vtk.vtkFloatArray()
    flux_array.SetName("Flux")
    for value in flux:
        flux_array.InsertNextValue(value)
    mesh.GetPointData().AddArray(flux_array)
    mesh.GetPointData().SetActiveScalars("Flux")
    logging.info("Flux values added to the mesh.")
    return mesh

def write_mesh(mesh: vtk.vtkPolyData, output_file: str):
    logging.info(f"Writing updated mesh to {output_file}")
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(mesh)
    writer.Write()
    logging.info("Mesh written successfully.")

def main():
    setup_logging()
    
    parser = argparse.ArgumentParser(description="Postprocess flux data and integrate with VTK mesh.")
    parser.add_argument("--file_name", type=str, required=True, help="Path to the flux.txt file.")
    parser.add_argument("--mesh", type=str, required=True, help="Path to the mesh.vtk file.")
    parser.add_argument("--output_file", type=str, required=False, help="Path to the output file.")
    args = parser.parse_args()
    
    flux = load_flux(args.file_name)
    mesh = load_mesh(args.mesh)
    
    # if len(flux) != mesh.GetNumberOfCells():
    #     logging.error("Number of flux values does not match number of mesh cells.")
    #     logging.error(f"Flux values: {len(flux)}, Mesh cells: {mesh.GetNumberOfCells()}")
    #     raise ValueError("Mismatch between flux values and mesh cells.")
    if args.output_file:
        output_file = args.output_file
    else:
        output_file = os.path.splitext(args.mesh)[0] + "_with_flux.vtk"
    
    updated_mesh = add_flux_to_mesh(mesh, flux)
    write_mesh(updated_mesh, output_file)

if __name__ == "__main__":
    main()