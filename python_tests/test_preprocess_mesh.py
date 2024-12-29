import sys
import os
import unittest
from unittest.mock import patch, MagicMock
import medcoupling as mc
import numpy as np

# Add the preprocess directory to the system path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'preprocess'))

import preprocess_mesh

class TestPreprocessMesh(unittest.TestCase):

    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    def test_export_faces_connectivity(self, mock_file):
        face_dict = { (1,2,3): [0,1], (4,5,6): [2] }
        preprocess_mesh.export_faces_connectivity(face_dict, 'output')
        mock_file.assert_called_with(os.path.join('output', 'faces.txt'), 'w')
        handle = mock_file()
        handle.write.assert_any_call("2\n")
        handle.write.assert_any_call("1 2 3 2 0 1\n")
        handle.write.assert_any_call("4 5 6 1 2\n")

    @patch('medcoupling.ReadMeshFromFile')
    @patch('preprocess_mesh.build_faces_cell_connectivity')
    @patch('preprocess_mesh.export_faces_connectivity')
    @patch('os.makedirs')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    def test_preprocess_mesh(self, mock_file, mock_makedirs, mock_export_faces, mock_build_faces, mock_read_mesh):
        # Create a mock mesh object
        mock_mesh = MagicMock()
        mock_mesh.getCoords.return_value.toNumPyArray.return_value = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
        mock_mesh.getNumberOfCells.return_value = 2  # Example number of cells

        # Create a mock for mesh_faces with getNodalConnectivity
        mock_mesh_faces = MagicMock()
        # Assume each cell is a tetrahedron with 4 nodes
        mock_mesh_faces.getNodalConnectivity.return_value.toNumPyArray.return_value = np.array([
            [0, 1, 2, 3],
            [4, 5, 6, 7]
        ])

        # Mock buildDescendingConnectivity to return mesh_faces and other necessary components
        mock_mesh.buildDescendingConnectivity.return_value = (mock_mesh_faces, None, None, None, None)

        # Mock build_faces_cell_connectivity to return a face dictionary and dummy arrays
        mock_build_faces.return_value = ({(1,2,3): [0,1]}, np.array([]), np.array([]))

        # Mock ReadMeshFromFile to return the mock_mesh
        mock_read_mesh.return_value = mock_mesh

        # Call the preprocess_mesh function with mock arguments
        preprocess_mesh.preprocess_mesh('mesh.med', 'field.med', 'output_dir')

        # Assertions to ensure functions are called correctly
        mock_read_mesh.assert_called_with('mesh.med')
        mock_makedirs.assert_called_with('output_dir', exist_ok=True)
        mock_export_faces.assert_called_with({(1,2,3): [0,1]}, 'output_dir')
        self.assertTrue(mock_file.called)

    @patch('medcoupling.ReadMeshFromFile')
    @patch('preprocess_mesh.build_faces_cell_connectivity')
    @patch('preprocess_mesh.export_faces_connectivity')
    @patch('os.makedirs')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    def test_preprocess_mesh_with_non_tetrahedral_mesh(self, mock_file, mock_makedirs, mock_export_faces, mock_build_faces, mock_read_mesh):
        # Create a mock mesh object with non-tetrahedral cells
        mock_mesh = MagicMock()
        # Assuming cells have 5 nodes instead of 4
        mock_mesh.getCoords.return_value.toNumPyArray.return_value = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1], [1,1,1]])
        mock_mesh.getNumberOfCells.return_value = 1  # One cell

        # Create a mock for mesh_faces with getNodalConnectivity
        mock_mesh_faces = MagicMock()
        # Cell with 5 nodes (invalid for tetrahedral)
        mock_mesh_faces.getNodalConnectivity.return_value.toNumPyArray.return_value = np.array([
            [0, 1, 2, 3, 4]
        ])

        # Mock buildDescendingConnectivity to return mesh_faces and other necessary components
        mock_mesh.buildDescendingConnectivity.return_value = (mock_mesh_faces, None, None, None, None)

        # Mock ReadMeshFromFile to return the mock_mesh
        mock_read_mesh.return_value = mock_mesh

        # Call the preprocess_mesh function and expect a ValueError
        with self.assertRaises(ValueError):
            preprocess_mesh.preprocess_mesh('mesh.med', 'field.med', 'output_dir')

        # Ensure that other functions are not called due to the exception
        mock_makedirs.assert_not_called()
        mock_export_faces.assert_not_called()
        self.assertFalse(mock_file.called)  # Changed from assertTrue to assertFalse

if __name__ == '__main__':
    unittest.main()