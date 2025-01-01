import sys
import os
import unittest
from unittest.mock import patch, mock_open

# Add the postprocess directory to the system path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'postprocess'))

import postprocess
import medcoupling as mc
import numpy as np
import vtk

class TestPostprocess(unittest.TestCase):

    @patch('builtins.open', new_callable=mock_open, read_data='1.0 2.0 3.0')
    def test_load_flux(self, mock_file):
        flux = postprocess.load_flux('flux.txt')
        self.assertEqual(flux, [1.0, 2.0, 3.0])
        mock_file.assert_called_with('flux.txt', 'r')

    @patch('medcoupling.ReadMeshFromFile')
    def test_load_mesh_success(self, mock_read_mesh):
        mock_mesh = mc.MEDCouplingUMesh()
        mock_read_mesh.return_value = mock_mesh
        mesh = postprocess.load_mesh('mesh.med')
        self.assertEqual(mesh, mock_mesh)
        mock_read_mesh.assert_called_with('mesh.med')

    @patch('medcoupling.ReadMeshFromFile', side_effect=Exception('Error'))
    def test_load_mesh_failure(self, mock_read_mesh):
        with self.assertRaises(FileNotFoundError):
            postprocess.load_mesh('invalid.med')

    def test_add_flux_to_mesh(self):
        mock_mesh = mc.MEDCouplingUMesh()
        mock_mesh.getNumberOfCells = lambda: 3
        flux = [1.0, 2.0, 3.0]
        field = postprocess.add_flux_to_mesh(mock_mesh, flux)
        self.assertEqual(field.getName(), 'FLUX')
        self.assertEqual(field.getNumberOfComponents(), 1)

if __name__ == '__main__':
    unittest.main()