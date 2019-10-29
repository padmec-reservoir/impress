"""This is the mpfad-3D preprocessor Tests."""
import unittest
import numpy as np
from mesh_preprocessor import Mesh


class PreMpfaDTest(unittest.TestCase):

    def setUp(self):
        K_1 = [1.0, 0.0, 0.0,
               0.0, 1.0, 0.0,
               0.0, 0.0, 1.0]
        self.mesh = Mesh('mesh/test_mesh_5_vols.h5m', dim=3)
        self.mesh.set_boundary_conditions('Dirichlet', {101: 0.0})
        self.mesh.set_material_prop('permeability', {1: K_1})

    def tearDown(self):
        del self.mesh

    def test_preprocessor_class_should_be_none(self):
        """Test class initiatilization."""
        mesh = self.mesh
        self.assertIsNotNone(mesh)

    # def test_get_all_volumes(self):
    #     """Test class get all volumes."""
    #     volumes = self.mesh.all_volumes
    #     self.assertEqual(len(volumes), 9)
    #
    # def test_get_all_internal_faces(self):
    #     """Test class get all internal faces."""
    #     in_faces = self.mesh.in_faces
    #     self.assertEqual(len(in_faces), 12)
    #
    # def test_get_all_boundary_faces(self):
    #     """Test class get all boundary faces."""
    #     b_faces = self.mesh.b_faces
    #     self.assertEqual(len(b_faces), 16)
    #
    # def test_get_all_internal_verts(self):
    #     """Test class get all internal verts."""
    #     in_verts = self.mesh.in_verts
    #     self.assertEqual(len(in_verts), 0)
    #
    # def test_get_all_boundary_verts(self):
    #     """Test class get all boundary verts."""
    #     b_verts = self.mesh.b_verts
    #     self.assertEqual(len(b_verts), 13)
    #
    # def test_screening_faces_gets_only_quadrilateral_intern_faces(self):
    #     """Test class to get only quadrilateral faces."""
    #     in_faces = self.mesh.in_faces
    #     _, in_quad_faces = self.mesh.screen_faces_by_verts(in_faces)
    #     self.assertEqual(len(in_quad_faces),  1)
    #
    # def test_screening_faces_gets_only_triangular_intern_faces(self):
    #     """Test class to get only quadrilateral faces."""
    #     in_faces = self.mesh.in_faces
    #     tri_faces, _ = self.mesh.screen_faces_by_verts(in_faces)
    #     self.assertEqual(len(tri_faces),  11)
    #
    # def test_screening_faces_gets_only_quadrilateral_boundary_faces(self):
    #     """Test class to get only quadrilateral faces."""
    #     in_faces = self.mesh.b_faces
    #     _, in_quad_faces = self.mesh.screen_faces_by_verts(in_faces)
    #     self.assertEqual(len(in_quad_faces),  6)
    #
    # def test_screening_faces_gets_only_triangular_boundary_faces(self):
    #     """Test class to get only quadrilateral faces."""
    #     in_faces = self.mesh.b_faces
    #     tri_faces, _ = self.mesh.screen_faces_by_verts(in_faces)
    #     self.assertEqual(len(tri_faces),  10)
    #
    # def test_if_permeability_tensor_is_assigned(self):
    #     """Test if permeability tensor is being assigned for all volumes."""
    #     all_volumes = self.mesh.all_volumes
    #     K_1 = [1.0, 0.0, 0.0,
    #            0.0, 1.0, 0.0,
    #            0.0, 0.0, 1.0]
    #     for volume in all_volumes:
    #         perm = self.mesh.M.permeability[[volume]][0]
    #         self.assertListEqual(list(perm), K_1)
    #
    # def test_get_dirichlet_faces(self):
    #     """Test if Dirichlet BC is implemented."""
    #     d_faces = self.mesh.M.faces.flag[101]
    #     dirichlet_faces_pressure = self.mesh.M.dirichlet_faces
    #     for face in d_faces:
    #         b_pressure = dirichlet_faces_pressure[[face]]
    #         self.assertEqual(b_pressure, 0.)
    #
    # def test_position_left_and_right_volumes_for_triangular_intern_faces(self):
    #     """Test if all volumes are properly oriented with the normal."""
    #     in_faces = self.mesh.in_faces
    #     tri_faces = self.mesh.screen_faces_by_verts(in_faces)[0]
    #     (left_volumes,
    #      right_volumes) = self.mesh.get_left_and_right_volumes(tri_faces)
    #     N_IJK = self.mesh.construct_face_vectors(tri_faces)[0]
    #     vector_left_to_right = (self.mesh.M.volumes.center[left_volumes]
    #                             - self.mesh.M.volumes.center[right_volumes])
    #     is_positive = np.sum(vector_left_to_right * N_IJK, axis=1)
    #     self.assertFalse(list(np.flatnonzero(is_positive < 0)))
    #
    # def test_position_left_and_right_volumes_for_quad_intern_faces(self):
    #     """Test if all volumes are properly oriented with the normal."""
    #     in_faces = self.mesh.in_faces
    #     quad_faces = self.mesh.screen_faces_by_verts(in_faces)[1]
    #     (left_volumes,
    #      right_volumes) = self.mesh.get_left_and_right_volumes(quad_faces)
    #     N_IJK = self.mesh.construct_face_vectors(quad_faces)[0]
    #     vector_left_to_right = (self.mesh.M.volumes.center[left_volumes]
    #                             - self.mesh.M.volumes.center[right_volumes])
    #     is_positive = np.sum(vector_left_to_right * N_IJK, axis=1)
    #     self.assertFalse(list(np.flatnonzero(is_positive < 0)))
    #
    # def test_get_normal_face_area_vector_points_outward(self):
    #     """Test if bondary volumes have their normal outward oriented."""
    #     quad_faces = self.mesh.screen_faces_by_verts(self.mesh.b_faces)[1]
    #     N_IJK = self.mesh.construct_face_vectors(quad_faces, boundary=True)[0]
    #     adj_vol = self.mesh.M.faces.bridge_adjacencies(quad_faces, 2, 3)
    #     i = self.mesh.get_position_IJK_verts(quad_faces)[0]
    #     vector_left_to_right = (self.mesh.M.nodes.coords[i]
    #                             - self.mesh.M.volumes.center[adj_vol])
    #     is_positive = np.sum(vector_left_to_right * N_IJK, axis=1)
    #     self.assertFalse(list(np.flatnonzero(is_positive < 0)))
    #
    # def test_calculate_area_vector(self):
    #     b_faces = self.mesh.b_faces
    #     tri_faces, quad_faces = self.mesh.screen_faces_by_verts(b_faces)
    #     face_area_quad = self.mesh.get_area(quad_faces)
    #     face_area_tri = self.mesh.get_area(tri_faces)
    #     self.assertEqual(np.sum(face_area_tri) + np.sum(face_area_quad), 6.)
    #
    # def test_calculate_additional_geometric_information(self):
    #     in_faces = self.mesh.in_faces
    #     quad_faces = self.mesh.screen_faces_by_verts(in_faces)[1]
    #     LR, h_L, h_R = self.mesh.get_additional_vectors_and_height(quad_faces)
    #     self.assertEqual(LR[0][0], h_L + h_R)
    #
    # @unittest.skip("Test not ready. Still missing a test case assertion")
    # def test_calculate_additional_geometric_information_boundary(self):
    #     b_faces = self.mesh.b_faces
    #     quad_faces = self.mesh.screen_faces_by_verts(b_faces)[1]
    #     LJ, h_L = self.mesh.get_additional_vectors_and_height(quad_faces,
    #                                                           boundary=True)
    #     pass
    #
    # def test_calculate_volumes(self):
    #     all_vol_volumes = self.mesh.get_volume(self.mesh.all_volumes)
    #     # print(all_vol_volumes)
    #     self.assertEqual(sum(all_vol_volumes), 1.)
    #
    # def test_run_preprocessor(self):
    #     try:
    #         self.mesh.run_preprocessor()
    #     except:
    #         self.fail("mesh.preprocessor did not run.")


if __name__ == "__main__":
    unittest.main()
