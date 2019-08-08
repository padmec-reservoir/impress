"""Universidade Federal de Pernambuco."""
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
# from PyTrilinos import Epetra, AztecOO
# import time

"""PRESTO - Python REservoir Simulation TOolbox
   ***************************************************************************
   Author: Ricardo Lira (April/2019)
   Project: Multipoint FLux Approximation with Diamond Stencil (MPFA-D) Solver
   Description: This class is designed for assembling the MPFA-D problem

   Depencencies: Build on top of the Intuitive Multilevel Preprocessor for
   Smart Simulation - IMPRESS"""


class MpfaD:
    def __init__(self, mesh, interpolation_method=None, x=None):
        self.mesh = mesh
        # self.comm = Epetra.PyComm()
        # std_map = Epetra.Map(len(mesh.all_volumes), 0, self.comm)
        # self.T = Epetra.CrsMatrix(Epetra.Copy, std_map, 0)
        self.volumes = mesh.all_volumes
        self.T = lil_matrix((len(self.volumes), len(self.volumes)),
                            dtype=np.float)
        # self.Q = Epetra.Vector(std_map)
        self.Q = lil_matrix((len(self.volumes), 1), dtype=np.float)
        if x is None:
            # self.x = Epetra.Vector(std_map)
            self.x = lil_matrix((len(self.volumes), 1), dtype=np.float)
        else:
            self.x = x
        self.mesh.run_preprocessor()


    def multiply(self, normal_vector, tensor, CD):
        vmv = np.sum(np.dot(normal_vector,
                            tensor) * CD, axis=1)
        area = np.sum(normal_vector * normal_vector, axis=1)
        return vmv / area

    def n_flux_term(self, K_L_n, h_L, face_area, K_R_n=0, h_R=0,
                    boundary=False):
        if boundary:
            K_eq = (1 / h_L)*(face_area * K_L_n)
            return K_eq
        K_eq = (K_R_n * K_L_n)/(K_R_n * h_L + K_L_n * h_R) * face_area
        return K_eq

    def d_flux_term(self, tan, vec, S, h1, Kn1, Kt1, h2=0, Kt2=0, Kn2=0,
                    boundary=False):
        if not boundary:
            mesh_anisotropy_term = (np.sum(tan * vec, axis=1)/(S ** 2))
            physical_anisotropy_term = -((1 / S) * (h1 * (Kt1 / Kn1)
                                         + h2 * (Kt2 / Kn2)))
            cross_diffusion_term = (mesh_anisotropy_term +
                                    physical_anisotropy_term)
            return cross_diffusion_term
        if boundary:
            dot_term = np.sum(-tan * vec, axis=1) * Kn1
            cdf_term = h1 * S * Kt1
            b_cross_difusion_term = (dot_term + cdf_term) / (2 * h1 * S)
            return b_cross_difusion_term

    def assemble_wells(self, source_term):
        self.Q += source_term

    def assemble_transmissibility_matrix(self, perm):
        is_boundary = True
        b_left_volumes_tri = self.mesh.b_left_volumes_tri
        K_L_n = self.mpfad.multiply(self.mesh.b_N_IJK_tri, perm,
                                    self.mesh.b_N_IJK_tri)
        K_L_JI = self.mpfad.multiply(self.mesh.b_N_IJK_tri, perm,
                                     self.mesh.b_tan_JI_tri)
        K_L_JK = self.mpfad.multiply(self.mesh.b_N_IJK_tri, perm,
                                     self.mesh.b_tan_JK_tri)
        D_JK = self.mpfad.d_flux_term(self.mesh.b_tan_JK_tri,
                                      self.mesh.b_LJ_tri, self.mesh.b_area_tri,
                                      self.mesh.b_h_L_tri, K_L_n, K_L_JK,
                                      boundary=is_boundary)
        D_JI = self.mpfad.d_flux_term(self.mesh.b_tan_JI_tri,
                                      self.mesh.b_LJ_tri, self.mesh.b_area_tri,
                                      self.mesh.b_h_L_tri, K_L_n, K_L_JI,
                                      boundary=is_boundary)
        K_eq_tri = self.mpfad.n_flux_term(K_L_n, self.mesh.b_h_L_tri,
                                          self.mesh.b_area_tri,
                                          boundary=is_boundary)

        # b_left_volumes_quad = self.mesh.b_left_volumes_quad
        # b_N_IJK_quad = self.mesh.b_N_IJK_quad
        # b_tan_JI_quad = self.mesh.b_tan_JI_quad
        # b_tan_JK_quad = self.mesh.b_tan_JK_quad
        # b_tan_J2I_quad = self.mesh.b_tan_J2I_quad
        # b_tan_J2K_quad = self.mesh.b_tan_J2K_quad
        # b_LJ_quad = self.mesh.b_LJ_quad
        # b_LJ2_quad = self.mesh.b_LJ2_quad
        # b_h_L_quad = self.mesh.b_h_L_quad
        # b_area_quad = self.mesh.b_area_quad / 2
        # b_i_quad = self.mesh.b_i_quad
        # b_j_quad = self.mesh.b_j_quad
        # b_k_quad = self.mesh.b_k_quad
        # b_j2_quad = self.mesh.b_j2_quad
        # K_L_n = self.mpfad.multiply(b_N_IJK_quad, perm, b_N_IJK_quad)
        # K_L_JI = self.mpfad.multiply(b_N_IJK_quad, perm, b_tan_JI_quad)
        # K_L_JK = self.mpfad.multiply(b_N_IJK_quad, perm, b_tan_JK_quad)
        # K_L_J2I = self.mpfad.multiply(b_N_IJK_quad, perm, b_tan_J2I_quad)
        # K_L_J2K = self.mpfad.multiply(b_N_IJK_quad, perm, b_tan_J2K_quad)
        # D_JK_quad = self.mpfad.d_flux_term(b_tan_JK_quad, b_LJ_quad,
        #                                    b_area_quad, b_h_L_quad,
        #                                    K_L_n, K_L_JK, boundary=is_boundary)
        # D_JI_quad = self.mpfad.d_flux_term(b_tan_JI_quad, b_LJ_quad,
        #                                    b_area_quad, b_h_L_quad,
        #                                    K_L_n, K_L_JI, boundary=is_boundary)
        # D_J2K_quad = self.mpfad.d_flux_term(b_tan_J2K_quad, b_LJ2_quad,
        #                                     b_area_quad, b_h_L_quad,
        #                                     K_L_n, K_L_J2K,
        #                                     boundary=is_boundary)
        # D_J2I_quad = self.mpfad.d_flux_term(b_tan_J2I_quad, b_LJ2_quad,
        #                                     b_area_quad, b_h_L_quad,
        #                                     K_L_n, K_L_J2I,
        #                                     boundary=is_boundary)
        # K_eq_quad = self.mpfad.n_flux_term(K_L_n, b_h_L_quad, b_area_quad,
        #                                    boundary=is_boundary)

        # b_ids_left = np.append(b_left_volumes_tri, b_left_volumes_quad)
        # b_LHS = np.append(K_eq_tri, K_eq_quad)
        # g_I, g_J, g_K = self.mesh.nodes_pressure(b_i_quad, b_j_quad,
        #                                          b_k_quad)
        # RHS_tri = (D_JK * (g_I - g_J) - K_eq_tri * g_J + D_JI * (g_J - g_K))
        # g_I, g_J, g_K, g_J2 = self.mesh.nodes_pressure(b_i_quad, b_j_quad,
        #                                                b_k_quad, b_j2_quad)
        # RHS_quad = ((D_JK_quad - D_J2K_quad) * g_I +
        #             (D_JI_quad - D_J2I_quad) * g_K + (D_JK_quad + D_JI_quad) *
        #             g_J - (D_J2K_quad + D_J2I_quad) * g_J2 - K_eq_quad *
        #             (g_J + g_J2))
        # b_RHS = np.append(RHS_tri, RHS_quad)
        # self.Q[b_ids_left] += -b_RHS
        # self.T.InsertGlobalValues(b_ids_left, b_ids_left, b_LHS)

        is_boundary = False

        # in_faces_tri, in_faces_quad = self.screen_faces_by_verts(self.in_faces)
        # self.in_left_volumes_tri, self.in_right_volumes_tri =\
        #     self.get_left_and_right_volumes(in_faces_tri)
        # self.in_N_IJK_tri, self.in_tan_JI_tri, self.in_tan_JK_tri =\
        #     self.construct_face_vectors(in_faces_tri)
        # self.in_LR_tri, self.in_h_L_tri, self.in_h_R_tri =\
        #     self.get_additional_vectors_and_height(in_faces_tri)
        # self.in_area_tri = self.get_area(in_faces_tri)
        # self.in_i_tri, self.in_j_tri, self.in_k_tri =\
        #     self.get_position_IJK_verts(in_faces_tri)
        #
        # in_faces_tri, in_faces_quad = self.screen_faces_by_verts(self.in_faces)
        # self.in_left_volumes_quad, self.in_right_volumes_quad =\
        #     self.get_left_and_right_volumes(in_faces_quad)
        # self.b_N_IJK_quad, self.b_tan_JI_quad, self.b_tan_JK_quad,\
        #     self.b_tan_J2I_quad, self.b_tan_J2K_quad =\
        #     self.construct_face_vectors(in_faces_quad)
        # self.in_LR_quad, self.in_h_L_quad, self.in_h_R_quad =\
        #     self.get_additional_vectors_and_height(in_faces_quad)
        # self.in_area_quad = self.get_area(in_faces_quad)
        # self.in_i_quad, self.in_j_quad, self.in_k_quad, self.in_j2_quad =\
        #     self.get_position_IJK_verts(in_faces_quad)


        # self.T.InsertGlobalValues(ids_cols, ids_rows, values)

    def _node_treatment(self, args):
        pass

        # go through vertex
        #   call interpolation method
        #   assemble(vertex)

        # self.T.InsertGlobalValues(ids_cols, ids_rows, values)


    def solve_linear_problem(self):
        pass
        # self.T.FillComplete()
        # linearProblem = Epetra.LinearProblem(self.T, self.x, self.Q)
        # solver = AztecOO.AztecOO(linearProblem)
        # solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_gmres)
        # solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_none)
        # solver.Iterate(2000, 1e-16)
        # return self.x
