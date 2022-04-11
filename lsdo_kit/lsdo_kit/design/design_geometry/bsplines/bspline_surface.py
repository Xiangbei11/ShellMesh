import numpy as np

import scipy.sparse as sps
from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
from lsdo_kit.cython.surface_projection_py import compute_surface_projection


class BSplineSurface:
    def __init__(self, name, order_u, order_v, knots_u, knots_v, shape, control_points):
        self.name = name
        self.order_u = order_u
        self.order_v = order_v
        self.knots_u = knots_u
        self.knots_v = knots_v
        self.shape = shape
        self.starting_geometry_index = None
        self.control_points = control_points

    def compute_eval_map_points(self, u_vec, v_vec):
        data = np.zeros(len(u_vec) * self.order_u * self.order_v)
        row_indices = np.zeros(len(data), np.int32)
        col_indices = np.zeros(len(data), np.int32)

        num_control_points = self.shape[0] * self.shape[1]

        get_basis_surface_matrix(self.order_u, self.shape[0], 0, u_vec, self.knots_u, 
            self.order_v, self.shape[1], 0, v_vec, self.knots_v, 
            len(u_vec), data, row_indices, col_indices)

        basis0 = sps.csc_matrix((data, (row_indices, col_indices)), shape=(len(u_vec), num_control_points) )
        
        return basis0

    def compute_eval_map_derivs1(self, u_vec, v_vec):
        data = np.zeros(len(u_vec) * self.order_u * self.order_v)
        row_indices = np.zeros(len(data), np.int32)
        col_indices = np.zeros(len(data), np.int32)        

        num_control_points = self.shape[0] * self.shape[1]

        get_basis_surface_matrix(self.order_u, self.control_points.shape[0], 1, u_vec, self.knots_u, 
            self.order_v, self.control_points.shape[1], 1, v_vec, self.knots_v, 
            len(u_vec), data, row_indices, col_indices)

        basis1 = sps.csc_matrix((data, (row_indices, col_indices)), shape=(len(u_vec), num_control_points) )
        
        return basis1

    def compute_eval_map_derivs2(self, u_vec, v_vec):
        data = np.zeros(len(u_vec) * self.order_u * self.order_v)
        row_indices = np.zeros(len(data), np.int32)
        col_indices = np.zeros(len(data), np.int32)      

        num_control_points = self.shape[0] * self.shape[1]
        
        get_basis_surface_matrix(self.order_u, self.control_points.shape[0], 2, u_vec, self.knots_u, 
            self.order_v, self.control_points.shape[1], 2, v_vec, self.knots_v, 
            len(u_vec), data, row_indices, col_indices)

        basis2 = sps.csc_matrix((data, (row_indices, col_indices)), shape=(len(u_vec), num_control_points) )
        
        return basis2

    def evaluate_points(self, u_vec, v_vec):       
        num_control_points = self.shape[0] * self.shape[1]
        
        basis0 = self.compute_eval_map_points(u_vec, v_vec)
        points = basis0.dot(self.control_points.reshape((num_control_points, 3)))

        return points

    def evaluate_der1(self, u_vec, v_vec):
        num_control_points = self.shape[0] * self.shape[1]
        
        basis1 = self.compute_eval_map_der1(u_vec, v_vec)
        derivs1 = basis1.dot(self.control_points.reshape((num_control_points, 3)))

        return derivs1 

    def evaluate_der2(self, u_vec, v_vec):
        num_control_points = self.shape[0] * self.shape[1]
        
        basis2 = self.compute_eval_map_der2(u_vec, v_vec)
        derivs2 = basis2.dot(self.control_points.reshape((num_control_points, 3)))

        return derivs2


    def project(self, points_to_project, max_iter=100):

        num_points = len(points_to_project)

        u_vec = np.zeros(num_points)
        v_vec = np.zeros(num_points)

        compute_surface_projection(
            self.order_u, self.control_points.shape[0],
            self.order_v, self.control_points.shape[1],
            num_points, max_iter,
            points_to_project, 
            self.control_points,
            u_vec, v_vec, 50
        )

        return u_vec, v_vec

    def compute_projection_eval_map(self, points_to_project, max_iter=100):
        u_vec, v_vec = self.project(points_to_project, max_iter)

        basis0 = self.compute_eval_map(u_vec, v_vec)
        
        return basis0
