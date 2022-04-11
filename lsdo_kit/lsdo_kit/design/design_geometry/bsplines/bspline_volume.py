import numpy as np

import scipy.sparse as sps
from lsdo_kit.cython.basis_matrix_volume_py import get_basis_volume_matrix
from lsdo_kit.cython.volume_projection_py import compute_volume_projection


class BSplineVolume:
    def __init__(self, name, order_u, order_v, order_w, knots_u, knots_v, knots_w, shape, control_points):
        self.name = name
        self.order_u = order_u
        self.knots_u = knots_u
        self.order_v = order_v
        self.knots_v = knots_v
        self.order_w = order_w
        self.knots_w = knots_w
        self.shape = shape
        self.starting_geometry_index = None
        self.control_points = control_points
        self.num_control_points = shape[0] * shape[1] * shape[2]

        # print('shape: ', shape)
        # print('num_control_points: ', self.num_control_points)


    def compute_eval_map_points(self, u_vec, v_vec, w_vec):
        data = np.zeros(len(u_vec) * self.order_u * self.order_v * self.order_w)
        row_indices = np.zeros(len(data), np.int32)
        col_indices = np.zeros(len(data), np.int32)

        # print('len(u_vec): ', len(u_vec))

        # print('order_u: ', self.order_u)
        # print('order_v: ', self.order_v)
        # print('order_w: ', self.order_w)

        # print('u_vec: ', u_vec)
        # print('v_vec: ', v_vec)
        # print('w_vec: ', w_vec)

        # print('data: ', data)
        # print('row_indices: ', row_indices)
        # print('col_indices: ', col_indices)

        get_basis_volume_matrix(self.order_u, self.shape[0], 0, u_vec, self.knots_u, 
            self.order_v, self.shape[1], 0, v_vec, self.knots_v,
            self.order_w, self.shape[2], 0, w_vec, self.knots_w, 
            len(u_vec), data, row_indices, col_indices)

        # print('Row Indices: ', np.where(row_indices<0))
        temp = np.where(col_indices<0)
        # print('Column Indices: ', temp)
        # print('Column Indices Values: ', col_indices[temp])
        basis0 = sps.csc_matrix((data, (row_indices, col_indices)), shape=(len(u_vec), self.num_control_points) )
        
        return basis0


    def compute_eval_map_der1(self, u_vec, v_vec, w_vec):
        data = np.zeros(len(u_vec) * self.order_u * self.order_v * self.order_w)
        row_indices = np.zeros(len(data), np.int32)
        col_indices = np.zeros(len(data), np.int32)

        get_basis_volume_matrix(self.order_u, self.shape[0], 1, u_vec, self.knots_u, 
            self.order_v, self.shape[1], 1, v_vec, self.knots_v,
            self.order_w, self.shape[2], 1, w_vec, self.knots_w, 
            len(u_vec), data, row_indices, col_indices)

        basis1 = sps.csc_matrix((data, (row_indices, col_indices)), shape=(len(u_vec), self.num_control_points) )
        
        return basis1


    def compute_eval_map_der2(self, u_vec, v_vec, w_vec):
        data = np.zeros(len(u_vec) * self.order_u * self.order_v * self.order_w)
        row_indices = np.zeros(len(data), np.int32)
        col_indices = np.zeros(len(data), np.int32)

        get_basis_volume_matrix(self.order_u, self.shape[0], 2, u_vec, self.knots_u, 
            self.order_v, self.shape[1], 2, v_vec, self.knots_v,
            self.order_w, self.shape[2], 2, w_vec, self.knots_w, 
            len(u_vec), data, row_indices, col_indices)

        basis2 = sps.csc_matrix((data, (row_indices, col_indices)), shape=(len(u_vec), self.num_control_points) )
        
        return basis2


    def evaluate_points(self, u_vec, v_vec, w_vec):

        basis0 = self.compute_eval_map_points(u_vec, v_vec, w_vec)
        points = basis0.dot(self.control_points.reshape((self.num_control_points, 3)))

        return points


    def evaluate_der1(self, u_vec, v_vec, w_vec):

        basis1 = self.compute_eval_map_der1(u_vec, v_vec, w_vec)
        derivs1 = basis1.dot(self.control_points.reshape((self.num_control_points, 3)))

        return derivs1

    def evaluate_der2(self, u_vec, v_vec, w_vec):

        basis2 = self.compute_eval_map_der2(u_vec, v_vec, w_vec)
        derivs2 = basis2.dot(self.control_points.reshape((self.num_control_points, 3)))

        return derivs2


    def project(self, points_to_project, max_iter=100):

        num_points = len(points_to_project)

        u_vec = np.zeros(num_points)
        v_vec = np.zeros(num_points)
        w_vec = np.zeros(num_points)

        compute_volume_projection(
            self.order_u, self.shape[0],
            self.order_v, self.shape[1],
            self.order_w, self.shape[2],
            num_points, max_iter,
            points_to_project.reshape(num_points * 3), 
            self.control_points.reshape(self.num_control_points * 3),
            self.knots_u, self.knots_v, self.knots_w,
            u_vec, v_vec, w_vec, 10, np.array([0.,0.,0.])
        )

        return u_vec, v_vec, w_vec

    def compute_projection_eval_map(self, points_to_project, max_iter=100):
        u_vec, v_vec, w_vec = self.project(points_to_project, max_iter)

        basis0 = self.compute_eval_map_points(u_vec, v_vec, w_vec)
        
        return basis0

