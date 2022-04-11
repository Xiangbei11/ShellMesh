from libc.stdlib cimport malloc, free

from lsdo_kit.cython.get_open_uniform cimport get_open_uniform
from lsdo_kit.cython.basis0 cimport get_basis0
from lsdo_kit.cython.basis1 cimport get_basis1
from lsdo_kit.cython.basis2 cimport get_basis2


cdef get_basis_volume_matrix(
    int order_u, int num_control_points_u, int u_der, double* u_vec, double* knot_vector_u,
    int order_v, int num_control_points_v, int v_der, double* v_vec, double* knot_vector_v,
    int order_w, int num_control_points_w, int w_der, double* w_vec, double* knot_vector_w,
    int num_points, double* data, int* row_indices, int* col_indices,
)