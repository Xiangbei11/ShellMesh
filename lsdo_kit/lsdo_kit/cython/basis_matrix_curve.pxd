from libc.stdlib cimport malloc, free

from lsdo_kit.cython.get_open_uniform cimport get_open_uniform
from lsdo_kit.cython.basis0 cimport get_basis0
from lsdo_kit.cython.basis1 cimport get_basis1
from lsdo_kit.cython.basis2 cimport get_basis2


cdef get_basis_curve_matrix(
    int order, int num_control_points, int u_der, double* u_vec, double* knot_vector,
    int num_points, double* data, int* row_indices, int* col_indices,
)