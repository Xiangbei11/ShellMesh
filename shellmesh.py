import numpy as np
import scipy.sparse as sps

import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit')
from lsdo_kit.simulation.mesh.mesh import Mesh
from lsdo_kit.design.design_geometry.core.pointset import ProjectedPointSet
from lsdo_kit.design.design_geometry.bsplines.bspline_surface import BSplineSurface
from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix



class ShellMesh(Mesh):
    def __init__(self, name, pointset_list=[]) -> None:
        super().__init__(name, pointset_list)
    
    def add_design_geometry(self, design_geometry_class = None):
        self.geo = design_geometry_class #TODO ShellMesh class should not include design_geometry_class?
        bspline_entity_list = []
        for bspline_entity_surface in self.geo.input_bspline_entity_dict.values():
            pointset = self.extract_pointset_from_bspline_entity(self.geo, bspline_entity_surface)
        
        
    def extract_pointset_from_bspline_entity(self, geo, bspline_entity): #TODO Move to lsdo_kit/lsdo_kit/design/design_geometry/core/pointset_functions ?
        order_u = bspline_entity.order_u
        order_v = bspline_entity.order_v
        num_control_points_u = bspline_entity.shape[0]
        num_control_points_v = bspline_entity.shape[1]
        num_points_u = 50   # TODO might want to pass these in as input
        num_points_v = 50
        num_points = num_points_u * num_points_v
        nnz = num_points * order_u * order_v
        data = np.zeros(nnz)
        row_indices = np.zeros(nnz, np.int32)
        col_indices = np.zeros(nnz, np.int32)

        knot_vector_u = bspline_entity.knots_u
        knot_vector_v = bspline_entity.knots_v
        u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u), np.ones(num_points_v)).flatten()
        v_vec = np.einsum('i,j->ij', np.ones(num_points_u), np.linspace(0., 1., num_points_v)).flatten()

        get_basis_surface_matrix(
            order_u, num_control_points_u, 0, u_vec, knot_vector_u,
            order_v, num_control_points_v, 0, v_vec, knot_vector_v,
            num_points, data, row_indices, col_indices,
        )

        basis0 = sps.csc_matrix(
            (data, (row_indices, col_indices)),
            shape=(num_points, num_control_points_u * num_control_points_v),
        )
        #pts = basis0.dot(bspline_entity.control_points) 

        relative_map = np.zeros((num_points,len(geo.total_cntrl_pts_vector)))
        linear_map = basis0.dot(np.identity(bspline_entity.shape[0] * bspline_entity.shape[1]))
        for surf in geo.input_bspline_entity_dict.values():
            j = 0
            if surf == bspline_entity:
                relative_map[:, j:j+surf.shape[0]*surf.shape[1]] = linear_map
            j = j + surf.shape[0]*surf.shape[1] 

        pointset = ProjectedPointSet( #TODO name: PointSet or ProjectedPointSet or another?
            pointset_id = geo.current_id,
            shape = np.append(num_points,3),
            output_starting_ind = geo.current_pointset_pts,
            parent_pointers = [],
            absolute_map = None,
            relative_map = relative_map,
            offset = None,            
            offset_absolute_map = None,
            physical_coordinates = None,
            permutation_matrix = None,
            )
        geo.pointsets_dict[geo.current_id] = pointset
        geo.current_id += 1
        geo.current_pointset_pts += np.cumprod(pointset.shape)[-2]
        geo.register_output(pointset)
        return pointset
    
    def identify_intersection(self, intersection_list):
        ''' You can assume that all intersections between structural features are between 
        an edge of one feature and a surface of another feature, 
        rather than a general surface-surface intersection. 
        To use an analogy, we can assume all are T intersections, rather than '+' intersections'''
        pass
        for intersection in intersection_list:
            if intersection[2] == 'T':
                pointset1 = intersection[0]
                parent_pointset = intersection[1]
                pointset2 = self.geo.extract_pointset(self.geo, parent_pointset, point_indices, shape, offset=np.array([0., 0., 0.]))
                self.geo.subtract_pointsets(self.geo, pointset1, pointset2, offset=np.array([0., 0., 0.]))
            if intersection[2] == '+':
                print("Sorry! This has not been implemented yet!")
    
    def create_triangulation(self):
        pass

if __name__ == '__main__':
    from shellmesh import ShellMesh
    shell_mesh = ShellMesh('shell_mesh')