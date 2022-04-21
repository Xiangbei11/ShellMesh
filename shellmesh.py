import numpy as np
import scipy.sparse as sps

import vedo
import triangle as tr

import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit')
from lsdo_kit.simulation.mesh.mesh import Mesh
from lsdo_kit.design.design_geometry.core.pointset import PointSet, ProjectedPointSet, DerivedPointSet, EvaluatedPointSets
from lsdo_kit.design.design_geometry.core.pointset import ProjectedPointSet
from lsdo_kit.design.design_geometry.bsplines.bspline_surface import BSplineSurface
from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
from lsdo_kit.design.design_geometry.core.pointset_functions import add_pointsets, subtract_pointsets, divide_pointset_by_scalar

class ShellMesh(Mesh):
    def __init__(self, name, pointset_list=[]) -> None:
        super().__init__(name, pointset_list)
        self.mapping = np.array([])
        self.con_indices = 0
        self.connectivity = np.array([],dtype=np.int32)
        self.triangulation_dict = {}

    def merge_OML(self, geo, merged_OML_list): 
        output_pointset = []
        for merged_list in merged_OML_list:
            pointset0 = merged_list[0]
            pointset1 = merged_list[1]  
            num_pts0 = pointset0.relative_map.shape[0]
            num_pts1 = pointset1.relative_map.shape[0]
            print('num_pts0',num_pts0)
            print('num_pts1',num_pts1)
            self.mapping = sps.vstack([pointset0.absolute_map,pointset1.absolute_map])
            print('self.mapping', self.mapping.shape)
            output_pointset = PointSet(
                pointset_id = geo.current_id,
                shape = np.append(num_pts0+num_pts1,3),
                output_starting_ind = geo.current_pointset_pts,
                parent_pointers = [],
                absolute_map = None,
                relative_map = self.mapping,
                offset = np.array([0., 0., 0.]),            
                offset_absolute_map = None,
                physical_coordinates = None,
                permutation_matrix = None,
                name = pointset0.name + '-' + pointset1.name
                )     
            geo.pointsets_dict[geo.current_id] = output_pointset
            geo.current_id += 1
            geo.current_pointset_pts += np.cumprod(output_pointset.shape)[-2]

            geo.assemble(pointset = pointset0)
            points0 = geo.evaluate(pointset = pointset0) 
            geo.assemble(pointset = pointset1)
            points1 = geo.evaluate(pointset = pointset1) 
            print('points0',points0.shape)
            print('points1',points1.shape)
            A = np.around(points0, decimals=8)
            B = np.around(points1, decimals=8)
            intersection_bool = (A[:, None] == B).all(-1).any(1)
            #print(np.where(np.invert((intersection_bool))))
            print(np.where((intersection_bool)))
            intersection_points = points0[intersection_bool]
            print('intersection', intersection_points.shape)


            points0_reduced_indices = np.where(np.invert((intersection_bool)))[0]
            point_indices = np.append(points0_reduced_indices, np.arange(num_pts1)+num_pts0)#
            # print('point_indices', len(point_indices), point_indices)
            # print('output_pointset', output_pointset.shape)
            output_pointset = geo.extract_pointset(output_pointset, point_indices, len(point_indices))#[1062    3]
            
            geo.assemble(pointset = output_pointset)
            points = geo.evaluate(pointset = output_pointset)
            num_points_u = 50
            num_points_v = 25
            u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u), np.ones(num_points_v)).flatten()
            v_vec = np.einsum('i,j->ij', np.ones(num_points_u), np.linspace(0., 1., num_points_v)).flatten()
            # print(u_vec[:10])
            # print(v_vec[:10]) 

            # vertices1 = np.concatenate((verts,edge,edge2))
            # print(np.shape(vertices1)) 
            # edge = np.array(([[0,3],[17,3],[17,14],[0,14],[2,6],[6,12],[12,16]]), dtype = np.int32)#,[21,22]
            # print(np.shape(edge)) 
            # A = dict(vertices=vertices1, segments=edge)#
            # B = tr.triangulate(A,'pc')
            # #print(B)
            # tr.compare(plt, A, B)
            # A = dict(vertices=vertices1)#, segments=edge
            # B = tr.triangulate(A)
            # tr.compare(plt, A, B)
            # plt.show()

            self.triangulation_dict[output_pointset] = edge
        return output_pointset, intersection_points

    def identify_intersection(self, geo, intersection_list):
        ''' You can assume that all intersections between structural features are between 
        an edge of one feature and a surface of another feature, 
        rather than a general surface-surface intersection. 
        To use an analogy, we can assume all are T intersections, rather than '+' intersections'''
        for intersection in intersection_list:
            if intersection[2] == '-':
                pass            
            if intersection[2] == 'T':
                print('T')
                # pointset1 = intersection[0]
                # parent_pointset = intersection[1]
                # pointset2 = self.geo.extract_pointset(self.geo, parent_pointset, point_indices, shape, offset=np.array([0., 0., 0.]))
                # self.geo.subtract_pointsets(self.geo, pointset1, pointset2, offset=np.array([0., 0., 0.]))
            if intersection[2] == '+':
                print("Sorry! This has not been implemented yet!")

    
    def create_triangulation(self):
        pass

if __name__ == '__main__':
    from shellmesh import ShellMesh
    shell_mesh = ShellMesh('shell_mesh')