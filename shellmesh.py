import numpy as np
import scipy.sparse as sps

import vedo
import triangle as tr
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'/Users/Sansara/Public/Code/Geomesh/ShellMesh/lsdo_kit')
from lsdo_kit.simulation.mesh.mesh import Mesh
from lsdo_kit.design.design_geometry.core.pointset import PointSet
from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
from lsdo_kit.cython.surface_projection_py import compute_surface_projection
from lsdo_kit.design.design_geometry.core.pointset import ProjectedPointSet
from lsdo_kit.design.design_geometry.bsplines.bspline_surface import BSplineSurface
from lsdo_kit.design.design_geometry.core.pointset_functions import add_pointsets, subtract_pointsets, divide_pointset_by_scalar
from lsdo_kit.cython.get_open_uniform_py import get_open_uniform

class ShellMesh(Mesh):
    def __init__(self, name, pointset_list=[]) -> None:
        super().__init__(name, pointset_list)
        self.pointest_dict = {}
        for pointset in pointset_list:
            self.pointest_dict[pointset.name] = pointset
        #self.mapping = np.array([])
        self.con_indices = 0
        self.connectivity = np.array([],dtype=np.int32)
        self.triangulation_dict = {}

    def extract_pointset_list_from_bspline_surface(self, geo, bspline_surface_list): 
        pointset_list = [] 
        for bspline_surface in bspline_surface_list:
            pointset = self.extract_pointset_from_bspline_surface(geo, bspline_surface)
            pointset_list.append(pointset)
        return pointset_list
    
    def extract_pointset_from_bspline_surface(self, geo, bspline_surface): #, num_points_u = 30, num_points_v = 18
        order_u = bspline_surface.order_u
        order_v = bspline_surface.order_v
        num_control_points_u = bspline_surface.shape[0]
        num_control_points_v = bspline_surface.shape[1]
        num_points_u = bspline_surface.shape[0] + 1 
        num_points_v = bspline_surface.shape[0] + 1 
        num_points = num_points_u * num_points_v
        nnz = num_points * order_u * order_v
        data = np.zeros(nnz)
        row_indices = np.zeros(nnz, np.int32)
        col_indices = np.zeros(nnz, np.int32)
        knot_vector_u = bspline_surface.knots_u
        knot_vector_v = bspline_surface.knots_v
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

        relative_map = np.zeros((num_points,len(geo.total_cntrl_pts_vector)))
        linear_map = basis0.dot(np.identity(num_control_points_u * num_control_points_v))
        j = 0
        for surf in geo.input_bspline_entity_dict.values():
            if surf == bspline_surface:
                relative_map[:, j:j+surf.shape[0]*surf.shape[1]] = linear_map
            j = j + surf.shape[0]*surf.shape[1] 
        relative_map = sps.csc_matrix(relative_map)
        pointset = PointSet(
            pointset_id = geo.current_id,
            shape = np.array([num_points_u*num_points_v, 3]),
            output_starting_ind = geo.current_pointset_pts,
            parent_pointers = [],
            absolute_map = None,
            relative_map = relative_map,
            offset = np.array([0., 0., 0.]),            
            offset_absolute_map = None,
            physical_coordinates = None,
            permutation_matrix = None,
            name = 'OML: ' + bspline_surface.name
            )
        geo.pointsets_dict[geo.current_id] = pointset
        geo.current_id += 1
        geo.current_pointset_pts += np.cumprod(pointset.shape)[-2]
        return pointset

    def merge_OML(self, geo, merged_OML_relationship, plot = False): #bspline_entity_list, 
        OML_ctrl_pointset_list = []
        for merged_list in merged_OML_relationship:
            pointset0 = merged_list[0]
            pointset1 = merged_list[1]  
            geo.assemble(pointset = pointset0)
            points0 = geo.evaluate(pointset = pointset0) 
            geo.assemble(pointset = pointset1)
            points1 = geo.evaluate(pointset = pointset1) 
            # print('points0',points0.shape)
            # print('points1',points1.shape)
            num_points0 = pointset0.shape[0]
            num_points1 = pointset1.shape[0]
            # print('num_points0', num_points0)
            # print('num_points1', num_points1)
            relative_map = sps.vstack([pointset0.absolute_map,pointset1.absolute_map])
            #print('relative_map', relative_map.shape)
            output_pointset = PointSet(
                pointset_id = geo.current_id,
                shape = np.append(num_points0+num_points1,3),
                output_starting_ind = geo.current_pointset_pts,
                parent_pointers = [],
                absolute_map = None,
                relative_map = relative_map,
                offset = np.array([0., 0., 0.]),            
                offset_absolute_map = None,
                physical_coordinates = None,
                permutation_matrix = None,
                name = None 
                )     
            geo.pointsets_dict[geo.current_id] = output_pointset
            geo.current_id += 1
            geo.current_pointset_pts += np.cumprod(output_pointset.shape)[-2]


            A = np.around(points0, decimals=8)
            B = np.around(points1, decimals=8)
            intersection_bool = (B[:, None] == A).all(-1).any(1)
            intersection_points = points1[intersection_bool]
            #print('intersection_points', intersection_points.shape)
            points1_reduced_indices = np.where(np.invert((intersection_bool)))[0]
            #print('points1_reduced_indices', points1_reduced_indices)

            point_indices = np.append(np.arange(num_points0), points1_reduced_indices+num_points1)#
            #print('point_indices', len(point_indices), point_indices)
            #print('output_pointset',output_pointset.shape)
            output_pointset = geo.extract_pointset(output_pointset, point_indices, len(point_indices))#[1062    3]
            #print('output_pointset',output_pointset.shape)

            if len(point_indices)%len(intersection_points) ==0:
                output_pointset.shape = np.array([len(point_indices)//len(intersection_points), len(intersection_points), 3]) 
            else: 
                print('Warning') 
            geo.assemble(pointset = output_pointset)
            geo.evaluate(pointset = output_pointset)
            merged_ctrl_pointsets = geo.fit_bspline_ctrl_pointsets([output_pointset])
            merged_ctrl_pointsets = merged_ctrl_pointsets[0]
            merged_ctrl_pointsets.name = 'ctrl_pts_' + merged_list[2]  
            
            OML_ctrl_pointset_list.append(merged_ctrl_pointsets)
            self.pointest_dict[merged_ctrl_pointsets.name] = merged_ctrl_pointsets

            if plot:
                vd_points0 = vedo.Points(points0, r=10, c='red',alpha=0.8)
                vd_points1 = vedo.Points(points1, r=15, c='green',alpha=0.5) 
                vd_test = vedo.Plotter(axes=1)
                vd_test.show(vd_points0, vd_points1,  'Test', viewup="z", interactive=True) 
        return OML_ctrl_pointset_list

    def identify_intersection(self, geo, intersection_list):
        ''' You can assume that all intersections between structural features are between 
        an edge of one feature and a surface of another feature, 
        rather than a general surface-surface intersection. 
        To use an analogy, we can assume all are T intersections, rather than '+' intersections'''
        pointset0 = geo.pointsets_dict[intersection_list[0][0]]
        edge = np.empty((0,2),dtype = int)
        print(edge)
        # exit()
        print(edge.shape)
        vd_points4 = []
        num_points_u_fitted = 40
        num_points_v_fitted = 20        
        surf_u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u_fitted), np.ones(num_points_v_fitted)).flatten()
        surf_v_vec = np.einsum('i,j->ij', np.ones(num_points_u_fitted), np.linspace(0., 1., num_points_v_fitted)).flatten()
        for intersection in intersection_list:
            if intersection[2] == '-':
                pass            
            if intersection[2] == 'T':
                print('T')            
                pointset1 = geo.pointsets_dict[intersection[1]]
                geo.assemble(pointset = pointset0)
                _ = geo.evaluate(pointset = pointset0)
                geo.assemble(pointset = pointset1)
                _ = geo.evaluate(pointset = pointset1)
                if intersection[1] == 64:
                    points = self.discritize_bspline_surface(pointset1, num_points_u = 35, num_points_v = 5)
                    A = dict(vertices=points)#
                    B = tr.triangulate(A,'pc')
                    tr.compare(plt, A, B)
                else:
                    points = self.discritize_bspline_surface(pointset1)
                num_points = len(points)
                #print(num_points)
                order_u = 4
                order_v = 4
                num_control_points_u = pointset0.shape[0]
                num_control_points_v = pointset0.shape[1]
                num_points_u = 40
                num_points_v = 20

                nnz = num_points_u * num_points_v * order_u * order_v
                data = np.zeros(nnz)
                knot_vector_u = np.zeros(num_control_points_u +order_u)
                knot_vector_v = np.zeros(num_control_points_v +order_v)
                get_open_uniform(order_u, num_control_points_u, knot_vector_u)
                get_open_uniform(order_v, num_control_points_v, knot_vector_v)
                cps = pointset0.physical_coordinates
                projection_direction=np.array([0., 0., 0.])
                projection_direction = np.array((projection_direction),ndmin=2)
                axis = np.repeat(projection_direction,num_points,axis = 0)
                surfs_cp = cps.reshape((1,num_control_points_u * num_control_points_v * 3))

                max_iter = 500
                u_vec = np.ones(num_points)
                v_vec = np.ones(num_points)
                surfs_index = np.zeros(num_points,dtype=int)
                surfs_u_order = np.empty((0,1),dtype = int)
                surfs_v_order = np.empty((0,1),dtype = int)
                surfs_num_control_points_u = np.empty((0,1),dtype = int)
                surfs_num_control_points_v = np.empty((0,1),dtype = int)
                surfs_cp = np.empty((0,num_control_points_u * num_control_points_v*3),dtype = int)
                for i in range(1): 
                    surfs_u_order = np.append(surfs_u_order,np.array((order_u),ndmin=2),axis = 0)
                    surfs_v_order = np.append(surfs_v_order,np.array((order_v),ndmin=2),axis = 0)
                    surfs_num_control_points_u = np.append(surfs_num_control_points_u,np.array(num_control_points_u, ndmin=2),axis = 0)
                    surfs_num_control_points_v = np.append(surfs_num_control_points_v,np.array(num_control_points_v, ndmin=2),axis = 0)
                    surfs_cp = np.append(surfs_cp,cps.reshape((1,num_control_points_u * num_control_points_v*3)),axis = 0)
                compute_surface_projection(
                    surfs_u_order.reshape((len(range(1)))),  surfs_num_control_points_u.reshape((len(range(1)))),#np.array((num_control_points_u,)).reshape((1,1)),
                    surfs_v_order.reshape((len(range(1)))),  surfs_num_control_points_v.reshape((len(range(1)))),#np.array((num_control_points_v,)).reshape((1,1)),
                    num_points, max_iter,
                    points.reshape(num_points * 3), 
                    cps.reshape((num_control_points_u * num_control_points_v * 3)),
                    knot_vector_u, knot_vector_v,
                    u_vec, v_vec, 50,
                    axis.reshape(num_points * 3),
                    surfs_index,
                    surfs_cp,
                )

                row_indices = np.zeros(nnz, np.int32)
                col_indices = np.zeros(nnz, np.int32)
                get_basis_surface_matrix(
                    order_u, num_control_points_u, 0, u_vec, knot_vector_u, 
                    order_v, num_control_points_v, 0, v_vec, knot_vector_v,
                    num_points, data, row_indices, col_indices,
                )
                basis0 = sps.csc_matrix(
                    (data, (row_indices, col_indices)), 
                    shape=(num_points, num_control_points_u * num_control_points_v),
                )           
                pts = basis0.dot(cps)
                vd_points4.append(vedo.Points(pts, r=15, c='red',alpha=0.5))
                print('u_vec',u_vec)
                print('v_vec',v_vec)

                num_points_u_fitted = num_points_u
                num_points_v_fitted = num_points_v 
                num_points_fitted = num_points_u_fitted*num_points_v_fitted          
                uv_vec = np.append(u_vec.reshape((num_points,1)), v_vec.reshape((num_points,1)), axis=1)
                surf_uv_vec = np.append(surf_u_vec.reshape((num_points_fitted,1)), surf_v_vec.reshape((num_points_fitted,1)), axis=1)
                # print(surf_uv_vec.shape)
                # print(np.shape(surf_uv_vec - uv_vec[0,:]))
                # print(len(np.linalg.norm(surf_uv_vec - uv_vec[0,:], axis = 1)))
                indices = []
                #print()
                for i in range(num_points):
                    dist = np.linalg.norm(surf_uv_vec - uv_vec[i,:], axis = 1)
                    #print(dist)
                    index = np.argmin(dist)
                    #print(index)
                    indices.append(index)
                    surf_uv_vec[index,:] = 1000
                    #index()
                
                #print(indices)
                for i in range(num_points):
                    #print(indices[i], )
                    surf_u_vec[indices[i]] = uv_vec[i,0]
                    surf_v_vec[indices[i]] = uv_vec[i,1]
                    surf_uv_vec[indices[i],:] = uv_vec[i,:]
                    if i != num_points - 1:
                        #print(np.array([[indices[i], indices[i+1]]]).shape)
                        edge = np.append(edge, np.array([[indices[i], indices[i+1]]]).reshape(1,2), axis = 0)                    
                # print(edge.shape)
                # print(surf_uv_vec.shape)
                # print(indices)
                # print(edge)
                # print()

                #plt.show()
               
                #surf_u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u_fitted), np.ones(num_points_v_fitted)).flatten()
                #surf_v_vec = np.einsum('i,j->ij', np.ones(num_points_u_fitted), np.linspace(0., 1., num_points_v_fitted)).flatten()
                order_u_fitted = 4
                order_v_fitted = 4
                num_control_points_u_fitted = pointset0.shape[0]
                num_control_points_v_fitted = pointset0.shape[1]
                nnz = num_points_u_fitted * num_points_v_fitted * order_u_fitted * order_v_fitted
                data = np.zeros(nnz)
                row_indices = np.zeros(nnz, np.int32)
                col_indices = np.zeros(nnz, np.int32)
                knot_vector_u = np.zeros(num_control_points_u_fitted+order_u_fitted)
                knot_vector_v = np.zeros(num_control_points_v_fitted+order_v_fitted)
                get_open_uniform(order_u_fitted, num_control_points_u_fitted, knot_vector_u)
                get_open_uniform(order_v_fitted, num_control_points_v_fitted, knot_vector_v)
                get_basis_surface_matrix(
                    order_u_fitted, num_control_points_u_fitted, 0, surf_u_vec, knot_vector_u,
                    order_v_fitted, num_control_points_v_fitted, 0, surf_v_vec, knot_vector_v,
                    num_points_u_fitted * num_points_v_fitted, data, row_indices, col_indices,
                        )
                basis0 = sps.csc_matrix(
                    (data, (row_indices, col_indices)), 
                    shape=(num_points_u_fitted * num_points_v_fitted, num_control_points_u_fitted * num_control_points_v_fitted),
                )
                points = basis0.dot(pointset0.physical_coordinates)
                print(points.shape)
            if intersection[2] == '+':
                print("Sorry! This has not been implemented yet!")

        A = dict(vertices=surf_uv_vec, segments=edge)#
        B = tr.triangulate(A,'pc')
        tr.compare(plt, A, B)
        #plt.show()

        return points, pts, B, vd_points4 #connectivity

    def discritize_bspline_surface(self, ctrl_pointset, num_points_u = 10, num_points_v = 8, direction = 'u'):
        order_u_fitted = 4
        order_v_fitted = 4
        num_control_points_u_fitted = ctrl_pointset.shape[0]
        num_control_points_v_fitted = ctrl_pointset.shape[1]
        num_points_u_fitted = num_points_u
        num_points_v_fitted = num_points_v
        nnz = num_points_u_fitted * num_points_v_fitted * order_u_fitted * order_v_fitted
        data = np.zeros(nnz)
        row_indices = np.zeros(nnz, np.int32)
        col_indices = np.zeros(nnz, np.int32)
        u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u_fitted), np.ones(num_points_v_fitted)).flatten()
        v_vec = np.einsum('i,j->ij', np.ones(num_points_u_fitted), np.linspace(0., 1., num_points_v_fitted)).flatten()
        knot_vector_u = np.zeros(num_control_points_u_fitted+order_u_fitted)
        knot_vector_v = np.zeros(num_control_points_v_fitted+order_v_fitted)
        get_open_uniform(order_u_fitted, num_control_points_u_fitted, knot_vector_u)
        get_open_uniform(order_v_fitted, num_control_points_v_fitted, knot_vector_v)
        get_basis_surface_matrix(
            order_u_fitted, num_control_points_u_fitted, 0, u_vec, knot_vector_u,
            order_v_fitted, num_control_points_v_fitted, 0, v_vec, knot_vector_v,
            num_points_u_fitted * num_points_v_fitted, data, row_indices, col_indices,
                )
        #pts = basis0.dot(bspline_entity.control_points) 
        basis0 = sps.csc_matrix(
            (data, (row_indices, col_indices)), 
            shape=(num_points_u_fitted * num_points_v_fitted, num_control_points_u_fitted * num_control_points_v_fitted),
        )
        #print(len(u_vec)) 
        # print('u_vec', u_vec)
        # print('v_vec', v_vec)
        points = basis0.dot(ctrl_pointset.physical_coordinates)
        #print(np.shape(points))
        if direction == 'u':
            indices = num_points_v* np.arange(num_points_u)
            points = points[indices, :]
        return points


    
    def create_triangulation(self):
        pass

if __name__ == '__main__':
    from shellmesh import ShellMesh
    shell_mesh = ShellMesh('shell_mesh')

'''    def mapping_from_OML(self): 
        relative_map = np.identity(len(self.total_cntrl_pts_vector))
        relative_map = sps.csc_matrix(relative_map)
        pointset = PointSet(
            pointset_id = self.current_id,
            shape = np.append(len(self.total_cntrl_pts_vector),3),
            output_starting_ind = self.current_pointset_pts,
            parent_pointers = [],
            absolute_map = None,
            relative_map = relative_map,
            offset = np.array([0., 0., 0.]),            
            offset_absolute_map = None,
            physical_coordinates = None,
            permutation_matrix = None,
            name = 'OML'
            )
        self.pointsets_dict[self.current_id] = pointset
        self.current_id += 1
        self.current_pointset_pts += np.cumprod(pointset.shape)[-2]
        
        self.assemble(pointset = pointset)
        points = self.evaluate(pointset = pointset)
        color = list(colors.colors.values())[self.current_id]
        vp_points = vedo.Points(points, r=5, c=color,alpha=0.3)
        vp_test = Plotter(axes=1)
        vp_test.show(vp_points, 'Test', viewup="z", interactive=True)
        pass

            
        ##pts = basis0.dot(bspline_entity.control_points) 
        # geo.assemble(pointset = pointset)
        # points = geo.evaluate(pointset = pointset)
        # color = list(colors.colors.values())[self.current_id]
        # vp_points = vedo.Points(points, r=5, c=color,alpha=0.3)
        # vp_test = Plotter(axes=1)
        # vp_test.show(vp_points, 'Test', viewup="z", interactive=False)

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
'''