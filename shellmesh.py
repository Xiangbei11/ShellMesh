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

from meshopt import meshopt
from member import Member

class ShellMesh(Mesh):
    def __init__(self, name, pointset_list=[]) -> None:
        super().__init__(name, pointset_list)
        self.constrained_node_indices = []
        self.num_of_nodes = 0
        self.connected_pointset_list = []
        self.connectivity_dict = {}
        self.members_dict = {}
        self.tri_connectivity = np.empty((0,3), dtype=np.int32) 
        # self.connectivity = np.array([],dtype=np.int32)
        # self.triangulation_dict = {}
        #self.mapping = np.array([])
        # self.pointest_dict = {}
        # for pointset in pointset_list:
        #     self.pointest_dict[pointset.name] = pointset

    def extract_pointset_list_from_bspline_surface(self, geo, bspline_surface_list): 
        
        self.mapping = sps.csc_matrix(np.empty((0,len(geo.total_cntrl_pts_vector))))
        self.total_cntrl_pts_vector = geo.total_cntrl_pts_vector

        pointset_list = [] 
        for bspline_surface in bspline_surface_list:
            pointset = self.extract_pointset_from_bspline_surface(geo, bspline_surface)
            pointset_list.append(pointset)
        return pointset_list
    
    def extract_pointset_from_bspline_surface(self, geo, bspline_surface):
        order_u = bspline_surface.order_u
        order_v = bspline_surface.order_v
        num_control_points_u = bspline_surface.shape[0]
        num_control_points_v = bspline_surface.shape[1]
        num_points_u = bspline_surface.shape[0] + 5 #Use 5 instead of 1 to keep more geometry 
        num_points_v = bspline_surface.shape[0] + 5       
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

    def merge_OML(self, geo, merged_OML_relationship_list, plot = False):
        OML_ctrl_pointset_list = []
        for merged_list in merged_OML_relationship_list:
            pointset0 = geo.pointsets_dict[merged_list[0]]
            pointset1 = geo.pointsets_dict[merged_list[1]]
            geo.assemble(pointset = pointset0)
            points0 = geo.evaluate(pointset = pointset0) 
            geo.assemble(pointset = pointset1)
            points1 = geo.evaluate(pointset = pointset1) 
            num_points0 = pointset0.shape[0]
            num_points1 = pointset1.shape[0]
            relative_map = sps.vstack([pointset0.absolute_map,pointset1.absolute_map])

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
            if len(intersection_points) != 0:               
                points1_reduced_indices = np.where(np.invert((intersection_bool)))[0]
                point_indices = np.append(np.arange(num_points0), points1_reduced_indices+num_points1)
                output_pointset = geo.extract_pointset(output_pointset, point_indices, len(point_indices))
                if len(point_indices)%len(intersection_points) ==0:
                    output_pointset.shape = np.array([len(point_indices)//len(intersection_points), len(intersection_points), 3]) 
                else: 
                    print('Warning0') 
                geo.assemble(pointset = output_pointset)
                geo.evaluate(pointset = output_pointset)
                merged_ctrl_pointsets = geo.fit_bspline_ctrl_pointsets([output_pointset])
                merged_ctrl_pointsets = merged_ctrl_pointsets[0]
                merged_ctrl_pointsets.name = 'ctrl_pts_' + merged_list[2]  
                OML_ctrl_pointset_list.append(merged_ctrl_pointsets)
            else:
                print('Warning1')
                print(merged_list[2],len(intersection_bool), len(intersection_points))
            
            if plot:
                vd_points0 = vedo.Points(points0, r=10, c='red',alpha=0.8)
                vd_points1 = vedo.Points(points1, r=15, c='green',alpha=0.5) 
                vd_test = vedo.Plotter(axes=0)
                vd_test.show(vd_points0, vd_points1,  'Test', viewup="z", interactive=True) 
        return OML_ctrl_pointset_list

    def identify_intersection_list(self, geo, intersection_list, plot = False):
        ''' Assume that all intersections between structural features are between 
        an edge of one feature and a surface of another feature, 
        rather than a general surface-surface intersection.(All T intersections, rather than '+' intersections)'''
        
        pointset_ini = geo.pointsets_dict[intersection_list[0][0]]
        geo.assemble(pointset = pointset_ini)
        _ = geo.evaluate(pointset = pointset_ini)
        num_points_u0, num_points_v0  = intersection_list[0][2], intersection_list[0][3]       
        u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
        v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
        u0_v0_vec = np.vstack((u0_vec, v0_vec)).T

        constrained_edges = np.empty((0,2),dtype = int)
        # constrained_points = np.empty((0,3)) 

        self.pointset_ini_check = 0
        if pointset_ini not in self.connected_pointset_list:
            self.connected_pointset_list.append(pointset_ini)
            self.num_of_nodes_ini = self.num_of_nodes
            self.num_of_nodes += len(u0_v0_vec)
            self.members_dict[pointset_ini.name] = Member(    
                name = pointset_ini.name,
                node_indices = [*range(self.num_of_nodes_ini, self.num_of_nodes)],
                )     
            
            _, basis0 = self.discritize_ctrl_pointset(pointset_ini, uv_vec = u0_v0_vec, num_points_u = num_points_u0, num_points_v = num_points_v0)
            self.mapping = sps.vstack([self.mapping, basis0.dot(pointset_ini.absolute_map)])
            self.mapping = sps.csc_matrix(self.mapping)
            self.pointset_ini_check = 1
        else:
            self.num_of_nodes_ini =  self.members_dict[pointset_ini.name].options['node_indices'][0]

        for intersection in intersection_list:
            pointset0 = geo.pointsets_dict[intersection[0]]
            if pointset0 != pointset_ini:
                print('Warning2')         
            pointset1 = geo.pointsets_dict[intersection[1]]
            num_points_u1, num_points_v1 = intersection[4], intersection[5] 
            geo.assemble(pointset = pointset1)
            _ = geo.evaluate(pointset = pointset1)  
            if self.pointset_ini_check:  
                u0_v0_vec, edge, constrained_node_indices = self.identify_intersection(geo, u0_v0_vec, pointset0, pointset1, num_points_u1,num_points_v1, edge_location = intersection[6], relationship = intersection[7])
            
            #if 'spar' in pointset1.name:
                constrained_edges = np.append(constrained_edges, edge, axis = 0)

                self.constrained_node_indices += constrained_node_indices
            else:
                u0_v0_vec = self.identify_intersection(geo, u0_v0_vec, pointset0, pointset1, num_points_u1,num_points_v1, edge_location = intersection[6], relationship = intersection[7])
            # constrained_points = np.append(constrained_points, points, axis = 0)
            # print('edge',edge.shape)
            # print('constrained_points', constrained_points.shape)
        
        if self.pointset_ini_check:
            _, basis0 = self.discritize_ctrl_pointset(pointset_ini, uv_vec = u0_v0_vec, num_points_u = num_points_u0, num_points_v = num_points_v0)
            self.mapping[self.num_of_nodes_ini:self.num_of_nodes_ini+len(u0_v0_vec),:] = sps.csc_matrix(basis0.dot(pointset_ini.absolute_map))
            
            #constrained_edges = np.delete(constrained_edges, indices, 0) constrained_edges[:-1,:]

            A = dict(vertices=u0_v0_vec, segments=constrained_edges)#
            B = tr.triangulate(A,'pc')#,'c'

            # tr.compare(plt, A, B)
            # plt.show(block = False)
            # print('len(constrained_edges)',len(constrained_edges))
            # print('len(B[vertices])', len(B['vertices']))
            # print('len(u0_v0_vec)',len(u0_v0_vec))
            # if len(u0_v0_vec) != len(B['vertices']):
            #     print('Warning4', len(u0_v0_vec), len(B['vertices']))
            #     check_conn = B['triangles']
            #     indices = []
            #     for i in range(len(check_conn)):
            #         if check_conn[i,0]>=2800 or check_conn[i,1]>=2800 or check_conn[i,2]>=2800:
            #             print(check_conn[i,:])
            #             indices.append(i)
            #     check_conn = np.delete(check_conn, indices, 0)
            #     B['triangles'] = check_conn          

            B['triangles'] += self.num_of_nodes_ini
            self.connectivity_dict[pointset_ini.name] = B['triangles']

            self.members_dict[pointset_ini.name].options['tri_connectivity'] = B['triangles']
            self.members_dict[pointset_ini.name].options['constrained_node_indices'] =  [x + self.num_of_nodes_ini for x in constrained_node_indices]
            self.constrained_node_indices += self.members_dict[pointset_ini.name].options['constrained_node_indices']
        else:
            print('Sorry')                                                            

        if plot:
            test_points = self.mapping.dot(self.total_cntrl_pts_vector)
            print('test_points', test_points.shape)
            #tr.compare(plt, A, B)
            #print(B['segments'])
            #plt.show(block = False)           
            #vd_points0 = vedo.Points(points, r=10, c='black',alpha=1.0)
            #vd_points1 = vedo.Points(constrained_points, r=15, c='red',alpha=0.5)

            #,rib_557,rib_556, rib_517, rib_477,
            # rib_557 = vedo.Points(test_points[557,:].reshape((1,3)), c = 'grey', r=30)
            # rib_556 = vedo.Points(test_points[556,:].reshape((1,3)), c = 'black', r=30)
            # rib_517 = vedo.Points(test_points[517,:].reshape((1,3)), c = 'blue', r=30)
            # rib_477 = vedo.Points(test_points[477,:].reshape((1,3)), c = 'green', r=30)
            # print(constrained_node_indices)
            # print(edge)
            vd_points1 = vedo.Points(test_points[self.constrained_node_indices,:], r=15, c='red',alpha=0.5)
            i = 0
            for conn in self.connectivity_dict.values():
                self.tri_connectivity = np.append(self.tri_connectivity, conn, axis = 0) 
                i+=1
            temp0 = []
            temp1 = []
            for i in range(12):
                temp0_p = test_points[3016+45*i:3016+3+45*i,:]
                temp0.append(vedo.Points(temp0_p, c = 'red', r=15))
                for j in range(3):
                    dist = np.linalg.norm(test_points[2680:2680+168,:] - temp0_p[j,:], axis = 1)
                    index = np.argmin(dist)
                    test_points[2680+index,:] = temp0_p[j,:]           
                if i == 0:
                    temp1_p = test_points[3061-3+45*i:3061+45*i,:]
                    temp1.append(vedo.Points(temp1_p, c = 'red', r=15,alpha=0.5))
                    test_points[2848+12:2848+15] = temp1_p
                elif i==1:
                    temp1_p = test_points[3061-3+45*i:3061+45*i,:]
                    temp1.append(vedo.Points(temp1_p, c = 'red', r=15,alpha=0.5))
                    test_points[2848+24:2848+27] = temp1_p
                elif i==2 or i ==3:
                    temp1_p = test_points[3061-3+45*i:3061+45*i,:]
                    temp1.append(vedo.Points(temp1_p, c = 'red', r=15,alpha=0.5))
                elif i == 4:
                    temp2_p = test_points[self.constrained_node_indices[186],:].reshape((1,3))
                    test_points[3061-1+45*i:3061+45*i,0] = temp2_p[0,0]-0.3
                    test_points[3061-2+45*i:3061+45*i-1,0] = temp2_p[0,0]-0.2
                    test_points[3061-3+45*i:3061+45*i-2,0] = temp2_p[0,0]-0.1
                    test_points[3061-3+45*i:3061+45*i,1] = temp2_p[0,1]

                    temp1_p = test_points[3061-3+45*i:3061+45*i,:]
                    temp1.append(vedo.Points(temp1_p, c = 'red', r=15,alpha=0.5))
                    temp2 = vedo.Points(test_points[2848+63:2848+66].reshape((3,3)), c = 'black', r=20)
                    test_points[2848+63:2848+66] = temp1_p
                elif i == 5:
                    temp2_p = test_points[self.constrained_node_indices[201],:].reshape((1,3))
                    #temp3 = vedo.Points(temp2_p, c = 'blue', r=20)
                    test_points[3061-1+45*i:3061+45*i,0] = temp2_p[0,0]-0.5
                    test_points[3061-2+45*i:3061+45*i-1,0] = temp2_p[0,0]-0.4
                    test_points[3061-3+45*i:3061+45*i-2,0] = temp2_p[0,0]-0.3
                    test_points[3061-3+45*i:3061+45*i,1] = temp2_p[0,1]

                    temp1_p = test_points[3061-3+45*i:3061+45*i,:]
                    temp1.append(vedo.Points(temp1_p, c = 'red', r=15,alpha=0.5))
                    temp2 = vedo.Points(test_points[2848+63+15:2848+66+15].reshape((3,3)), c = 'black', r=20)
                    test_points[2848+63+15:2848+66+15] = temp1_p
                elif i == 6:
                    temp2_p = test_points[self.constrained_node_indices[201+15],:].reshape((1,3))
                    temp3 = vedo.Points(temp2_p, c = 'blue', r=20)
                    test_points[3061-1+45*i:3061+45*i,0] = temp2_p[0,0]-0.3
                    test_points[3061-2+45*i:3061+45*i-1,0] = temp2_p[0,0]-0.2
                    test_points[3061-3+45*i:3061+45*i-2,0] = temp2_p[0,0]-0.1
                    test_points[3061-3+45*i:3061+45*i,1] = temp2_p[0,1]

                    temp1_p = test_points[3061-3+45*i:3061+45*i,:]
                    temp1.append(vedo.Points(temp1_p, c = 'red', r=15,alpha=0.5))
                    temp2 = vedo.Points(test_points[2848+63+15+12:2848+66+15+12].reshape((3,3)), c = 'black', r=20)
                    test_points[2848+63+15+12:2848+66+15+12] = temp1_p
                else:
                    temp1_p = test_points[3061-3+45*i:3061+45*i,:]
                    temp1.append(vedo.Points(temp1_p, c = 'red', r=15,alpha=0.5))
                    temp2 = vedo.Points(test_points[2848+63+15+24+12*(i-7):2848+66+15+24+12*(i-7)].reshape((3,3)), c = 'black', r=20)
                    test_points[2848+63+15+24+12*(i-7):2848+66+15+24+12*(i-7)] = temp1_p
                # temp1_p = test_points[3061-3+45*i:3061+45*i,:]
                # temp1.append(vedo.Points(temp1_p, c = 'black', r=20))
                # indices =[]
                # for j in range(3):
                #     dist = np.linalg.norm(test_points[2848:2848+168,:] - temp1_p[j,:], axis = 1)
                #     index = np.argmin(dist)
                #     if index in indices:
                #         temp = dist[index] 
                #         dist[index] = 1000
                #         index_new = np.argmin(dist)
                #         dist[index] = temp
                #         indices.append(index_new)
                #     else:
                #         indices.append(index)
                #     if j == 0: 
                #         xx= test_points[2848+index,0]
                #     test_points[2848+index,:] = temp1_p[j,:]                 
            #temp2 = vedo.Points(test_points[2848:2848+168,:], c = 'blue', r=20)
            mesh = vedo.Mesh([test_points, self.tri_connectivity])#points, np.array(B['triangles'].tolist()).reshape((length,3))
            mesh.backColor().lineColor('green').lineWidth(3)
            vd_test = vedo.Plotter(axes=0)#, mesh
            vd_test.show(mesh, vd_points1, temp0, temp1, 'Test',viewup="z", interactive=True)   #, vd_points0        
    
    def identify_intersection(self, geo, u0_v0_vec, pointset0, pointset1, num_points_u1, num_points_v1, edge_location, relationship):
        if relationship == '-':
            print("Sorry! This has not been implemented yet!") 
        if relationship == '+':
            print("Sorry! This has not been implemented yet!")                       
        if relationship == 'T':
            print('T')   
            projecting_points, pointset1_indices = self.discritize_ctrl_pointset_edge(pointset1, num_points_u = num_points_u1, num_points_v = num_points_v1, edge_location = edge_location)
            #points = self.discritize_ctrl_pointset(pointset0, uv_vec = u0_v0_vec, num_points_u = num_points_u0, num_points_v = num_points_v0)
            num_points = len(projecting_points)
            constrained_u0_v0_vec, constrained_points = self.project_points(projecting_points, num_points, pointset0)
            #print('constrained_u0_v0_vec',num_points, constrained_u0_v0_vec.shape, constrained_u0_v0_vec)#
            indices = []
            constrained_edge = np.empty((0,2),dtype = int)
            #print('self.constrained_indices', self.constrained_indices)

            for i in range(num_points):
                #dist = np.linalg.norm(points - constrained_points[i,:], axis = 1)
                dist = np.linalg.norm(u0_v0_vec - constrained_u0_v0_vec[i,:], axis = 1)
                index = np.argmin(dist)
                if not (i == num_points - 1 or i==0):
                    while index in self.constrained_node_indices and np.min(dist) > 1e-6: #
                        print('xxxindex',index, np.min(dist))
                        dist[index] = 1000
                        index = np.argmin(dist)
                        #print('TTTindex', index)
                indices.append(index)
                u0_v0_vec[index,:] = 1000
            #print('indices', indices)            
            constrained_node_indices = indices
            for i in range(num_points):
                #print(indices[i], )
                u0_v0_vec[indices[i],:] = constrained_u0_v0_vec[i,:]

            if self.pointset_ini_check:
                for i in range(num_points):
                    if i != num_points - 1:
                        constrained_edge = np.append(constrained_edge, np.array([[indices[i], indices[i+1]]]).reshape(1,2), axis = 0)
                connectivity = self.connectivity_dict[pointset1.name]  
                sequ = range(self.num_of_nodes, self.num_of_nodes + num_points_u1*num_points_v1) 
                pointset1_indices = list(pointset1_indices)  
                
                j = 0   
                for i in range(num_points_u1*num_points_v1):
                    if sequ[i] in pointset1_indices:
                        connectivity[connectivity == sequ[i]] = indices[pointset1_indices.index(sequ[i])] + self.num_of_nodes_ini
                    else: 
                        connectivity[connectivity == sequ[i]] = range(self.num_of_nodes, self.mapping.shape[0])[j]
                        j+=1
                
                self.connectivity_dict[pointset1.name] = connectivity 
                self.members_dict[pointset1.name].options['node_indices'] = constrained_node_indices +[*range(self.num_of_nodes, self.mapping.shape[0])]
                self.num_of_nodes += len(range(self.num_of_nodes, self.mapping.shape[0]))
                return u0_v0_vec, constrained_edge, constrained_node_indices
            else:
                print('Sorry')
                return u0_v0_vec

            #print('indices',indices)         
            #print('constrained_edge', constrained_edge)
        


    def discritize_ctrl_pointset_edge(self, pointset1, num_points_u = 10, num_points_v = 8, edge_location = 'u0', order_u_fitted = 4, order_v_fitted = 4):
        num_control_points_u_fitted = pointset1.shape[0]
        num_control_points_v_fitted = pointset1.shape[1]
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
        basis0 = sps.csc_matrix(
            (data, (row_indices, col_indices)), 
            shape=(num_points_u_fitted * num_points_v_fitted, num_control_points_u_fitted * num_control_points_v_fitted),
        )

        points = basis0.dot(pointset1.physical_coordinates)   
        indices_u0 = num_points_v * np.arange(num_points_u)  
        indices_v0 = np.arange(num_points_v)
        indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1
        indices_v1 = np.arange(num_points_v)+ (num_points_u*num_points_v) - num_points_v            
        if edge_location == 'u0':            
            indices = indices_u0
        elif edge_location == 'v0': 
            indices = indices_v0  
        elif edge_location == 'u1': 
            indices = indices_u1
        elif edge_location == 'v1':
            indices = indices_v1      
        else:
            print('Warning4')          
        
        u1_v1_vec = np.vstack((u_vec, v_vec)).T
        A = dict(vertices=u1_v1_vec)
        B = tr.triangulate(A,'pc')
        points = points[indices, :]   
        if pointset1 not in self.connected_pointset_list:
            self.connected_pointset_list.append(pointset1)

            self.members_dict[pointset1.name] = Member(    
                name = pointset1.name,   
                tri_connectivity = B['triangles'],
                constrained_node_indices = [0],
                )   
            #points, np.array(B['triangles'].tolist()).reshape((length,3))
 
            B['triangles'] += self.num_of_nodes
            self.connectivity_dict[pointset1.name] = B['triangles']
            #print('mapping', mapping.shape)
            mapping = basis0.dot(pointset1.absolute_map) 
            mapping = sps.csc_matrix(np.delete(mapping.toarray(), indices, 0))#
            # print('mapping', mapping.shape)                
            # print('self.num_of_nodes', self.num_of_nodes)
            # print('mapping', type(mapping))
            # print('indices', len(indices))
            self.mapping = sps.vstack([self.mapping, mapping])
            self.mapping = sps.csc_matrix(self.mapping)  
            indices += self.num_of_nodes 
            #self.members_dict[pointset1.name].options['node_indices'] = list(indices)    
            # print('self.num_of_nodes', self.num_of_nodes)                        
        else:
            print('sorry')
            # vd_points0 = vedo.Points(points, r=10, c='black',alpha=1.0) 
            # points = points[indices, :]              
            # vd_points1 = vedo.Points(points, r=15, c='red',alpha=0.5)
            # vd_test = vedo.Plotter(axes=1)
            # vd_test.show(vd_points0, vd_points1, 'Test', viewup="z", interactive=True)  #mesh,  
            # exit()


        return points, indices

    def discritize_ctrl_pointset(self, ctrl_pointset, uv_vec, num_points_u = 10, num_points_v = 8, order_u_fitted = 4, order_v_fitted = 4):
        num_control_points_u_fitted = ctrl_pointset.shape[0]
        num_control_points_v_fitted = ctrl_pointset.shape[1]
        nnz = num_points_u * num_points_v * order_u_fitted * order_v_fitted
        data = np.zeros(nnz)
        row_indices = np.zeros(nnz, np.int32)
        col_indices = np.zeros(nnz, np.int32)
        knot_vector_u = np.zeros(num_control_points_u_fitted+order_u_fitted)
        knot_vector_v = np.zeros(num_control_points_v_fitted+order_v_fitted)
        surf_u_vec = uv_vec[:,0]
        surf_v_vec = uv_vec[:,1]   
        get_open_uniform(order_u_fitted, num_control_points_u_fitted, knot_vector_u)
        get_open_uniform(order_v_fitted, num_control_points_v_fitted, knot_vector_v)
        get_basis_surface_matrix(
            order_u_fitted, num_control_points_u_fitted, 0, surf_u_vec, knot_vector_u,
            order_v_fitted, num_control_points_v_fitted, 0, surf_v_vec, knot_vector_v,
            num_points_u * num_points_v, data, row_indices, col_indices,
                )
        basis0 = sps.csc_matrix(
            (data, (row_indices, col_indices)), 
            shape=(num_points_u * num_points_v, num_control_points_u_fitted * num_control_points_v_fitted),
        )
        points = basis0.dot(ctrl_pointset.physical_coordinates)
        return points, basis0  

    def project_points(self, points, num_points, pointset0):
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
        uv_vec = np.append(u_vec.reshape((num_points,1)), v_vec.reshape((num_points,1)), axis=1)
        return uv_vec, pts
 

    def create_triangulation(self):
        pass

    def pymeshopt(self):
        fixedvert = np.array([])
        print(len(self.projected))
        for i in range(len(self.projected)):
            fixedvert = np.concatenate((fixedvert,self.projected[i]))
        fixedvert = self.constrained_node_indices.astype(np.int32)
        itr=[0,1,2,1,0,1,0,1,3]
        
        # m = meshopt(self.pointCoords.astype(np.float32),self.conne.astype(np.int32),self.quadlist.astype(np.int32),w1=0.5, w2=0.5,w3=0.45,w4=0.9,itr =itr,plot = False,fixedvert,self.dupvert, ref_geo)
        # m.optimization()
        # self.skinpts = m.vertexCoords
        # self.skinquads = m.quadlist
        # m.saveasvtk(name)   

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