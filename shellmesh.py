from typing import Set
import numpy as np
import scipy.sparse as sps

import vedo
import triangle as tr
import matplotlib.pyplot as plt
import meshio

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
        self.num_of_nodes = 0
        self.connected_pointset_list = []
        self.members_dict = {}
        self.tri_connectivity = np.empty((0,3), dtype=np.int32) 
        self.num_of_members = 0

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
        basis_matrix = sps.csc_matrix(
            (data, (row_indices, col_indices)),
            shape=(num_points, num_control_points_u * num_control_points_v),
        )
        relative_map = np.zeros((num_points,len(geo.total_cntrl_pts_vector)))
        linear_map = basis_matrix.dot(np.identity(num_control_points_u * num_control_points_v))
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
        self.pointset_ini_check = 0
        self.ini_starting_ind = self.num_of_nodes

        num_points_u0, num_points_v0  = intersection_list[0][2], intersection_list[0][3]       
        u0_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u0), np.ones(num_points_v0)).flatten()
        v0_vec = np.einsum('i,j->ij', np.ones(num_points_u0), np.linspace(0., 1., num_points_v0)).flatten()            
        u0_v0_vec = np.vstack((u0_vec, v0_vec)).T
          
        pointset_ini = geo.pointsets_dict[intersection_list[0][0]] 
        if pointset_ini not in self.connected_pointset_list:
            geo.assemble(pointset = pointset_ini)
            _ = geo.evaluate(pointset = pointset_ini)   
            self.pointset_ini_check = 1
            self.connected_pointset_list.append(pointset_ini)
            self.num_of_members += 1
            self.num_of_nodes += len(u0_v0_vec)    
            basis_matrix = self.discritize_ctrl_pointset(pointset_ini, uv_vec = u0_v0_vec)

            constrained_edges = np.empty((0,2),dtype = int)
            indices_u0 = num_points_v0 * np.arange(num_points_u0)  
            indices_v0 = np.arange(num_points_v0)
            indices_u1 = num_points_v0 * np.arange(num_points_u0)+num_points_v0-1
            indices_v1 = np.arange(num_points_v0)+ (num_points_u0*num_points_v0) - num_points_v0
            for i in range(num_points_u0):
                if i != num_points_u0 - 1:
                    constrained_edges = np.append(constrained_edges, np.array([[indices_u0[i], indices_u0[i+1]]]).reshape(1,2), axis = 0) 
                    constrained_edges = np.append(constrained_edges, np.array([[indices_u1[i], indices_u1[i+1]]]).reshape(1,2), axis = 0)
            for i in range(num_points_v0):
                if i != num_points_v0 - 1:
                    constrained_edges = np.append(constrained_edges, np.array([[indices_v0[i], indices_v0[i+1]]]).reshape(1,2), axis = 0) 
                    constrained_edges = np.append(constrained_edges, np.array([[indices_v1[i], indices_v1[i+1]]]).reshape(1,2), axis = 0)

            self.members_dict[pointset_ini.name] = Member(    
                name = pointset_ini.name[9:],
                id = self.num_of_members - 1,
                node_indices = [*range(self.ini_starting_ind, self.num_of_nodes)],
                constrained_node_indices = [],
                u_v_vec = np.copy(u0_v0_vec),
                constrained_edges = constrained_edges,
                mapping = sps.csc_matrix(basis_matrix.dot(pointset_ini.absolute_map)),
                constrained_boundary_node_indices =[]
                )                             
        else:
            if not np.array_equal(self.members_dict[pointset_ini.name].options['u_v_vec'], u0_v0_vec) and intersection_list[0][7] == 'T':
                print('Warning5')
 
        for intersection in intersection_list:
            pointset0 = geo.pointsets_dict[intersection[0]]
            if pointset0 != pointset_ini:
                print('Warning2')         
            pointset1 = geo.pointsets_dict[intersection[1]]
            if pointset1 not in self.connected_pointset_list:
                geo.assemble(pointset = pointset1)
                _ = geo.evaluate(pointset = pointset1) 
            self.identify_intersection(
                np.copy(self.members_dict[pointset0.name].options['u_v_vec']),
                pointset0, pointset1, 
                num_points_u0 = intersection[2], num_points_v0 = intersection[3], 
                num_points_u1 = intersection[4], num_points_v1 = intersection[5], 
                edge_location = intersection[6], relationship = intersection[7]
                )
        
        if self.pointset_ini_check:
            basis_matrix = self.discritize_ctrl_pointset(pointset_ini, uv_vec = self.members_dict[pointset_ini.name].options['u_v_vec'])            
            self.members_dict[pointset_ini.name].options['mapping'] = sps.csc_matrix(basis_matrix.dot(pointset_ini.absolute_map))
        else:   
            pass                                          
        
        A = dict(vertices=self.members_dict[pointset0.name].options['u_v_vec'], segments=self.members_dict[pointset_ini.name].options['constrained_edges'])
        B = tr.triangulate(A,'p')  
        self.members_dict[pointset_ini.name].options['tri_connectivity'] = B['triangles']
        connectivity_check = np.copy(B['triangles'])
        connectivity_check_list = list(np.copy(B['triangles']).flatten())
        if np.amax(connectivity_check) >= len(self.members_dict[pointset0.name].options['u_v_vec']):
            print('WWWWD')
            # print(len(connectivity_check))
            # print(len(connectivity_check_list))
            # print(max(connectivity_check_list))
            m = max(connectivity_check_list)
            max_flatten_indices_list = [i for i, j in enumerate(connectivity_check_list) if j == m]
            #print(max_flatten_indices_list)
            max_indices_list = []
            for i in max_flatten_indices_list:
                print(i, (i)//3, connectivity_check[(i)//3,:]) 
                max_indices_list.append((i)//3)
                #print()
            connectivity_check = np.delete(connectivity_check, max_indices_list,axis = 0)
            print(len(connectivity_check))
            if pointset1.name != 'ctrl_pts_OML_lower_wing':
                print(pointset1.name)
                connectivity_check = np.append(connectivity_check, np.array([[1125,1095,1096], [1095,1065,1066], [1097,1096,1066], [1066,1095,1096]], dtype = np.int32), axis = 0)#
            else:
                print(pointset1.name)
                #plot = 1
                connectivity_check = np.append(connectivity_check, np.array([[1125,1095,1096], [1095,1065,1066], [1097,1095,1066], [1097,1095,1096]], dtype = np.int32), axis = 0)  #
            self.members_dict[pointset_ini.name].options['tri_connectivity'] = connectivity_check
            
        connectivity_check_list = list(np.copy(connectivity_check).flatten())    
        if np.amax(connectivity_check) >= len(self.members_dict[pointset0.name].options['u_v_vec']):
            print('RRRRD')
            # print(len(connectivity_check))
            # print(len(connectivity_check_list))
            # print(max(connectivity_check_list))
            m = max(connectivity_check_list)
            max_flatten_indices_list = [i for i, j in enumerate(connectivity_check_list) if j == m]
            #print(max_flatten_indices_list)
            max_indices_list = []
            for i in max_flatten_indices_list:
                print(i, (i)//3, connectivity_check[(i)//3,:]) 
                max_indices_list.append((i)//3)
                #print()
            connectivity_check = np.delete(connectivity_check, max_indices_list,axis = 0)
            print(len(connectivity_check))
            if pointset1.name != 'ctrl_pts_OML_lower_wing':
                print(pointset1.name)
                connectivity_check = np.append(connectivity_check, np.array([[404,434,435], [374,404,405], [435,405,406], [404,405,435]], dtype = np.int32), axis = 0)#
                #plot = 1
            else:
                print(pointset1.name)
                #plot = 1
                connectivity_check = np.append(connectivity_check, np.array([[404,434,435], [374,404,405], [404,405,406], [404,435,406]], dtype = np.int32), axis = 0)  #
            self.members_dict[pointset_ini.name].options['tri_connectivity'] = connectivity_check
    
        if plot:
            tr.compare(plt, A, B)         
            plt.show(block = False) 
            # plt.show()
            # exit()
            test_points1 = self.members_dict[pointset_ini.name].options['mapping'].dot(self.total_cntrl_pts_vector)
            mesh = vedo.Mesh([test_points1, connectivity_check])
            mesh.backColor().lineColor('green').lineWidth(3)
            vd_points1 = vedo.Points(test_points1, r=15, c='red',alpha=0.5)  
            vd_points2 = vedo.Points(test_points1[[404,405,435],:].reshape(3,3), r=20, c='black')  #1379, 1340
            vd_test = vedo.Plotter(axes=1)
            vd_test.show(mesh, vd_points1,vd_points2, 'Test', viewup="z", interactive=True) 
            #exit()
   
    def identify_intersection(self, u0_v0_vec, pointset0, pointset1, num_points_u0, num_points_v0, num_points_u1, num_points_v1, edge_location, relationship):
        if relationship == '-':
            print('-', pointset0.name, pointset1.name) 
            basis_matrix0 = self.discritize_ctrl_pointset(pointset0, uv_vec = self.members_dict[pointset0.name].options['u_v_vec'])
            basis_matrix1 = self.discritize_ctrl_pointset(pointset1, uv_vec = self.members_dict[pointset1.name].options['u_v_vec'])
            mapping0 = basis_matrix0.dot(pointset0.absolute_map)#self.members_dict[pointset0.name].options['mapping']
            mapping1 = basis_matrix1.dot(pointset1.absolute_map)#self.members_dict[pointset1.name].options['mapping']
            indices_list = [0,0]
            pointset = [pointset0, pointset1]
            for i in range(2):
                constrained_edge_location = edge_location[i]
                if i == 0:
                    num_points_u = num_points_u0
                    num_points_v = num_points_v0
                else:
                    num_points_u = num_points_u1
                    num_points_v = num_points_v1       
                if constrained_edge_location == 'u0':            
                    indices_list[i] = list(num_points_v * np.arange(num_points_u))
                elif constrained_edge_location == 'v0': 
                    indices_list[i] = list(np.arange(num_points_v))
                elif constrained_edge_location == 'u1': 
                    indices_list[i] = list(num_points_v * np.arange(num_points_u)+num_points_v-1)
                elif constrained_edge_location == 'v1':
                    indices_list[i] = list(np.arange(num_points_v)+ (num_points_u*num_points_v) - num_points_v)      
                else:
                    print('Warning9') 

            points0 = mapping0.dot(self.total_cntrl_pts_vector)[indices_list[0],:]
            points1 = mapping1.dot(self.total_cntrl_pts_vector)[indices_list[1],:]
            diff = np.linalg.norm(points0-points1, axis =1)
            if len(diff[diff>1e-1]) != 0:
                print('Warning10')

            for i in range(2):
                self.members_dict[pointset[i].name].options['constrained_node_indices'] += indices_list[i]
                #self.members_dict[pointset[i].name].options['constrained_node_indices'] = list(set(self.members_dict[pointset[i].name].options['constrained_node_indices'])) 
                self.members_dict[pointset[i].name].options['constrained_boundary_node_indices'] += indices_list[i]  
                self.members_dict[pointset[i].name].options['constrained_boundary_node_indices'] = list(set(self.members_dict[pointset[i].name].options['constrained_boundary_node_indices']))                 
            mapping = sps.csc_matrix(np.delete(mapping1.toarray(), self.members_dict[pointset1.name].options['constrained_boundary_node_indices'], 0))
            self.members_dict[pointset1.name].options['mapping'] = mapping 

            pointset0_node_indices = np.copy(np.array(self.members_dict[pointset0.name].options['node_indices']))
            pointset0_indices = pointset0_node_indices[indices_list[0]]
            pointset1_node_indices = np.copy(np.array(self.members_dict[pointset1.name].options['node_indices']))
            pointset1_indices = pointset1_node_indices[indices_list[1]] 
 
            new_list = []
            old_list = []   
            for mem in self.members_dict.values():
                if mem.options['id'] != self.members_dict[pointset1.name].options['id']:
                    node_indices = np.copy(np.array(mem.options['node_indices']))
                    node_indices_old = np.copy(np.array(mem.options['node_indices']))
                    for i in pointset1_indices:  
                        node_indices[node_indices_old>i] -= 1         
                        temp = node_indices[node_indices_old == i]
                        if len(temp) ==1:
                            new_list.append(temp[0])
                            old_list.append(i)
                        elif len(temp) ==0:
                            pass
                        else:
                            print('Warning7')
                            exit()                            
                    mem.options['node_indices'] = list(node_indices)  
            
            starting_pointset1_index = 0
            for mem in self.members_dict.values():               
                if mem.options['id'] < self.members_dict[pointset1.name].options['id']:
                    starting_pointset1_index += len(mem.options['mapping'].toarray())                 
            pointset1_length = self.members_dict[pointset1.name].options['mapping'].shape[0]
            
            j = 0
            for i in range(num_points_u1*num_points_v1):
                if i in indices_list[1]:                        
                    pointset1_node_indices[i] = pointset0_indices[indices_list[1].index(i)]  
                if i in self.members_dict[pointset1.name].options['constrained_node_indices']:
                    if i in self.members_dict[pointset1.name].options['constrained_boundary_node_indices']:
                        pass
                    else:
                        pointset1_node_indices[i] = range(starting_pointset1_index, starting_pointset1_index + pointset1_length)[j]
                        j+=1  
                    if i in indices_list[1]:                        
                        pointset1_node_indices[i] = pointset0_indices[indices_list[1].index(i)]  
                else:
                    pointset1_node_indices[i] = range(starting_pointset1_index, starting_pointset1_index + pointset1_length)[j]
                    j+=1  
            self.members_dict[pointset1.name].options['node_indices'] = list(pointset1_node_indices)

            for mem in self.members_dict.values():
                if mem.options['id'] > self.members_dict[pointset1.name].options['id']:
                    node_indices = np.copy(np.array(mem.options['node_indices']))
                    for i in new_list:
                        node_indices[node_indices==i] = pointset0_indices[pointset1_indices.index(old_list[new_list.index(i)])]
                    mem.options['node_indices'] = list(node_indices)  

            for i in range(2):
                if i == 0:
                    num_points_u = num_points_u0
                    num_points_v = num_points_v0
                else:
                    num_points_u = num_points_u1
                    num_points_v = num_points_v1  
                if edge_location[0] == 'u1':
                    self.members_dict[pointset[i].name].options['constrained_node_indices'] += list(np.arange(num_points_v))
                    self.members_dict[pointset[i].name].options['constrained_node_indices'] += list(np.arange(num_points_v)+ (num_points_u*num_points_v) - num_points_v)
                    #self.members_dict[pointset[i].name].options['constrained_node_indices'] = list(set(self.members_dict[pointset[i].name].options['constrained_node_indices']))                     
                    self.members_dict[pointset[i].name].options['constrained_boundary_node_indices'] += list(np.arange(num_points_v)) 
                    self.members_dict[pointset[i].name].options['constrained_boundary_node_indices'] += list(np.arange(num_points_v)+ (num_points_u*num_points_v) - num_points_v)
                    self.members_dict[pointset[i].name].options['constrained_boundary_node_indices'] = list(set(self.members_dict[pointset[i].name].options['constrained_boundary_node_indices']))                                  
        elif relationship == '+':
            print("Sorry! This has not been implemented yet!")                       
        elif relationship == 'T':
            print('T', pointset1.name) 
            u0_v0_vec_ini = np.copy(u0_v0_vec)
            pointset1_points, indices_pointset1 = self.construct_constrained_triangulation(pointset1, num_points_u = num_points_u1, num_points_v = num_points_v1, constrained_edge_location = edge_location)        
            indices_pointset1 = list(indices_pointset1)  
            num_points = len(pointset1_points)            
            constrained_u0_v0_vec = self.project_points(pointset1_points, num_points, pointset0)
            starting_pointset1_index = 0
            starting_pointset0_index = 0
            for mem in self.members_dict.values():
                if mem.options['id'] < self.members_dict[pointset0.name].options['id']:
                    starting_pointset0_index += len(mem.options['mapping'].toarray())                
                if mem.options['id'] < self.members_dict[pointset1.name].options['id']:
                    starting_pointset1_index += len(mem.options['mapping'].toarray())            

            pointset0_indices = []
            indices_pointset0 = []                        
            pointset1_indices = []
            for i in range(num_points):
                dist = np.linalg.norm(u0_v0_vec - constrained_u0_v0_vec[i,:], axis = 1)
                index = np.argmin(dist)                
                pointset0_index = self.members_dict[pointset0.name].options['node_indices'][index]
                
                if not (i == num_points - 1 or i==0):
                    while index in self.members_dict[pointset0.name].options['constrained_node_indices']:
                        #print('TESTindex',i, pointset0_index, np.min(dist))
                        dist[index] = 1000
                        index = np.argmin(dist)
                        pointset0_index = self.members_dict[pointset0.name].options['node_indices'][index]
                else:
                    '''Special case for test_eVTOL_shellmesh_0.py
                    if i==num_points - 1 and 'rib5' in pointset1.name and 'rear' in pointset0.name:'''
                    if not self.pointset_ini_check:
                        while index in self.members_dict[pointset0.name].options['constrained_node_indices'] and index not in self.members_dict[pointset0.name].options['constrained_boundary_node_indices']:
                            print('XTTTESTindex',i, pointset0_index, np.min(dist), pointset0.name, pointset1.name)
                            dist[index] = 1000
                            index = np.argmin(dist)
                            pointset0_index = self.members_dict[pointset0.name].options['node_indices'][index]

                pointset0_indices.append(pointset0_index)
                indices_pointset0.append(index)
                pointset1_indices.append(self.members_dict[pointset1.name].options['node_indices'][indices_pointset1[i]])
                u0_v0_vec[index,:] = 1000           
            
            constrained_edges = np.empty((0,2),dtype = int)
            for i in range(num_points):
                if i != num_points - 1:
                    constrained_edges = np.append(constrained_edges, np.array([[indices_pointset0[i], indices_pointset0[i+1]]]).reshape(1,2), axis = 0)
            self.members_dict[pointset0.name].options['constrained_edges'] = np.append(self.members_dict[pointset0.name].options['constrained_edges'], constrained_edges, axis = 0)

            pointset1_length = self.members_dict[pointset1.name].options['mapping'].shape[0]
            if pointset1 not in self.connected_pointset_list:
                self.connected_pointset_list.append(pointset1)
                
                pointset1_check = 0
                self.num_of_nodes += pointset1_length 
                for i in range(num_points):
                    u0_v0_vec[indices_pointset0[i],:] = constrained_u0_v0_vec[i,:]
                self.members_dict[pointset0.name].options['u_v_vec'] = u0_v0_vec                   
            else:
                pointset1_check = 1
                indices_pointset1_reduced = []
                pointset1_indices_reduced = []
                indices_pointset0_reduced = []
                pointset0_indices_reduced = []
                for i in indices_pointset0:
                    if i not in self.members_dict[pointset0.name].options['constrained_node_indices']:
                        indices_pointset1_reduced.append(indices_pointset1[indices_pointset0.index(i)]) 
                        pointset1_indices_reduced.append(pointset1_indices[indices_pointset0.index(i)]) 
                        indices_pointset0_reduced.append(i)  
                        pointset0_indices_reduced.append(pointset0_indices[indices_pointset0.index(i)])            
                self.num_of_nodes -= len(indices_pointset1_reduced)
                for i in range(num_points):
                    if indices_pointset0[i] in indices_pointset0_reduced:
                        u0_v0_vec[indices_pointset0[i],:] = constrained_u0_v0_vec[i,:]
                    else:
                        u0_v0_vec[indices_pointset0[i],:] = u0_v0_vec_ini[indices_pointset0[i],:] 
                self.members_dict[pointset0.name].options['u_v_vec'] = u0_v0_vec

                new_list = []
                old_list = []
                for mem in self.members_dict.values():
                    if mem.options['id'] != self.members_dict[pointset1.name].options['id']:
                        node_indices = np.copy(np.array(mem.options['node_indices']))
                        node_indices_old = np.copy(np.array(mem.options['node_indices']))
                        for i in pointset1_indices_reduced:  
                            node_indices[node_indices_old>i] -= 1         

                            temp = node_indices[node_indices_old == i]
                            if len(temp) ==1:
                                new_list.append(temp[0])
                                old_list.append(i)
                            elif len(temp) ==0:
                                pass
                            else:
                                print('Warning7')
                                exit()                            
                        mem.options['node_indices'] = list(node_indices)      
            
            self.members_dict[pointset0.name].options['constrained_node_indices'] += indices_pointset0
            #self.members_dict[pointset0.name].options['constrained_node_indices'] = list(set(self.members_dict[pointset0.name].options['constrained_node_indices']))
            
            if not self.pointset_ini_check: 
                basis_matrix = self.discritize_ctrl_pointset(pointset0, uv_vec = self.members_dict[pointset0.name].options['u_v_vec']) 
                constrained_indices_pointset0 = self.members_dict[pointset0.name].options['constrained_node_indices']     
                mapping = self.members_dict[pointset0.name].options['mapping']
                node_indices = np.array(self.members_dict[pointset0.name].options['node_indices']) 
                j = 0
                for i in range(num_points_u0*num_points_v0): 
                    if i in constrained_indices_pointset0 and node_indices[i]:
                        if node_indices[i] > starting_pointset0_index and node_indices[i] < starting_pointset0_index + len(mapping.toarray()):
                            if i in indices_pointset0_reduced: 
                                mapping[j,:] = basis_matrix.dot(pointset0.absolute_map)[i,:] 
                            j += 1                
                    else:
                        j += 1 
                self.members_dict[pointset0.name].options['mapping'] = mapping
                
            pointset0_indices = []
            for index in indices_pointset0:                
                pointset0_index = self.members_dict[pointset0.name].options['node_indices'][index]
                pointset0_indices.append(pointset0_index)

            pointset1_node_indices = np.copy(np.array(self.members_dict[pointset1.name].options['node_indices']))                  
            constrained_indices_pointset1 = self.members_dict[pointset1.name].options['constrained_node_indices']
            constrained_boundary_indices_pointset1 = self.members_dict[pointset1.name].options['constrained_boundary_node_indices']
            j = 0
            for i in range(num_points_u1*num_points_v1):
                if i in constrained_indices_pointset1:
                    if i in constrained_boundary_indices_pointset1:
                        pass
                    else:
                        pointset1_node_indices[i] = range(starting_pointset1_index, starting_pointset1_index + pointset1_length)[j]
                        j+=1  
                    if i in indices_pointset1:                        
                        pointset1_node_indices[i] = pointset0_indices[indices_pointset1.index(i)]  
                else: 
                    pointset1_node_indices[i] = range(starting_pointset1_index, starting_pointset1_index + pointset1_length)[j]
                    j+=1  
            self.members_dict[pointset1.name].options['node_indices'] = list(pointset1_node_indices)

            if pointset1_check and len(new_list)!=0:
                for mem in self.members_dict.values():
                    if mem.options['id'] > self.members_dict[pointset1.name].options['id']:
                        node_indices = np.copy(np.array(mem.options['node_indices']))
                        for i in new_list:
                            node_indices[node_indices==i] = pointset0_indices[pointset1_indices.index(old_list[new_list.index(i)])]
                        mem.options['node_indices'] = list(node_indices)               

        else:
            print('Warning4')
        
    def construct_constrained_triangulation(self, pointset1, num_points_u, num_points_v, constrained_edge_location, order_u_fitted = 4, order_v_fitted = 4, plot = False):
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
        u_v_vec = np.vstack((u_vec, v_vec)).T
        knot_vector_u = np.zeros(num_control_points_u_fitted+order_u_fitted)
        knot_vector_v = np.zeros(num_control_points_v_fitted+order_v_fitted)
        get_open_uniform(order_u_fitted, num_control_points_u_fitted, knot_vector_u)
        get_open_uniform(order_v_fitted, num_control_points_v_fitted, knot_vector_v)
        get_basis_surface_matrix(
            order_u_fitted, num_control_points_u_fitted, 0, u_vec, knot_vector_u,
            order_v_fitted, num_control_points_v_fitted, 0, v_vec, knot_vector_v,
            num_points_u_fitted * num_points_v_fitted, data, row_indices, col_indices,
                )
        if pointset1 not in self.connected_pointset_list:
            basis_matrix = sps.csc_matrix(
                (data, (row_indices, col_indices)), 
                shape=(num_points_u_fitted * num_points_v_fitted, num_control_points_u_fitted * num_control_points_v_fitted),
            )  
        else:
            basis_matrix = self.discritize_ctrl_pointset(pointset1, uv_vec = self.members_dict[pointset1.name].options['u_v_vec'])
        indices_u0 = num_points_v * np.arange(num_points_u)  
        indices_v0 = np.arange(num_points_v)
        indices_u1 = num_points_v * np.arange(num_points_u)+num_points_v-1
        indices_v1 = np.arange(num_points_v)+ (num_points_u*num_points_v) - num_points_v            
        if constrained_edge_location == 'u0':            
            indices = indices_u0
        elif constrained_edge_location == 'v0': 
            indices = indices_v0  
        elif constrained_edge_location == 'u1': 
            indices = indices_u1
        elif constrained_edge_location == 'v1':
            indices = indices_v1      
        else:
            print('Warning3')            
        
        pts = basis_matrix.dot(pointset1.physical_coordinates)  
        points = pts[indices, :]
            
        if pointset1 not in self.connected_pointset_list:
            constrained_edges = np.empty((0,2),dtype = int)  
            for i in range(num_points_u):
                if i != num_points_u - 1:
                    constrained_edges = np.append(constrained_edges, np.array([[indices_u0[i], indices_u0[i+1]]]).reshape(1,2), axis = 0) 
                    constrained_edges = np.append(constrained_edges, np.array([[indices_u1[i], indices_u1[i+1]]]).reshape(1,2), axis = 0)
            for i in range(num_points_v):
                if i != num_points_v - 1:
                    constrained_edges = np.append(constrained_edges, np.array([[indices_v0[i], indices_v0[i+1]]]).reshape(1,2), axis = 0) 
                    constrained_edges = np.append(constrained_edges, np.array([[indices_v1[i], indices_v1[i+1]]]).reshape(1,2), axis = 0)                                 
            A = dict(vertices=u_v_vec, segments=constrained_edges)
            B = tr.triangulate(A,'p') 

            self.num_of_members += 1 
            self.members_dict[pointset1.name] = Member(    
                name = pointset1.name[9:],   
                id = self.num_of_members - 1,
                node_indices = [*range(self.num_of_nodes, self.num_of_nodes+len(u_v_vec))], 
                tri_connectivity = np.copy(B['triangles']),
                constrained_edges = constrained_edges,
                u_v_vec = u_v_vec,

                )           
            self.members_dict[pointset1.name].options['constrained_node_indices'] = list(indices)
            self.members_dict[pointset1.name].options['constrained_boundary_node_indices'] = list(indices) 
                                         
        else:
            self.members_dict[pointset1.name].options['constrained_node_indices'] += list(indices)
            #self.members_dict[pointset1.name].options['constrained_node_indices'] = list(set(self.members_dict[pointset1.name].options['constrained_node_indices'])) 
            self.members_dict[pointset1.name].options['constrained_boundary_node_indices'] += list(indices)   
            self.members_dict[pointset1.name].options['constrained_boundary_node_indices'] = list(set(self.members_dict[pointset1.name].options['constrained_boundary_node_indices']))  
                
        mapping = basis_matrix.dot(pointset1.absolute_map) 
        mapping = sps.csc_matrix(np.delete(mapping.toarray(), self.members_dict[pointset1.name].options['constrained_boundary_node_indices'], 0))
        self.members_dict[pointset1.name].options['mapping'] = mapping                 

        if plot: 
            print(np.shape(B['triangles']))
            tr.compare(plt, A, B)
            plt.show(block = False)
            vd_points0 = vedo.Points(pts, r=10, c='black',alpha=1.0)              
            vd_points1 = vedo.Points(points, r=15, c='red',alpha=0.5)
            vd_test = vedo.Plotter(axes=1)
            vd_test.show(vd_points0, vd_points1, 'Test', viewup="z", interactive=True)
        
        return points, indices

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
        surfs_u_order = np.append(surfs_u_order,np.array((order_u),ndmin=2),axis = 0)
        surfs_v_order = np.append(surfs_v_order,np.array((order_v),ndmin=2),axis = 0)
        surfs_num_control_points_u = np.append(surfs_num_control_points_u,np.array(num_control_points_u, ndmin=2),axis = 0)
        surfs_num_control_points_v = np.append(surfs_num_control_points_v,np.array(num_control_points_v, ndmin=2),axis = 0)
        surfs_cp = np.append(surfs_cp,cps.reshape((1,num_control_points_u * num_control_points_v*3)),axis = 0)
        compute_surface_projection(
            surfs_u_order.reshape((len(range(1)))),  surfs_num_control_points_u.reshape((len(range(1)))),
            surfs_v_order.reshape((len(range(1)))),  surfs_num_control_points_v.reshape((len(range(1)))),
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
        uv_vec = np.append(u_vec.reshape((num_points,1)), v_vec.reshape((num_points,1)), axis=1)
        return uv_vec
 

    def discritize_ctrl_pointset(self, ctrl_pointset, uv_vec, order_u_fitted = 4, order_v_fitted = 4):
        num_control_points_u_fitted = ctrl_pointset.shape[0]
        num_control_points_v_fitted = ctrl_pointset.shape[1]
        nnz = len(uv_vec) * order_u_fitted * order_v_fitted
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
            len(uv_vec), data, row_indices, col_indices,
                )
        basis_matrix = sps.csc_matrix(
            (data, (row_indices, col_indices)), 
            shape=(len(uv_vec), num_control_points_u_fitted * num_control_points_v_fitted),
        )
        return basis_matrix  

    def construct_whole_structure_mesh(self, plot = False):
        constrained = []   
        print()
        for memb in self.members_dict.values():
            #print(memb.options['id'], len(inde))
            if memb.options['id'] == 1 or memb.options['id'] == 2:
                num_u_spar = 56
                num_v_spar = 4
                indices_v0 = np.arange(num_v_spar)
                indices_v1 = np.arange(num_v_spar)+ (num_u_spar*num_v_spar) - num_v_spar  
                indices_list = [indices_v0, indices_v1]
                for indices in indices_list:
                    memb.options['constrained_node_indices'] += list(indices)
                    memb.options['constrained_boundary_node_indices'] += list(indices)   
                    memb.options['constrained_boundary_node_indices'] = list(set(memb.options['constrained_boundary_node_indices']))                 
            conn = np.copy(memb.options['tri_connectivity']) 
            inde = np.copy(memb.options['node_indices'])
            cons = list(np.copy(memb.options['constrained_node_indices']))
            
            # # # if memb.options['id'] == 5:
            # # #     # print(inde[0])
            # # #     # print(cons)
            # # #     # print(inde[cons])
            # # #     inde[0]= 404
            # # #     memb.options['node_indices'] = list(inde)
            # # # if memb.options['id'] == 11:
            # # #     # print(inde[0])
            # # #     # print(cons)
            # # #     # print(inde[cons])                
            # # #     inde[0]= 1095
            # # #     memb.options['node_indices'] = list(inde)

            for i in range(len(conn)):
                for j in range(3):                    
                    conn[i,j] = inde[conn[i,j]] 
            for i in range(len(cons)):
                cons[i] = inde[cons[i]] 

            self.tri_connectivity = np.append(self.tri_connectivity, conn, axis = 0)      
            self.mapping = sps.csc_matrix(sps.vstack([self.mapping, memb.options['mapping']]))
            constrained += cons 
            #print(memb.options['id'], len(memb.options['mapping'].toarray()))                    
        
        constrained = list(set(constrained))
        self.total_points = self.mapping.dot(self.total_cntrl_pts_vector)
        print('constrained', len(constrained))
        print('self.tri_connectivity', len(self.tri_connectivity))
        print('total_points', len(self.total_points))

        total_test = np.copy(self.total_points)
        for memb in self.members_dict.values():
            inde = list(np.copy(memb.options['node_indices']))
            memb.options['coordinates'] = self.total_points[inde,:]
            total_test[inde,:] = memb.options['coordinates']

        if plot:
            mesh = vedo.Mesh([self.total_points, self.tri_connectivity], alpha=0.3)
            mesh.backColor().lineColor('green').lineWidth(3)
            vd_points1 = vedo.Points(self.total_points[constrained,:], r=20, c='red',alpha = 0.7)  
            temp = [1095]
            vd_points2 = vedo.Points(self.total_points[temp,:], r=25, c='black') 
            vd_test = vedo.Plotter(axes=1)
            vd_test.show(mesh, 'Test', vd_points1, vd_points2, viewup="z", interactive=False)#True
            #exit()
            # mesh_test = vedo.Mesh([total_test, self.tri_connectivity], alpha=0.3)
            # mesh_test.backColor().lineColor('green').lineWidth(3) 
            # vd_test = vedo.Plotter(axes=1)
            # vd_test.show(mesh_test, 'Test', vd_points1, vd_points2, viewup="z", interactive=True)
    
    def save_tri_vtk(self, save_name, points, trilist):        
        cells = []
        for i in range(len(trilist)):
            tuple_tri = ("triangle",trilist[i,:].reshape(1,3))
            cells.append(tuple_tri)
        meshio.write_points_cells(
            save_name+"_tri.vtk",
            points,
            cells
            )    
        print('save as vtk finished',save_name) 

    def optimizie_mesh(self, plot = False):
        for memb in self.members_dict.values():
            #if memb.options['id'] ==1:
            print(memb.options['id'], memb.options['name'])
            quadlist = np.empty((0,4),dtype=np.int32)
            vertexCoords = np.copy(memb.options['coordinates'].astype('float32'))
            trilist = np.copy(memb.options['tri_connectivity'])
            fixedvert = np.copy(np.array(memb.options['constrained_node_indices'], dtype='int32'))
            if memb.options['id'] ==3:
                itr=[0,1,2,3]                              
                m = meshopt(vertexCoords,trilist,quadlist, fixedvert=fixedvert ,itr=itr, w1=1.,w2=1.,w3=1., plot = 0)#,w4 =1.
                m.optimization()
                m.saveasvtk('CAD/'+memb.options['name']) 
                memb.options['opt_coordinates'] = m.vertexCoords
                memb.options['opt_connectivity'] = m.quadlist

            if True:#plot:#
                mesh = vedo.Mesh([vertexCoords, trilist], alpha=0.3)
                mesh.backColor().lineColor('green').lineWidth(3)
                vd_points1 = vedo.Points(vertexCoords[fixedvert,:], r=20, c='red',alpha = 0.7)  
                vd_test = vedo.Plotter(axes=1)
                vd_test.show(mesh, 'optimizie_mesh', vd_points1, viewup="z", interactive=True)               
    
    def construct_whole_structure_optmesh(self, vtk_file_name):
        print('Construct whole structure optmesh')
        pts = np.empty((0,3)) 
        quads = np.empty((0,4), dtype=np.int32) 
        for memb in self.members_dict.values():
            quads = np.concatenate((quads,np.copy(memb.options['opt_connectivity']+len(pts))),axis=0)
            pts = np.concatenate((pts, np.copy(memb.options['opt_coordinates'])),axis=0)
        import pymeshopt
        opt = pymeshopt.Pymeshopt(pts.astype(np.float32),np.array([[0,1,2]],dtype=np.int32),quads.astype(np.int32),1.,1.,1.)
        self.uni = opt.pymergeduppts()
        self.saveasvtk(vtk_file_name)

    def saveasvtk(self, vtk_file_name):
        points = self.uni.vertlist
        quadlist = self.uni.quadlist
        print('Starting to save as vtk', len(points), len(quadlist))
        triindex1 = np.array([0,1,2],dtype=np.int32)
        triindex2 = np.array([2,3,0],dtype=np.int32)  
        cells = []  
        for i in range(len(quadlist)):
            tuple1 = ("triangle",quadlist[i,triindex1].reshape(1,3))
            tuple2 = ("triangle",quadlist[i,triindex2].reshape(1,3))
            cells.append(tuple1)
            cells.append(tuple2)
        meshio.write_points_cells(
            vtk_file_name+'_whole_model'+'_'+str(len(points))+'_tri_'+str(len(cells))+'.vtk',
            points,
            cells
            )        
        cells = []
        for i in range(len(quadlist)):
            tuple = ("quad",quadlist[i,:].reshape(1,4))
            cells.append(tuple)
        meshio.write_points_cells(
            vtk_file_name+'_whole_model'+'_'+str(len(points))+'_quad_'+str(len(cells))+'.vtk',
            points,
            cells
            ) 



if __name__ == '__main__':
    from shellmesh import ShellMesh
    shell_mesh = ShellMesh('shell_mesh')
