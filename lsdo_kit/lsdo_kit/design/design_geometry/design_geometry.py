import numpy as np
from numpy.core.fromnumeric import ndim
from numpy.core.shape_base import vstack
from scipy import sparse, interpolate
import scipy.sparse as sps
import pandas as pd
#from pandas.core.indexes.base import ensure_index
from dataclasses import dataclass
from toposort import toposort
import re
from typing import List
import copy

#import pyiges

from geomdl import fitting

# from lsdo_kit.DesignGeometry.bspline_entities import BSplineCurve, BSplineSurface, BSplineVolume, BSplineEntity
from lsdo_kit.design.design_geometry.bsplines.bspline_curve import BSplineCurve
from lsdo_kit.design.design_geometry.bsplines.bspline_surface import BSplineSurface
from lsdo_kit.design.design_geometry.bsplines.bspline_volume import BSplineVolume
from lsdo_kit.design.design_geometry.utils.io.step_io import read_openvsp_stp, write_step, read_gmsh_stp
from lsdo_kit.design.design_geometry.utils.io.iges_io import read_iges, write_iges
from lsdo_kit.design.design_geometry.core.pointset_functions import define_linear_combination, _define_linear_combination
from lsdo_kit.design.design_geometry.core.pointset_functions import perform_linear_interpolation, perform_bilinear_interpolation, perform_trilinear_interpolation
from lsdo_kit.design.design_geometry.core.pointset_functions import perform_2d_transfinite_interpolation, perform_3d_transfinite_interpolation
from lsdo_kit.design.design_geometry.core.pointset_functions import extract_pointset
from lsdo_kit.design.design_geometry.core.pointset_functions import add_pointsets, subtract_pointsets, divide_pointset_by_scalar

# from lsdo_kit.splines.basis_matrix_surface_py import get_basis_surface_matrix
# from lsdo_kit.splines.surface_projection_py import compute_surface_projection
# from lsdo_kit.splines.get_open_uniform_py import get_open_uniform

from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
from lsdo_kit.cython.surface_projection_py import compute_surface_projection
from lsdo_kit.cython.get_open_uniform_py import get_open_uniform

from lsdo_kit.design.design_geometry.core.pointset import ProjectedPointSet, DerivedPointSet, EvaluatedPointSets, PointSet


import matplotlib.pyplot as plt
from vedo import Points, Plotter, colors, LegendBox, show
import vedo # for using Mesh

class DesignGeometry:

    def __init__(self, file_name=None, plot=True):
        self.input_bspline_entity_dict = {}
        self.initial_input_bspline_entity_dict = {}
        self.total_cntrl_pts_vector = np.array([])
        self.cntrl_pts_unique = np.array([])
        self.output_geo_control_points_dict = {}
        self.components = {}
        self.components_ffd_dict = {}
        self.num_points = 0
        self.read_file(file_name)

        self.current_pointset_pts = 0
        self.registered_names = [] 
        self.current_id = 0
        self.parent_pointers_dict = {}
        self.pointsets_dict = {}
        self.eval_map = None
        self.output_bspline_entity_dict = self.initial_input_bspline_entity_dict

        self.evaluated_pointsets = EvaluatedPointSets(
            pointset_id = self.current_id,
            shape = None,
            parent_pointers = [],
            absolute_map = None,
            relative_map = None,
            offset = np.array([0., 0., 0.]),            
            offset_absolute_map = None,
            physical_coordinates = None,
            output_starting_ind = 0,
            permutation_matrix = None,
            )
        self.current_id += 1
        self.mesh_list = {}
        self.bspline_mesh_list = {}
        self.geometric_outputs = None
        print('num_surf',len(self.input_bspline_entity_dict.values()))
        if plot == True: #There are only 239 colors in vedo.colors.colors dictionary. If number of bspline surfaces is larger than 238*2, reduce the index of the color agian
            vps_ini = []
            vps = []
            for i in range(len(self.input_bspline_entity_dict.values())):  
                if i < 10 or i > 57:# 
                    pass
                else: 
                    print(i)     
                    if i > 238:
                        color = list(colors.colors.values())[i-239]#239
                    else:
                        color = list(colors.colors.values())[i]
                    surf_ini = list(self.initial_input_bspline_entity_dict.values())[i] 
                    vps_ini.append(Points(surf_ini.control_points, r=8, c = color))
                    surf = list(self.input_bspline_entity_dict.values())[i]
                    vps.append(Points(surf.control_points, r=8, c = color))
            print()
            vp_init_out = Plotter()
            vp_init_out.show(vps_ini, 'Initial control points', axes=1, viewup="z", interactive = False)
            vp_init = Plotter()
            vp_init.show(vps, 'Interpolated control points', axes=1, viewup="z", interactive = False)       


    def project_points(self, points_to_be_projected, projection_targets_names=[], projection_direction=np.array([0., 0., 0.]), offset=np.array([0., 0., 0.]), plot=False):
        if len(np.shape(points_to_be_projected))== 1:
            shape = 1 
            num_points = 1
        elif len(np.shape(points_to_be_projected))== 2:
            shape = np.shape(points_to_be_projected)[:-1]
            num_points = np.shape(points_to_be_projected)[0]
        else:
            shape = np.shape(points_to_be_projected)[:-1]
            num_points = np.shape(points_to_be_projected)[0]*np.shape(points_to_be_projected)[1]

        if len(np.shape(projection_direction)) == 1:
            projection_direction = np.array((projection_direction),ndmin=2)
            projection_direction = np.repeat(projection_direction,num_points,axis = 0)
        
        if projection_targets_names == []:#if no input target is passed
            projection_targets = list(self.input_bspline_entity_dict.values()) 
            # print(projection_targets) 
            # print(type(projection_targets))
        else:
            projection_targets = []
            for target in projection_targets_names:# create matrix of projection target identifiers
                if target in self.components.keys():
                    for entity in self.components[target].embedded_entities_objects:
                        projection_targets.append(entity)
                elif target in self.input_bspline_entity_dict.keys():
                    projection_targets.append(self.input_bspline_entity_dict[target])
                else:
                    raise Exception("When specifying target names, specify a component name or a geometry bpsline name.")

            # print(projection_targets)
        # print('Components Geo: ', self.components)
        # print('Projection Targets: ', projection_targets)
        # print('Projection Targets Names: ', projection_targets_names)
        # print('Bspline Entities Dict: ', self.input_bspline_entity_dict)
        relative_map = np.zeros((num_points,len(self.total_cntrl_pts_vector)))
        temp = 1e16 * np.ones((num_points,))
        surf_final = [None] * num_points
        cps_final = [None] * num_points
        linear_map_final = [None] * num_points

        surfs_u_order = np.empty((0,1),dtype = int)
        surfs_v_order = np.empty((0,1),dtype = int)
        surfs_num_control_points_u = np.empty((0,1),dtype = int)
        surfs_num_control_points_v = np.empty((0,1),dtype = int)

        surfs_cp = np.empty((0,projection_targets[0].shape[0]*projection_targets[0].shape[1]*3),dtype = int)
        for surf in projection_targets: 
            surfs_u_order = np.append(surfs_u_order,np.array((surf.order_u),ndmin=2),axis = 0)
            surfs_v_order = np.append(surfs_v_order,np.array((surf.order_v),ndmin=2),axis = 0)
            surfs_num_control_points_u = np.append(surfs_num_control_points_u,np.array(surf.shape[0], ndmin=2),axis = 0)
            surfs_num_control_points_v = np.append(surfs_num_control_points_v,np.array(surf.shape[1], ndmin=2),axis = 0)
            surfs_cp = np.append(surfs_cp,surf.control_points.reshape((1,projection_targets[0].shape[0]*projection_targets[0].shape[1]*3)),axis = 0)

        points = points_to_be_projected
        axis = projection_direction
        for surf in projection_targets:  
            num_control_points_u = surf.shape[0]
            num_control_points_v = surf.shape[1]
            u_order = surf.order_u
            v_order = surf.order_v
            knot_vector_u = surf.knots_u
            knot_vector_v = surf.knots_v
            cps = surf.control_points
            max_iter = 500
            u_vec = np.ones(num_points)
            v_vec = np.ones(num_points)
            surfs_index = np.zeros(num_points,dtype=int)

        compute_surface_projection(
            surfs_u_order.reshape((len(projection_targets))), surfs_num_control_points_u.reshape((len(projection_targets))),
            surfs_v_order.reshape((len(projection_targets))), surfs_num_control_points_v.reshape((len(projection_targets))),
            num_points, max_iter,
            points.reshape(num_points * 3), 
            cps.reshape((num_control_points_u * num_control_points_v * 3)),
            knot_vector_u, knot_vector_v,
            u_vec, v_vec, 50,
            axis.reshape(num_points * 3),
            surfs_index,
            surfs_cp,
        )
        #print('u_vec',u_vec)
        #print('v_vec',v_vec)
        #print('surfs_index', surfs_index)

        for ns,surf in zip(range(len(projection_targets)),projection_targets):
            num_control_points_u = surf.shape[0]
            num_control_points_v = surf.shape[1]
            u_order = surf.order_u
            v_order = surf.order_v
            knot_vector_u = surf.knots_u
            knot_vector_v = surf.knots_v
            cps = surf.control_points
            nnz = num_points * u_order * v_order
            data = np.zeros(nnz)
            row_indices = np.zeros(nnz, np.int32)
            col_indices = np.zeros(nnz, np.int32)
            get_basis_surface_matrix(
                u_order, num_control_points_u, 0, u_vec, knot_vector_u, 
                v_order, num_control_points_v, 0, v_vec, knot_vector_v,
                num_points, data, row_indices, col_indices,
            )
            basis0 = sps.csc_matrix(
                (data, (row_indices, col_indices)), 
                shape=(num_points, num_control_points_u * num_control_points_v),
            )           
            pts = basis0.dot(cps)
            pts = np.reshape(pts, (num_points, 3))
            points = np.reshape(points, (num_points, 3))
            linear_map = basis0.dot(np.identity(surf.shape[0] * surf.shape[1]))
            for i in range(np.shape(projection_direction)[0]):
                if surfs_index[i] == ns:
                    surf_final[i] = surf
                    linear_map_final[i] = linear_map[i,:]
                    cps_final[i] = cps
                    
        for i in range(np.shape(projection_direction)[0]):
            j = 0
            for surf in self.input_bspline_entity_dict.values():
                if surf == surf_final[i]:
                    relative_map[i, j:j+surf.shape[0]*surf.shape[1]] = linear_map_final[i]
                j = j + surf.shape[0]*surf.shape[1] 
        #  - Want to directly build sparse basis matrix. Must starting index of surface projection target
        #     - For now, we can loop over the shapes of the surfaces to find the starting index of the proj target
        #       - In the future, we could probably use Anugrah's Vector Class (Array Manager)
        relative_map = sps.csc_matrix(relative_map)
        pointset = ProjectedPointSet(
            pointset_id = self.current_id,
            shape = np.append(shape,3),
            output_starting_ind = self.current_pointset_pts,
            parent_pointers = [],
            absolute_map = None,
            relative_map = relative_map,
            offset = offset,            
            offset_absolute_map = None,
            physical_coordinates = None,
            permutation_matrix = None,
            )
        self.pointsets_dict[self.current_id] = pointset
        self.current_id += 1
        self.current_pointset_pts += np.cumprod(shape)[-1]
        self.register_output(pointset)
        offset = np.array(offset,ndmin=2)
        offset = np.repeat(offset,np.shape(projection_direction)[0],axis = 0)
        point_test = relative_map.dot(self.total_cntrl_pts_vector) + offset#test
        if plot == True: 
            cps1 = [] 
            cps2 = [] 
            point_test = relative_map.dot(self.total_cntrl_pts_vector) + offset
            for surf, cps in zip(surf_final, cps_final):
                cps1.append(Points(surf.control_points, r=8, c = 'b', alpha = 0.5).legend('Contol points of projected surface'))
                cps2.append(Points(cps.reshape((num_control_points_u * num_control_points_v, 3)), r=8, c='b', alpha = 0.5).legend('Interpolated contol points of projected surface'))
            projecting_point = Points(points.reshape((num_points, 3)), r=15, c='g').legend('Projecting curve')
            projected_point = Points(point_test.reshape((num_points, 3)), r=15, c='r').legend('Projected curve')
            vp_project_curve = Plotter(N=3, axes=1) #TODO legend
            vp_project_curve.show(cps1, 'Control points of surface to be projected', at=0, viewup="z")#, lb1
            vp_project_curve.show(cps2, projecting_point, 'Surface + projecting curve', at=1, viewup="z")#, lb2
            vp_project_curve.show(cps2, projected_point, 'Surface + projected curve', at=2, viewup="z",interactive=False)#, lb3
        return pointset, point_test


    def perform_linear_combination(self, parent_pointset_list, coefficients_matrix, shape, offset=np.array([0., 0., 0.])):
        """ 
        Perform any arbitrary linear combination between PointSets and creates a new PointSet

        Parameters
        ----------
        parent_pointset_list : list
            A list of the parent pointsets in the order they are used in the relative map
        relative_map : csc_matrix     
            A sprase matrix defining the linear combination s.t. the result is a flattened shape (n,3)
        shape : np.ndarray            
            A numpy array defining the desired unflattened shape of the output
        offset : np.ndarray           
            A numpy array of shape (n,3) defining how much to shift each point

        Returns
        -------
        pointset : PointSet           
            A PointSet object that is the result of the linear combination
        """
        return _define_linear_combination(self, parent_pointset_list=parent_pointset_list, relative_map=coefficients_matrix, shape=shape, offset=offset)
        

    def _define_linear_combination(self, parent_pointset_list, relative_map, shape, offset=np.array([0., 0., 0.])):
        return _define_linear_combination(self, parent_pointset_list, relative_map, shape, offset)


    def perform_linear_interpolation(self, pointset_start, pointset_end, shape, output_parameters=np.array([0]), offset=np.array([0., 0., 0.])):       
        return perform_linear_interpolation(self, pointset_start, pointset_end,shape,output_parameters)

    def perform_bilinear_interpolation(self, point_00, point_10, point_01, point_11, shape, output_parameters=np.array([0., 0., 0.]), offset=np.array([0., 0., 0.])):
        return perform_bilinear_interpolation(self, point_00, point_10, point_01, point_11, shape, output_parameters, offset)

    def perform_trilinear_interpolation(self, point_000, point_100, point_010, point_110, point_001, point_101, point_011, point_111, 
        shape, output_parameters=np.array([0., 0., 0.]), offset=np.array([0., 0., 0.])):        
        
        return perform_trilinear_interpolation(self, point_000, point_100, point_010, point_110, point_001, point_101, point_011, point_111, 
        shape, output_parameters, offset)

    def perform_2d_transfinite_interpolation(self, u_curve0, u_curve1, v_curve0, v_curve1, output_parameters=np.array([0., 0., 0.]), offset=np.array([0., 0., 0.])):
        return perform_2d_transfinite_interpolation(self, u_curve0, u_curve1, v_curve0, v_curve1, output_parameters, offset)

    def perform_3d_transfinite_interpolation(self, output_parameters=None, offset=np.array([0., 0., 0.])):
        return perform_3d_transfinite_interpolation(self, output_parameters, offset)


    def extract_pointset(self, parent_pointset, point_indices, shape, offset=np.array([0., 0., 0.])):
        return extract_pointset(self, parent_pointset, point_indices, shape, offset)

    def add_pointsets(self, pointset1, pointset2, offset=np.array([0., 0., 0.])):
        return add_pointsets(self, pointset1, pointset2, offset=offset)

    def subtract_pointsets(self, pointset1, pointset2, offset=np.array([0., 0., 0.])):
        return subtract_pointsets(self, pointset1, pointset2, offset=offset)

    def divide_pointset_by_scalar(self, pointset1, scalar, offset=np.array([0., 0., 0.])):
        return divide_pointset_by_scalar(self, pointset1, scalar, offset=offset)

    def register_output(self, pointset, name=None, mesh_names=[], bspline_mesh_names=[]):
        '''
        Create/add on to the single instance of evaluated_pointsets
        '''
        if name is None:
            name = f'pointset_{pointset.pointset_id}'

        self.evaluated_pointsets.parent_pointers.append(pointset.pointset_id)
        pointset.name = name

        # for mesh_name in mesh_names:
        #     self.mesh_list[mesh_name].register_output(pointset, name)
        # for bspline_mesh_name in bspline_mesh_names:
        #     self.bspline_mesh_list[bspline_mesh_name].register_output(pointset, name)

    def assemble_all_absolute_maps(self):
        '''
            construct appropriate DAG representation
            call a topological sorting algorithm to get a topological ordering
            call the pointsets in this order
            ask each to evaluate its absolute maps
        '''
        ordered_pointers = list(toposort(self.parent_pointers_dict))
        for set in ordered_pointers:
            for id in set:
                self.evaluate_my_absolute_map(self.pointsets_dict[id])

    # method in PointSet
    def assemble_my_absolute_map(self):
        '''
        calculates absolute map using source PointSetss' abssolute maps
        '''

    def assemble_absolute_map(self, point_set):
        #print('DesignGeometry Pointset Name:', point_set.name)
        #print('DesignGeometry Pointset: ', point_set.absolute_map)
        
        if point_set.absolute_map != None and point_set.pointset_id != 0:
            return
        else:
            input_map = None
            offset_input_map = None
            for parent_id in point_set.parent_pointers:
                # calc absolute map of parent pointsets if they aren't already calculated
                self.assemble_absolute_map(self.pointsets_dict[parent_id])
                # self.assemble_absolute_map(self.pointsets_dict[f'{parent_id}'])

                # build input map
                if input_map == None:
                    input_map = self.pointsets_dict[parent_id].absolute_map
                    offset_input_map = self.pointsets_dict[parent_id].offset_absolute_map
                    # input_map = self.pointsets_dict[f'{parent_id}'].absolute_map
                else:
                    input_map = sps.vstack((input_map, self.pointsets_dict[parent_id].absolute_map))
                    offset_input_map = np.vstack((offset_input_map, self.pointsets_dict[parent_id].offset_absolute_map))
                    # input_map = np.vstack(input_map, self.pointsets_dict[f'{parent_id}'].absolute_map)
            if input_map == None:
                # print('type pointset: ', type(point_set))
                # print('type pointset map: ', type(point_set.relative_map))
                
                input_map = sps.eye(point_set.relative_map.shape[-1])    # assuming shape is (num_outputs, num_inputs) where we want columns
                offset_input_map = np.zeros((point_set.relative_map.shape[-1], 3))
            if point_set.pointset_id == 0:
                point_set.relative_map = sps.eye(input_map.shape[0])
            if point_set.offset.shape[0] == 1:
                relative_offset = np.tile(point_set.offset, [point_set.relative_map.shape[0], 1])
            else:
                relative_offset = point_set.offset
            relative_offset = point_set.offset.reshape(np.prod(point_set.offset.shape)//3, 3)

            # print('Input map shape: ', input_map.shape)
            # print('Relative map shape: ', point_set.relative_map.shape)
            # calc absolute map of this point_set
            # print('P'point_set.relative_map.shape)
            # print(input_map.shape)
            point_set.absolute_map = point_set.relative_map.dot(input_map)
            point_set.offset_absolute_map = np.dot(point_set.relative_map.todense(), offset_input_map) + relative_offset
            return


    def assemble(self, pointset = None):#first evalution 
        '''
        Concatenates vertically all matrices. 
        '''
        if pointset is None:
            self.assemble_absolute_map(self.evaluated_pointsets)
            self.eval_map = self.evaluated_pointsets.absolute_map
            self.offset_eval_map = self.evaluated_pointsets.offset_absolute_map
        else:
            self.assemble_absolute_map(pointset)
        return

    def evaluate(self, pointset=None):
        ''' 
        Output allows user to see the physical coordinates of pointsets.
        i.e allows user to call the physical_coordinates attribute of current pointset instances.

        (geometry.assemble() must be run first)
        '''
        if pointset == None:
            points = self.eval_map.dot(self.total_cntrl_pts_vector)# + self.offset_eval_map 
            points += self.offset_eval_map
            for pointset in self.pointsets_dict.values():
                #print(pointset.pointset_id,pointset.name)
                start_index = pointset.output_starting_ind
                end_index = start_index + np.cumprod(pointset.shape)[-2]
                pointset.physical_coordinates = points[start_index:end_index]
                # TODO: Add interface between the output vector and the output dictionary to mesh_evaluation model 

        else:
            points = pointset.absolute_map.dot(self.total_cntrl_pts_vector)
            points += pointset.offset_absolute_map
            pointset.physical_coordinates = points
        #exit()
        return points

    def fit_bspline_entities(self, pointset_list, output_vec=None):
        '''Least square fit b-spline surface'''
        # for pointer_id in self.evaluated_pointsets.parent_pointers:
        for pointset in pointset_list:
            #print(pointset.pointset_id,pointset.name)
            component_shape = pointset.shape
            #print(pointset.pointset_id,pointset.name)#9 primary_spar
            if output_vec is None:
                entity_points_be_fitted = pointset.physical_coordinates
            else:
                start_index = pointset.output_starting_ind
                end_index = start_index + np.cumprod(component_shape)[-2]
                entity_points_be_fitted = output_vec[start_index:end_index]
            #print(component_shape)
            #print(component_shape[:-1])
            #print(entity_starting_point)
            if len(component_shape[:-1]) == 0:  # is point  
                print('fitting points has not been implemented yet')
                pass        #is point
            elif len(component_shape[:-1]) == 1:  # is curve
                curve = fitting.approximate_curve(entity_points_be_fitted, 3)   #TODO hardcoded cubic bspline
                bspline_entity_curve = BSplineCurve(
                    name=pointset.name,
                    order_u=curve.order,
                    shape=np.array[len(curve.ctrlpts), 3],
                    control_points=curve.ctrlpts,
                    knots_u=curve.knotvector)
                self.output_bspline_entity_dict[pointset.name] = bspline_entity_curve
            elif len(component_shape[:-1]) == 2:  # is surface
                order_u_fitted = 4
                order_v_fitted = 4
                num_control_points_u_fitted = component_shape[0] - 1
                num_control_points_v_fitted = component_shape[1] - 1
                num_points_u_fitted = component_shape[0]
                num_points_v_fitted = component_shape[1]

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

                a = np.matmul(basis0.toarray().T, basis0.toarray())
                if np.linalg.det(a) == 0:
                    print('lstsq')
                    cps_fitted,_,_,_ = np.linalg.lstsq(a, np.matmul(basis0.toarray().T, entity_points_be_fitted), rcond=None)
                else: 
                    #print('linalg')
                    cps_fitted = np.linalg.solve(a, np.matmul(basis0.toarray().T, entity_points_be_fitted))            
            
                bspline_entity_surface = BSplineSurface(
                    name=pointset.name,
                    order_u=order_u_fitted,
                    order_v=order_v_fitted,
                    shape=np.array([num_control_points_u_fitted, num_control_points_v_fitted, 3]),
                    control_points=np.array(cps_fitted).reshape((num_control_points_u_fitted*num_control_points_v_fitted,3)),
                    knots_u=np.array(knot_vector_u),
                    knots_v=np.array(knot_vector_v))
                self.output_bspline_entity_dict[pointset.name] = bspline_entity_surface

            elif len(component_shape[:-1]) == 3:  # is volume
                print('fitting BSplineVolume has not been implemented yet')
                pass
        pass 


    def remove_multiplicity(self, bspline_entity):
        # TODO allow it to be curves or volumes too
        component_shape = np.array(bspline_entity.shape)
        if len(component_shape[:-1]) == 0:  # is point
            print('fitting points has not been implemented yet')
            pass        #is point
        elif len(component_shape[:-1]) == 1:  # is curve
            print('fitting curves has not been implemented yet')
            pass 
        elif len(component_shape[:-1]) == 2: 
            order_u = bspline_entity.order_u
            order_v = bspline_entity.order_v
            num_control_points_u = bspline_entity.shape[0]
            num_control_points_v = bspline_entity.shape[1]
            num_points_u = 20   # TODO might want to pass these in as input
            num_points_v = 20

            nnz = num_points_u * num_points_v * order_u * order_v
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
                num_points_u * num_points_v, data, row_indices, col_indices,
            )

            basis0 = sps.csc_matrix(
                (data, (row_indices, col_indices)),
                shape=(num_points_u * num_points_v, num_control_points_u * num_control_points_v),
            )

            pts = basis0.dot(bspline_entity.control_points)
            #print(pts.shape)
            order_u_fitted = 4
            order_v_fitted = 4

            num_control_points_u_fitted = 6
            num_control_points_v_fitted = 15
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

            a = np.matmul(basis0.toarray().T, basis0.toarray())
            if np.linalg.det(a) == 0:
                cps_fitted,_,_,_ = np.linalg.lstsq(a, np.matmul(basis0.toarray().T, pts), rcond=None)
            else: 
                cps_fitted = np.linalg.solve(a, np.matmul(basis0.toarray().T, pts))            
        
            bspline_entity_surface = BSplineSurface(
                name=bspline_entity.name,
                order_u=order_u_fitted,
                order_v=order_v_fitted,
                shape=np.array([num_control_points_u_fitted, num_control_points_v_fitted, 3]),
                control_points=np.array(cps_fitted).reshape((num_control_points_u_fitted*num_control_points_v_fitted,3)),
                knots_u=np.array(knot_vector_u),
                knots_v=np.array(knot_vector_v))
            #print(bspline_entity.name) 
            bspline_entity_surface.starting_geometry_index = self.num_points       
            self.input_bspline_entity_dict[bspline_entity.name] = bspline_entity_surface
            self.num_points += np.cumprod(bspline_entity_surface.shape)[-2]
        elif len(component_shape[:-1]) == 3:  # is volume
            print('fitting BSplineVolume has not been implemented yet')
            pass
        # return self.fit_bspline_entities(pts)
        return

    def read_file(self, file_name):
        if file_name != None:
            if ((file_name[-4:].lower() == '.stp') or (file_name[-5:].lower() == '.step')):
                with open(file_name, 'r') as f:
                    if 'CASCADE' in f.read():#Not sure, could be another string to identify
                        self.read_gmsh_stp(file_name)
                    else: 
                        self.read_openvsp_stp(file_name)
                
            elif ((file_name[-4:].lower() == '.igs') or (file_name[-5:].lower() == '.iges')):
                print('warning, not implemented yet!!')
                self.read_iges(file_name) #TODO
            else:
                print("Please input an iges file or a stp file from openvsp.")

        for bspline in self.input_bspline_entity_dict.values():
            self.output_geo_control_points_dict[f'{bspline.name}'] = np.reshape(bspline.control_points, (-1,3))

    def read_openvsp_stp(self, file_name):
        read_openvsp_stp(self, file_name)

    def read_gmsh_stp(self, file_name):
        read_gmsh_stp(self, file_name)

    def read_iges(self,file_name):
        read_iges(self, file_name)

    def write_step(self, file_name, plot=False):
        write_step(self, file_name, plot)

    def write_iges(self, file_name, plot = False):
        write_iges(self, file_name, plot)

    def generate_bspline_list(self, names_list):
        pointer_list = []
        for name in names_list:
            pointer_list.append(self.initial_input_bspline_entity_dict[name])
        return pointer_list

    def add_component(self, comp_obj):
        ''' 
        Extract control points, objects and names from the step file and import into the component object
        Uses these points to create a bounding box around the object to create a ffd box around it.
        Store component and FFD block into design geometry object.

        '''
        comp_obj._extract_embedded_entity_data(self)
        comp_obj._auto_generate_ffd()

        self.components[comp_obj.name] = comp_obj
        self.components_ffd_dict[comp_obj.name] = comp_obj.ffd

    def add_geometric_outputs(self, geo_outputs):
        self.geometric_outputs = geo_outputs

    def fit_bspline_ctrl_pointsets(self, fitting_pointset_list, output_vec = None, plot = False):
        output_pointset_list = [] 
        for pointset_fitted in fitting_pointset_list:
            component_shape = pointset_fitted.shape
            #print(pointset.pointset_id,pointset.name)
            if output_vec is None:
                entity_points_be_fitted = pointset_fitted.physical_coordinates
            else:
                start_index = pointset_fitted.output_starting_ind
                end_index = start_index + np.cumprod(component_shape)[-2]
                entity_points_be_fitted = output_vec[start_index:end_index]
            if len(component_shape[:-1]) == 0:  # is point  
                print('fitting points has not been implemented yet')
                pass        #is point
            elif len(component_shape[:-1]) == 1:  # is curve
                print('fitting curves has not been implemented yet')
                pass  
            elif len(component_shape[:-1]) == 2:  # is surface
                order_u_fitted = 4
                order_v_fitted = 4
                num_control_points_u_fitted = component_shape[0] - 1
                num_control_points_v_fitted = component_shape[1] - 1
                num_points_u_fitted = component_shape[0]
                num_points_v_fitted = component_shape[1]

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

                a = np.matmul(basis0.toarray().T, basis0.toarray())
                if np.linalg.det(a) == 0:
                    print('TODO: lstsq')
                    print('np.linalg.det(a) == 0', np.linalg.det(a) == 0)                      
                    cps_fitted = entity_points_be_fitted                        
                    ind = np.lexsort((cps_fitted[:,0],cps_fitted[:,1],cps_fitted[:,2])) 
                    print(ind.shape)
                    points_sorted = cps_fitted[ind]

                    order_u_fitted = 4
                    order_v_fitted = 4
                    num_control_points_u_fitted = 6
                    num_control_points_v_fitted = 180
                    num_points_u_fitted = component_shape[1]
                    num_points_v_fitted = component_shape[0]

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

                    a = np.matmul(basis0.toarray().T, basis0.toarray())
                    print('np.linalg.det(a) == 0', np.linalg.det(a) == 0)
                    print('num_control_points_u_fitted', num_control_points_u_fitted, 'num_control_points_v_fitted',num_control_points_v_fitted, 'num_points_u_fitted',num_points_u_fitted, 'num_points_v_fitted',num_points_v_fitted)
                    cps_fitted = np.linalg.solve(a, np.matmul(basis0.toarray().T, points_sorted))  
                    vp_points0 = vedo.Points(points_sorted, r=10, c='red',alpha=0.8)#cps_fitted[vector_indices,:]
                    vp_points1 = vedo.Points(cps_fitted, r=15, c='green',alpha=0.3)
                    vp_test = Plotter(axes=1)
                    vp_test.show(vp_points0,vp_points1, 'Test', viewup="z", interactive=True)  
                    exit()                      
                    #relative_map = sps.csc_matrix(relative_map)
                else: 
                    #print('num_control_points_u_fitted', num_control_points_u_fitted, 'num_control_points_v_fitted',num_control_points_v_fitted, 'num_points_u_fitted',num_points_u_fitted, 'num_points_v_fitted',num_points_v_fitted)
                    cps_fitted = np.linalg.solve(a, np.matmul(basis0.toarray().T, entity_points_be_fitted))  
                    relative_map = np.matmul(np.linalg.inv(a),basis0.toarray().T)
                relative_map = sps.csc_matrix(relative_map)
                pointset = PointSet(
                    pointset_id = self.current_id,
                    shape = np.array([num_control_points_u_fitted, num_control_points_v_fitted, 3]),#np.append(num_control_points_u_fitted * num_control_points_v_fitted,3),
                    output_starting_ind = self.current_pointset_pts,
                    parent_pointers = [pointset_fitted.pointset_id],
                    absolute_map = None,
                    relative_map = relative_map,
                    offset = np.array([0., 0., 0.]),            
                    offset_absolute_map = None,
                    physical_coordinates = None,
                    permutation_matrix = None,
                    name = 'ctrl_pts_' + pointset_fitted.name
                    )

                self.pointsets_dict[self.current_id] = pointset
                self.current_id += 1
                self.current_pointset_pts += np.cumprod(pointset.shape)[-2]  
                output_pointset_list.append(pointset)  

                if plot:#
                    self.assemble(pointset = pointset)
                    points = self.evaluate(pointset = pointset)
                    #points = relative_map.dot(entity_points_be_fitted)
                    vp_points0 = vedo.Points(entity_points_be_fitted, r=10, c='red',alpha=0.8)
                    vp_points1 = vedo.Points(cps_fitted, r=15, c='green',alpha=0.3)
                    vp_test = Plotter(axes=1)
                    vp_test.show(vp_points0,vp_points1, 'Test', viewup="z", interactive=True)

            else:  # is volume
                print('fitting BSplineVolume has not been implemented yet')
                pass 
        return output_pointset_list    
    
if __name__ == "__main__": 
    geo = DesignGeometry('CAD/eVTOL.stp')
