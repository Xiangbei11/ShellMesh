from lsdo_kit.cython.get_open_uniform_py import get_open_uniform
import numpy as np
import matplotlib.pyplot as plt
from vedo import Points, Plotter, LegendBox

from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 

import lsdo_kit
from lsdo_kit.design.design_geometry.bsplines.bspline_curve import BSplineCurve
from lsdo_kit.design.design_geometry.bsplines.bspline_surface import BSplineSurface
from lsdo_kit.design.design_geometry.bsplines.bspline_volume import BSplineVolume

from lsdo_kit.design.design_geometry.utils.calculate_rotation_mat import calculate_rotation_mat
# import os
# os.chdir("../lsdo_geo/lsdo_kit/design/design_geometry/core")

xglobal = np.array([1,0,0])
yglobal = np.array([0,1,0])
zglobal = np.array([0,0,1])

property_names_list = ['rot_x', 'rot_y', 'rot_z', 'trans_x', ' trans_y', 'trans_z', 'scale_y', 'scale_z', 'shape']

class FFD(object):
    def __init__(self, name: str, control_points, embedded_entities_pointers = [], local_axes=dict(xprime=None, yprime=None, zprime=None),
                    origin=None, local_control_points=None, shape=None, order_u=4, order_v=4, order_w=4,
                    knots_u=None, knots_v=None, knots_w=None):
        
        self.nxp = control_points.shape[0]
        self.nyp = control_points.shape[1]
        self.nzp = control_points.shape[2]
        
        self.embedded_entities_pointers = embedded_entities_pointers
        self.embedded_points = [] 
        self.embedded_entities_indices = []

        self.name = name

        self.properties_list = []
        self.properties_dict = {}
        
        if shape is None:
            shape = control_points.shape
        else:
            control_points = control_points.reshape(shape)

        if knots_u is None:
            knots_u = np.zeros(shape[0] + order_u)
            get_open_uniform(order_u, shape[0], knots_u)
        
        if knots_v is None:
            knots_v = np.zeros(shape[1] + order_v)
            get_open_uniform(order_v, shape[1], knots_v)

        if knots_w is None:
            knots_w = np.zeros(shape[2] + order_w)
            get_open_uniform(order_w, shape[2], knots_w)

        # if embedded_entities_pointers == []:
        #     pass
        # else:
        #     self._generate_embedded_coordinates()

        self.b_spline_volume = BSplineVolume(name, order_u, order_v, order_w, knots_u, knots_v, knots_w, shape, control_points)
        self.control_points = self.b_spline_volume.control_points

        self._generate_origin(origin)
        self._generate_local_coordinate_directions(local_axes)
        self._generate_local_control_points_coordinates()

    def add_shape_parameter(self, property_name: str, parameter_name: str, order: int , num_cp: int, dv: str, val = None):
        # if val is None:
        #     val = [None] * len(parameter_names)
        if property_name in property_names_list:
            if property_name not in self.properties_dict.keys():
                self.param_list = []
                self.param_dict = {}
                self.properties_dict[property_name] = {}
                self.properties_list.append(property_name)

                self.param_dict[parameter_name] = {'dv' : dv, 'val': val, 'order' : order, 'num' : num_cp}
                # self.param_list.append(self.param_dict) 
                # print('param_list first loop: ', self.param_list)
                self.properties_dict[property_name].update({'parameters': self.param_dict})
                # print('First loop:', self.properties_dict)
            
            elif property_name in self.properties_dict.keys():
                self.param_dict[parameter_name] = {'dv' : dv, 'val': val, 'order' : order, 'num' : num_cp}
                # print('param_dict: ', self.param_dict)
                # print('param_dict 1: ', self.param_dict['linear'])
                # print('param_dict 1: ', self.param_dict['quadratic'])
                
                # print('param_list before: ', self.param_list)
                # self.param_list.append(self.param_dict) 
                # print('param_list after: ', self.param_list)
                
                self.properties_dict[property_name].update({'parameters': self.param_dict})

                # print('properties_dict: ', self.properties_dict)
                # print('2nd loop:', self.param_list)
                # print('1st param_list: ', self.param_list[0])


        else:
            raise Exception('Must specify one of the following properties: rot_x, rot_y, rot_z, trans_x,  \
                trans_y, trans_z, scale_y, scale_z, shape')

    def _generate_midlines(self):
        u_vec = np.linspace(0, 1, self.nxp)

        v_vec_left = np.zeros(self.nxp)
        v_vec_mid  = np.ones(self.nxp) * 0.5
        v_vec_right = np.ones(self.nxp)
        
        w_vec_bot = np.zeros(self.nxp)
        w_vec_mid = np.ones(self.nxp) * 0.5
        w_vec_top = np.ones(self.nxp)

        left_mid_pt = self.b_spline_volume.evaluate_points(u_vec, v_vec_left, w_vec_mid)
        right_mid_pt = self.b_spline_volume.evaluate_points(u_vec, v_vec_right, w_vec_mid)
        bot_mid_pt = self.b_spline_volume.evaluate_points(u_vec, v_vec_mid, w_vec_bot)
        top_mid_pt = self.b_spline_volume.evaluate_points(u_vec, v_vec_mid, w_vec_top)

        self.top_bot_vec = top_mid_pt - bot_mid_pt
        self.right_left_vec = right_mid_pt - left_mid_pt

        #TODO: Discuss the direction of the midlines, which determines where the normal points 
    def _generate_origin(self, origin):
        if origin == None:

            u_vec = np.linspace(0, 1, self.control_points.shape[0])
            v_vec = np.ones(self.control_points.shape[0]) * 0.5
            w_vec = np.ones(self.control_points.shape[0]) * 0.5
        
            self.origin = self.b_spline_volume.evaluate_points(u_vec, v_vec, w_vec)
        
        elif origin.shape == (self.nxp, 3):
            self.origin = origin

        else:
            raise Exception('Origin shape must be equal to (nxp, 3)')


    def _generate_local_control_points_coordinates(self):
        self.local_control_points = np.zeros((self.nxp, self.nyp, self.nzp, 3))
        for i in range(self.nxp):

            self.local_control_points[i, :, :, :] = self.control_points[i, :, :, :] - self.origin[i, :]

            global_2_loc_rotmat_x = calculate_rotation_mat(xglobal, self.xprime[i, :])
            global_2_loc_rotmat_y = calculate_rotation_mat(yglobal, self.yprime[i, :])
            global_2_loc_rotmat_z = calculate_rotation_mat(zglobal, self.zprime[i, :])  

            # TODO: Change to MATMUL
            rot_mat = global_2_loc_rotmat_x * global_2_loc_rotmat_y * global_2_loc_rotmat_z

            self.initial_scale_y = np.max(self.local_control_points[i,:,:,1]) - np.min(self.local_control_points[i,:,:,1])
            self.initial_scale_z = np.max(self.local_control_points[i,:,:,2]) - np.min(self.local_control_points[i,:,:,2])

            self.local_control_points[i,:,:,:] = np.matmul(self.local_control_points[i,:,:,:], rot_mat)

            # Normalizing all of the y-coordinates
            self.local_control_points[i,:,:,1] = self.local_control_points[i,:,:,1] / self.initial_scale_y

            # Normalizing all of the z-coordinates
            self.local_control_points[i,:,:,2] = self.local_control_points[i,:,:,2] / self.initial_scale_z


    def _generate_local_coordinate_directions(self, local_axes):
        for k, v in local_axes.items():
            if k[0].lower() == 'x':
                if v is None:
                    self.xprime = np.tile(np.array([1,0,0]), (self.nxp,1))
                elif v.shape == (3,):
                    self.xprime = np.tile(v, (self.nxp,1))
                elif v.shape == (self.nxp, 3):
                    self.xprime = v
                else:
                    raise Exception('Xprime must either be a: (3,) array, (nxp, 3) array, or None')

            elif k[0].lower() == 'y':
                if v is None:
                    self.yprime = np.tile(np.array([0,1,0]), (self.nxp,1))
                elif v.shape == (3,):
                    self.yprime = np.tile(v, (self.nxp,1))
                elif v.shape == (self.nxp, 3):
                    self.yprime = v
                else:
                    raise Exception('Yprime must either be a: (3,) array, (nxp, 3) array, or None')
                    
            elif k[0].lower() == 'z':
                if v is None:
                    self.zprime = np.tile(np.array([0,0,1]), (self.nxp,1))
                elif v.shape == (3,):
                    self.zprime = np.tile(v, (self.nxp,1))
                elif v.shape == (self.nxp, 3):
                    self.zprime = v
                else:
                    raise Exception('Zprime must either be a: (3,) array, (nxp, 3) array, or None')
                    
            else:
                raise Exception('Keys need to begin with either an: x,y,z character to denote local axes!')


    def _generate_ffd_origin(self):
        # It doesn't really matter, so let's take the middle
        u_vec = 0.5
        v_vec = 0.5
        w_vec = 0.5
    
        self.ffd_origin = self.b_spline_volume.evaluate_points(u_vec, v_vec, w_vec)

    def _generate_ffd_axes(self):
        # It doesn't particularly matter, so let's ensure it's a normal basis by using the Earth-fixed frame
        self.ffd_xprime = np.array([1, 0, 0])
        self.ffd_yprime = np.array([0, 1, 0])
        self.ffd_zprime = np.array([0, 0, 1])


    def _generate_exterior_points(self, nu, nv, nw):

        v_vec_front, w_vec_front = np.mgrid[0:1:nv*1j, 0:1:nw*1j]
        u_vec_front = np.zeros(v_vec_front.shape) 

        v_vec_back, w_vec_back = np.mgrid[0:1:nv*1j, 0:1:nw*1j]
        u_vec_back = np.ones(v_vec_front.shape)

        u_vec_bot, v_vec_bot = np.mgrid[0:1:nu*1j, 0:1:nv*1j]
        w_vec_bot = np.zeros(v_vec_bot.shape)

        u_vec_top, v_vec_top = np.mgrid[0:1:nu*1j, 0:1:nv*1j]
        w_vec_top = np.ones(v_vec_bot.shape)

        u_vec_left, w_vec_left = np.mgrid[0:1:nu*1j, 0:1:nw*1j]
        v_vec_left = np.zeros(u_vec_left.shape)

        u_vec_right, w_vec_right = np.mgrid[0:1:nu*1j, 0:1:nw*1j]
        v_vec_right = np.ones(u_vec_right.shape)

        u_points = np.concatenate((u_vec_front.flatten(), u_vec_back.flatten(), u_vec_bot.flatten(), u_vec_top.flatten(), u_vec_left.flatten(), u_vec_right.flatten()))
        v_points = np.concatenate((v_vec_front.flatten(), v_vec_back.flatten(), v_vec_bot.flatten(), v_vec_top.flatten(), v_vec_left.flatten(), v_vec_right.flatten()))
        w_points = np.concatenate((w_vec_front.flatten(), w_vec_back.flatten(), w_vec_bot.flatten(), w_vec_top.flatten(), w_vec_left.flatten(), w_vec_right.flatten()))
        
        exterior_points = self.b_spline_volume.evaluate_points(u_points, v_points, w_points)

        return exterior_points
    
    def add_embedded_entities(self, entities):
        for i in entities:
            self.embedded_entities_pointers.append(i)

    def _generate_embedded_coordinates(self):

        for i in self.embedded_entities_pointers:
            # print('i: ', i)
            if isinstance(i, lsdo_kit.design.design_geometry.bsplines.bspline_curve.BSplineCurve) or isinstance(i, lsdo_kit.design.design_geometry.bsplines.bspline_surface.BSplineSurface) or isinstance(i, lsdo_kit.design.design_geometry.bsplines.bspline_volume.BSplineVolume):
                # print(i.control_points.shape)
                # print(len(i.control_points))
                self.embedded_points.append(i.control_points)  # Control points are given in Cartesian Coordinates
                self.embedded_entities_indices.append(len(i.control_points))
            else:
                self.embedded_points.append(i.physical_coordinates)  # Here i is a PointSet, make sure that the PointSet is evalauted
                self.embedded_entities_indices.append(len(i.physical_coordinates))

    def plot(self, nu, nv, nw):
        exterior_points = self._generate_exterior_points(nu, nv, nw)
        vp_init = Plotter()
        vps = []
        vps1 = Points(exterior_points, r=8, c = 'red')
        vps.append(vps1)       

        if self.embedded_entities_pointers == []:
            pass
        
        elif self.embedded_points == []:
            self._generate_embedded_coordinates()

        
        for i in self.embedded_points:
            vps2 =  Points(i, r=8, c='blue')
            vps.append(vps2)

        vp_init.show(vps, 'Bspline Volume', axes=1, viewup="z", interactive = True)

    def project_points_FFD(self):

        self._generate_embedded_coordinates()

        embedded_points = np.concatenate(self.embedded_points, axis=0 )
        
        self.ffd_application_map = self.b_spline_volume.compute_projection_eval_map(embedded_points)

        return self.ffd_application_map

if __name__ == "__main__":

    from lsdo_kit.geometry.utils.generate_ffd import create_ffd
    nxp = 5
    nyp = 5
    nzp = 5

    point000 = np.array([170. ,0. ,100.])
    point010 = np.array([130., 230., 100.])
    point001 = np.array([170., 0., 170.])
    point011 = np.array([130., 230., 170.])
    
    point100 = np.array([240. ,0. ,100.])
    point101 = np.array([240. ,0. ,170.])
    point110 = np.array([200. ,230. ,100.])
    point111 = np.array([200. ,230. ,170.])

    control_points = np.zeros((2,2,2,3))
    
    control_points[0,0,0,:] = point000
    control_points[0,0,1,:] = point001

    control_points[0,1,0,:] = point010
    control_points[0,1,1,:] = point011
    
    control_points[1,0,0,:] = point100
    control_points[1,0,1,:] = point101
    
    control_points[1,1,0,:] = point110
    control_points[1,1,1,:] = point111

    ffd_control_points = create_ffd(control_points, nxp, nyp, nzp)

    # ''' Camber surface creation script for this case '''
    # path_name = '../examples/CAD/'
    # file_name = 'eVTOL.stp'
    # geo = DesignGeometry(path_name + file_name)

    # wing_surface_names = [
    # 'Surf_WFWKRQIMCA, Wing, 0, 12', 'Surf_WFWKRQIMCA, Wing, 0, 13', 
    # 'Surf_WFWKRQIMCA, Wing, 0, 14', 'Surf_WFWKRQIMCA, Wing, 0, 15', 
    # ]

    # bspline_entities = [geo.input_bspline_entity_dict[wing_surface_names[0]],
    #    geo.input_bspline_entity_dict[wing_surface_names[1]], 
    #    geo.input_bspline_entity_dict[wing_surface_names[2]],
    #    geo.input_bspline_entity_dict[wing_surface_names[3]]]


    # local_axes = {'xprime': np.array([1,0,0]), 'yprime': np.array([0,1,0]), 'zprime': np.array([0,0,1]) }

    # test_ffd = FFD('test', ffd_control_points, embedded_entities_pointers=bspline_entities)
    test_ffd = FFD('test', ffd_control_points)
    
    test_ffd.add_shape_parameter('rot_x', 'linear', 2, 3, False, val=1.0)
    test_ffd.add_shape_parameter('rot_x', 'quadratic', 3, 4, False, val=1.0)

    test_ffd.add_shape_parameter('rot_y', 'linear', 3,4, True, val=2.0)

    # print('full dict: ', test_ffd.properties_dict)
    # print('\n')
    # print('rot_x: ', test_ffd.properties_dict['rot_x']['parameters'])   
    # print('\n')
    # print('rot_y: ',test_ffd.properties_dict['rot_y']['parameters'])
    # print('\n')

    # print('properties_list: ', test_ffd.properties_dict.keys())
    for property_name, property_dict in test_ffd.properties_dict.items():
        print('property_name: ', property_name)
        print('property_dict: ', property_dict)
        
        # property_var = ffd_block_csdl_model.create_input(property_name)
        # parameters_list = test_ffd.properties_dict[property_name]
        # print('parameters_list: ', parameters_list)
        # print(property_dict['parameters'])
        for parameter_name, parameter_info in property_dict['parameters'].items():
                print('parameter_name: ', parameter_name)
                # for parameter_info in parameter_dict():
                    # print('parameter_name: ', parameter_name)
                print('parameter_info: ', parameter_info)
    
    # print(test_ffd.properties_dict['rot_x']['parameters'][1])    
    

    # print(test_ffd.embedded_entities_indices)

    # print('ORIGINAL CONTROL PTS: ', test_ffd.control_points[0,:,:,:])
    # print('ORIGIN: ', test_ffd.origin[0,:])
    # print('LOCAL CONTROL PTS: ', test_ffd.control_points[0,:,:,:] - test_ffd.origin[0,:])

    # print(test_ffd.local_control_points)
    # print(test_ffd.BSplineVolume.control_points)

    test_ffd.plot(nxp, nyp, nzp)

    # test_ffd.translate_control_points(offset=np.array([10., 50., 100.]))

    # print(test_ffd.BSplineVolume.control_points)

    # test_ffd.plot(nu, nv, nw)


