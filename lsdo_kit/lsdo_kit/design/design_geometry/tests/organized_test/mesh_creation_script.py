import csdl
import numpy as np
from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
from lsdo_kit.old_files.mesh import Mesh
from lsdo_kit.design.design_geometry.core.component import Component

# from lsdo_kit.tests.organized_test.test_inner_opt import geo

def generate_meshes(geo):

    ''' Creating the chamber line mesh for mesh eval '''

    ''' Using points from rect_wing 10.17.2021 '''

    # top_wing_surface_names = [
    #     'Surf_QRMTIHUVEP, RectWing, 0, 3', 
    #     'Surf_QRMTIHUVEP, RectWing, 1, 9', 
    #     ]

    # bot_wing_surface_names = [
    #     'Surf_QRMTIHUVEP, RectWing, 0, 2',
    #     'Surf_QRMTIHUVEP, RectWing, 1, 8', 
    #     ]

    top_wing_surface_names = [
        'RectWing, 0, 3', 
        'RectWing, 1, 9', 
        ]

    bot_wing_surface_names = [
        'RectWing, 0, 2',
        'RectWing, 1, 8', 
        ]

    # wing_comp = Component(stp_entity_names=['RectWing'], name='wing')
    # # c._add_GeometryPrimitve_to_vehicle()
    # geo.add_component(wing_comp)

    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    left_lead_point = np.array([0., -9000., 2000.])/1000
    left_trail_point = np.array([4000.0, -9000.0, 2000.])/1000
    right_lead_point = np.array([0.0, 9000.0, 2000.])/1000
    right_trail_point = np.array([4000.0, 9000.0, 2000.])/1000

    '''Project points'''
    wing_lead_left, wing_lead_left_coord = geo.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_left ,wing_trail_left_coord = geo.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_lead_right ,wing_lead_right_coord = geo.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_right ,wing_trail_right_coord = geo.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

    print('WING LEAD LEFT COORD: ', wing_lead_left_coord)
    print('WING LEAD RIGHT COORD: ', wing_lead_right_coord)


    num_pts2 = [15]#20
    num_pts1 = [3]#5
    surface_shape1 = np.append(num_pts1,num_pts2)

    left_chord = geo.perform_linear_interpolation(wing_lead_left, wing_trail_left, num_pts1)
    right_chord = geo.perform_linear_interpolation(wing_lead_right, wing_trail_right, num_pts1)

    # lead_edge_curve = geo.perform_linear_interpolation(wing_lead_left, wing_lead_right, num_pts1)
    # trail_edge_curve = geo.perform_linear_interpolation(wing_trail_left, wing_trail_right, num_pts1)

    chord_surface = geo.perform_linear_interpolation(left_chord, right_chord, surface_shape1)
    geo.assemble(chord_surface)
    chord_surface_points = geo.evaluate(chord_surface)

    '''Project from chord to wing surfaces'''
    top_points, top_plot = geo.project_points(chord_surface_points, projection_direction=up_direction, projection_targets_names=top_wing_surface_names)
    bot_points, bot_plot = geo.project_points(chord_surface_points, projection_direction=down_direction, projection_targets_names=bot_wing_surface_names)

    '''Create Camber Surface'''
    surface_shape2 = np.append(surface_shape1,1)
    output_parameters = 0.5*np.ones([num_pts1[0]*num_pts2[0],1])
    camber_surface_mesh = geo.perform_linear_interpolation(top_points, bot_points, surface_shape2, output_parameters = output_parameters)
    camber_surface_mesh.name = 'camber_surface_mesh'

    camber_surface_mesh.reshape((num_pts1[0],num_pts2[0],3))
    print("CAMBER SURFACE MESH SHAPE: ", camber_surface_mesh.shape)
    return left_chord, right_chord, camber_surface_mesh