import numpy as np
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
# from lsdo_kit.design.design_geometry.bsplines.bspline_curve import BSplineCurve
# from lsdo_kit.design.design_geometry.bsplines.bspline_surface import BSplineSurface
# from lsdo_kit.design.design_geometry.bsplines.bspline_volume import BSplineVolume

import matplotlib.pyplot as plt

from vedo import Points, Plotter, LegendBox

import scipy.sparse as sps

import lsdo_geo.utils.helper_functions as hf

path_name = 'CAD/'
file_name = 'eVTOL.stp'
geo = DesignGeometry(path_name + file_name)

up_direction = np.array([0., 0., 1.])
down_direction = np.array([0., 0., -1.])

wing_surface_names = [
    'Surf_WFWKRQIMCA, Wing, 0, 12', 'Surf_WFWKRQIMCA, Wing, 0, 13', 
    'Surf_WFWKRQIMCA, Wing, 0, 14', 'Surf_WFWKRQIMCA, Wing, 0, 15', 
    ]
top_wing_surface_names = [ 
    'Surf_WFWKRQIMCA, Wing, 0, 13', 
    'Surf_WFWKRQIMCA, Wing, 0, 15', 
    ]
bot_wing_surface_names = [
    'Surf_WFWKRQIMCA, Wing, 0, 12',
    'Surf_WFWKRQIMCA, Wing, 0, 14', 
    ]

# ''' Points to be projected'''
# lead_root_point = np.array([170., 15., 130.])
# lead_tip_point = np.array([140., 215., 150.])
# trail_tip_point = np.array([190., 215., 150.])
# trail_root_point = np.array([230., 15., 130.])

# '''Project points'''
# wing_lead_root, wing_lead_root_coord = geo.project_point(lead_root_point, top_wing_surface_names, projection_direction = down_direction)
# wing_lead_tip, wing_lead_tip_coord = geo.project_point(lead_tip_point, top_wing_surface_names, projection_direction = down_direction)
# wing_trail_tip, wing_trail_tip_coord = geo.project_point(trail_tip_point, top_wing_surface_names, projection_direction = down_direction)
# wing_trail_root, wing_trail_root_coord = geo.project_point(trail_root_point, top_wing_surface_names, projection_direction = down_direction)

# '''Create chord surface using linear interpolation of lead/trail curves'''
# # num_points_along_edge = [7]
# num_points_along_chord = [10]
# rib_front_back = np.array([0.3,0.8])
# # surface_shape1 = np.append(num_points_along_edge,num_points_along_chord)
# # lead_edge_curve = geo.perform_linear_interpolation(wing_lead_root, wing_lead_tip, num_points_along_edge)
# # trail_edge_curve = geo.perform_linear_interpolation(wing_trail_root, wing_trail_tip, num_points_along_edge)

# output_parameters0 = np.linspace(rib_front_back[0],rib_front_back[1],num_points_along_chord[0])
# # print(output_parameters0)
# root_curve = geo.perform_linear_interpolation(wing_lead_root,wing_trail_root,num_points_along_chord, output_parameters = output_parameters0)
# tip_curve = geo.perform_linear_interpolation(wing_lead_tip,wing_trail_tip, num_points_along_chord, output_parameters = output_parameters0)

# # chord_surface = geo.perform_linear_interpolation(lead_edge_curve, trail_edge_curve, surface_shape1)

# # geo.assemble(chord_surface)
# # chord_surface_points = geo.evaluate(chord_surface)

# '''Create Rib'''

# # test1 = chord_surface_points[2*num_points_along_chord[0]:3*num_points_along_chord[0],:]
# # test2 = chord_surface_points[5*num_points_along_chord[0]:6*num_points_along_chord[0],:]

# # Method 1: Interpolate root and tip to rib_location along wing
# rib_location1 = 0.7
# shape1 = np.append(num_points_along_chord,1)
# output_parameters1 = rib_location1*np.ones([num_points_along_chord[0],1])
# chord_at_rib1 = geo.perform_linear_interpolation(root_curve,tip_curve,shape1, output_parameters = output_parameters1)

# geo.assemble(chord_at_rib1)
# care1 = geo.evaluate(chord_at_rib1)

# # #Method 2: Round rib_location up and down to get chord line indeces, interpolate between them
# # rib_location2 = 0.2345
# # rib_location_index = rib_location2 * num_points_along_edge[0]
# # if abs(rib_location_index-round(rib_location_index)) < 0.00000001:# What number do I compare it to
# #     chord_at_rib2 = chord_surface_points[int(rib_location_index*num_points_along_chord[0]):int((rib_location_index+1)*num_points_along_chord[0]),:]
# # else:
# #     index1 = int(np.floor(rib_location_index))
# #     index2 = int(np.ceil(rib_location_index))
# #     weight = rib_location_index-index1
# #     chord_line1 = chord_surface_points[index1*num_points_along_chord[0]:index2*num_points_along_chord[0],:]#np arrays, so can't perform linear interpolation
# #     chord_line2 = chord_surface_points[index2*num_points_along_chord[0]:(index2+1)*num_points_along_chord[0],:]
# #     print(type(chord_line1))
# #     print(type(chord_line2))
# #     output_parameters2 = weight*np.ones([num_points_along_chord[0],1])
# #     chord_at_rib2 = geo.perform_linear_interpolation(chord_line1,chord_line2,shape1,output_parameters = output_parameters2)
# #     print(type(chord_at_rib2))

# # geo.assemble(chord_at_rib2)
# # care2 = geo.evaluate(chord_at_rib2)

# '''Project from chord to wing surfaces'''
# top_points, top_plot = geo.project_surface(care1, top_wing_surface_names, up_direction)
# bot_points, bot_plot = geo.project_surface(care1, bot_wing_surface_names, down_direction)

# '''Interpolate between wing surfaces'''
# num_points_height = [5]
# surface_shape2 = np.append(shape1,num_points_height)
# rib_surface_mesh = geo.perform_linear_interpolation(top_points, bot_points, surface_shape2)

rib_surface_mesh = hf.create_element_from_chord(geo, top_wing_surface_names,bot_wing_surface_names, np.array([170., 15., 130.]),  np.array([140., 215., 150.]), np.array([190., 215., 150.]),np.array([230., 15., 130.]),
    shape = np.array([7,1,5]), loc_front_back=np.array([0.3, 0.2, 0.7]))

geo.register_output(rib_surface_mesh, name = 'Rib')
geo.assemble(rib_surface_mesh)
points_to_be_fitted = geo.evaluate(rib_surface_mesh)
# geo.fit_bspline_entities(points_to_be_fitted)
# print(points_to_be_fitted)

path_name = 'CAD/'
file_name = 'eVTOL_rib_from_chord.igs'
geo.write_iges(path_name + file_name)


'''Plotting'''
corner_projecting = []
corner_projected = []
vp_test = Plotter(N=1, axes=1)
top_bot = []
side = []
cps = []
surface_points = []
# chord_points = Points(chord_surface_points, r = 5, c = 'blue').legend('Chord_Surface')
# top_points_plot = Points(top_points, r = 12, c = 'plum').legend('test')
# bot_points_plot = Points(bot_points, r = 12, c = 'green').legend('test')
# care1_points = Points(care1, r = 12, c = 'brown').legend('Chord at Rib 1')
rib = Points(points_to_be_fitted, r = 12, c = 'black').legend('Ribs')
# care2_points = Points(care2, r = 12, c = 'violet').legend('Chord at Rib 2')
# corner_projecting.append(Points(np.vstack((lead_root_point,lead_tip_point,trail_tip_point ,trail_root_point)),
    # r=15, c='crimson').legend('Corner projecting'))
# corner_projected.append(Points(np.vstack((wing_lead_root_coord,wing_lead_tip_coord,wing_trail_tip_coord,wing_trail_root_coord)),
#     r=15, c='darkgoldenrod').legend('Corner projected'))
for target in wing_surface_names:
    cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))
    bspline_entity = geo.input_bspline_entity_dict[target]

    num_points_u = 50   # TODO might want to pass these in as input
    num_points_v = 50

    u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u), np.ones(num_points_v)).flatten()
    v_vec = np.einsum('i,j->ij', np.ones(num_points_u), np.linspace(0., 1., num_points_v)).flatten()

    pts = bspline_entity.evaluate_points(u_vec, v_vec)
    surface_points.append(Points(pts, r=5, c='seagreen',alpha=0.3).legend('Surface points'))
    cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))        
# for surf in geo.output_bspline_entity_dict.values():
#     if surf.name == 'Rib':
#         bspline_fitted_cps = Points(surf.control_points, r=9, c='plum').legend('Fitted bspline')


# vp_test.show(cps, surface_points, corner_projecting, corner_projected, 'Control points of surface to be projected', at=0, viewup="z", interactive = True)#, lb1
vp_test.show(cps, surface_points, rib, 'Control points of surface to be projected', at=0, viewup="z", interactive = True)#, lb1

