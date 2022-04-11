import numpy as np
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 

import matplotlib.pyplot as plt

from vedo import Points, Plotter, LegendBox

from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
import scipy.sparse as sps

''' Camber surface creation script '''
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

''' Points to be projected'''
lead_root_point = np.array([170., 15., 130.])
lead_tip_point = np.array([140., 215., 150.])
trail_tip_point = np.array([190., 215., 150.])
trail_root_point = np.array([230., 15., 130.])

'''Project points'''
wing_lead_root, wing_lead_root_coord = geo.project_points(lead_root_point, top_wing_surface_names, projection_direction = down_direction)
wing_lead_tip ,wing_lead_tip_coord = geo.project_points(lead_tip_point, top_wing_surface_names, projection_direction = down_direction)
wing_trail_tip ,wing_trail_tip_coord = geo.project_points(trail_tip_point, top_wing_surface_names, projection_direction = down_direction)
wing_trail_root ,wing_trail_root_coord = geo.project_points(trail_root_point, top_wing_surface_names, projection_direction = down_direction)

print(type(wing_lead_root))


'''Create chord surface using linear interpolation of lead/trail curves'''
num_pts1 = [20]#20
num_pts2 = [20]#5
surface_shape1 = np.append(num_pts1,num_pts2)
lead_edge_curve = geo.perform_linear_interpolation(wing_lead_root, wing_lead_tip, num_pts1)
trail_edge_curve = geo.perform_linear_interpolation(wing_trail_root, wing_trail_tip, num_pts1)

chord_surface = geo.perform_linear_interpolation(lead_edge_curve, trail_edge_curve, surface_shape1)

'''Create chord surface using transfinite interpolation'''
# root_curve = geo.perform_linear_interpolation(wing_lead_root, wing_trail_root, num_pts1)
# tip_curve = geo.perform_linear_interpolation(wing_lead_tip, wing_trail_tip, num_pts1)
# chord_surf = geo.perform_2d_transfinite_interpolation(lead_edge_curve, trail_edge_curve, root_curve, tip_curve)

'''Find absolute map of chord surface'''
geo.assemble(chord_surface)
chord_surface_mesh = geo.evaluate(chord_surface)
# print(chord_surface_mesh)

'''Translate chord surface up and down'''
# chord_surface_trans_up = np.zeros([np.shape(chord_surface_mesh)[0],3])
# chord_surface_trans_up[:,:2] = chord_surface_mesh[:,:2]
# chord_surface_trans_up[:,2]= chord_surface_mesh[:,2]+8

# chord_surface_trans_down = np.zeros([np.shape(chord_surface_mesh)[0],3])
# chord_surface_trans_down[:,:2] = chord_surface_mesh[:,:2]
# chord_surface_trans_down[:,2]= chord_surface_mesh[:,2]-8

'''Project from translated surfaces to wing surfaces'''
# top_points, top_plot = geo.project_surface(chord_surface_trans_up, top_wing_surface_names, down_direction)
# bot_points, bot_plot = geo.project_surface(chord_surface_trans_down, bot_wing_surface_names, up_direction)

'''Project from chord to wing surfaces'''
top_points, top_plot = geo.project_points(chord_surface_mesh, top_wing_surface_names, up_direction)
bot_points, bot_plot = geo.project_points(chord_surface_mesh, bot_wing_surface_names, down_direction)

'''Create Camber Surface'''
surface_shape2 = np.append(surface_shape1,1)
output_parameters = 0.5*np.ones([num_pts1[0]*num_pts2[0],1])
camber_surface_mesh = geo.perform_linear_interpolation(top_points, bot_points, surface_shape2, output_parameters = output_parameters)



# num_interpolations1 = 20#20
# num_interpolations2 = 20#5

# lead_edge_curve = geo.perform_linear_interpolation(wing_lead_root, wing_lead_tip, num_interpolations1)
# trail_edge_curve = geo.perform_linear_interpolation(wing_trail_root, wing_trail_tip, num_interpolations1)

# chord_surface = geo.perform_linear_interpolation(lead_edge_curve, trail_edge_curve, num_interpolations2)

# '''Create chord surface using transfinite interpolation'''
# # root_curve = geo.perform_linear_interpolation(wing_lead_root, wing_trail_root, num_pts1)
# # tip_curve = geo.perform_linear_interpolation(wing_lead_tip, wing_trail_tip, num_pts1)
# # chord_surf = geo.perform_2d_transfinite_interpolation(lead_edge_curve, trail_edge_curve, root_curve, tip_curve)

# '''Find absolute map of chord surface'''
# geo.assemble(chord_surface)
# chord_surface_mesh = geo.evaluate(chord_surface)
# # print(chord_surface_mesh)

# '''Translate chord surface up and down'''
# # chord_surface_trans_up = np.zeros([np.shape(chord_surface_mesh)[0],3])
# # chord_surface_trans_up[:,:2] = chord_surface_mesh[:,:2]
# # chord_surface_trans_up[:,2]= chord_surface_mesh[:,2]+8

# # chord_surface_trans_down = np.zeros([np.shape(chord_surface_mesh)[0],3])
# # chord_surface_trans_down[:,:2] = chord_surface_mesh[:,:2]
# # chord_surface_trans_down[:,2]= chord_surface_mesh[:,2]-8

# '''Project from translated surfaces to wing surfaces'''
# # top_points, top_plot = geo.project_surface(chord_surface_trans_up, top_wing_surface_names, down_direction)
# # bot_points, bot_plot = geo.project_surface(chord_surface_trans_down, bot_wing_surface_names, up_direction)

# '''Project from chord to wing surfaces'''
# top_points, top_plot = geo.project_surface(chord_surface_mesh, top_wing_surface_names, up_direction)
# bot_points, bot_plot = geo.project_surface(chord_surface_mesh, bot_wing_surface_names, down_direction)

# '''Create Camber Surface'''
# num_interpolations3 = 1

# output_parameters = 0.5*np.ones([num_interpolations1*num_interpolations2,1])
# camber_surface_mesh = geo.perform_linear_interpolation(top_points, bot_points, num_interpolations3, output_parameters = output_parameters)





geo.register_output(camber_surface_mesh, name = 'Camber Surface')
geo.assemble(camber_surface_mesh)

'''Evaluate the physical coornadites of points to be fitted'''
camber_surface = geo.evaluate(camber_surface_mesh)
geo.fit_bspline_entities([camber_surface_mesh])
geo.write_iges('CAD/eVTOL_camber.igs')


'''Plot'''
corner_projecting = []
corner_projected = []
vp_test = Plotter(N=1, axes=1)
top_bot = []
side = []
cps = []
surface_points = []
chord_points = Points(chord_surface_mesh, r = 7, c = 'plum').legend('Chord_Surface')
#up_chord_points = Points(chord_surface_trans_up, r = 5, c = 'green').legend('Chord Translated Up')
#down_chord_points = Points(chord_surface_trans_down, r = 5, c = 'green').legend('Chord Translated Down')
top_surf_points = Points(top_points.physical_coordinates, r = 7, c = 'orange').legend('Top Surface Projection')
camber_points = Points(camber_surface, r = 7, c = 'blue').legend('Camber_Surface')
projected_surfaces = Points(np.vstack((top_plot, bot_plot)), r = 15, c = 'pink').legend('Top and Bottom')
corner_projecting.append(Points(np.vstack((lead_root_point,lead_tip_point,trail_tip_point ,trail_root_point)),
    r=15, c='crimson').legend('Corner projecting'))
corner_projected.append(Points(np.vstack((wing_lead_root_coord,wing_lead_tip_coord,wing_trail_tip_coord,wing_trail_root_coord)),
    r=15, c='darkgoldenrod').legend('Corner projected'))
for target in wing_surface_names:
    cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))
    bspline_entity = geo.input_bspline_entity_dict[target]
    order_u = bspline_entity.order_u
    order_v = bspline_entity.order_v
    num_control_points_u = bspline_entity.shape[0]
    num_control_points_v = bspline_entity.shape[1]
    num_points_u = 50   # TODO might want to pass these in as input
    num_points_v = 50

    nnz = num_points_u * num_points_v * order_u * order_v
    data = np.zeros(nnz)
    row_indices = np.zeros(nnz, np.int32)
    col_indices = np.zeros(nnz, np.int32)

    knot_vector_u = bspline_entity.knots_u
    knot_vector_v = bspline_entity.knots_u
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
    surface_points.append(Points(pts, r=5, c='seagreen',alpha=0.3).legend('Surface points'))
    cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))        
# for surf in geo.output_bspline_entity_dict.values():
#     if surf.name == 'Rib':
#         bspline_fitted_cps = Points(surf.control_points, r=9, c='plum').legend('Fitted bspline')


# vp_test.show(cps, surface_points, corner_projecting, corner_projected, 'Control points of surface to be projected', at=0, viewup="z", interactive = True)#, lb1
vp_test.show(cps, surface_points, corner_projecting, corner_projected, chord_points, camber_points, 'Control points of surface to be projected', at=0, viewup="z", interactive = True)#, lb1
