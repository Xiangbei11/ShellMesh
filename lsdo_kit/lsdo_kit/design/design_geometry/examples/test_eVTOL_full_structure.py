import numpy as np
import pyiges

from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry
from lsdo_kit.geometry.utils.helper_functions import create_rib_or_spar, create_uniform_ribs_and_spars

import matplotlib.pyplot as plt
import time

from vedo import Points, Plotter, LegendBox

from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
import scipy.sparse as sps

time_start = time.time()

''' Spars and ribs creation script '''
path_name = 'CAD/'#lsdo_geo/DesignGeometry/
file_name = 'eVTOL.stp'
geo = DesignGeometry(path_name + file_name)

'''
Wing structures
'''
print('Starting wing structures')
wing_surface_names_sb = [
    'Surf_WFWKRQIMCA, Wing, 0, 12',
    'Surf_WFWKRQIMCA, Wing, 0, 13',
    'Surf_WFWKRQIMCA, Wing, 0, 14',
    'Surf_WFWKRQIMCA, Wing, 0, 15',
    ]
top_wing_surface_names_sb = [ 
    'Surf_WFWKRQIMCA, Wing, 0, 13', 
    'Surf_WFWKRQIMCA, Wing, 0, 15', 
    ]
bot_wing_surface_names_sb = [
    'Surf_WFWKRQIMCA, Wing, 0, 12',
    'Surf_WFWKRQIMCA, Wing, 0, 14', 
    ]

wing_surf_names_po = []
# wing_surf_names_mid = []
for surf_name in geo.input_bspline_entity_dict.keys():
    if ', Wing' in surf_name and surf_name not in wing_surface_names_sb:
        wing_surf_names_po.append(surf_name)

# print(wing_surf_names_po)

down_direction = np.array([0., 0., -1.])
up_direction = np.array([0., 0., 1.])
# num_ribs = 5

# Project the points onto the surface
top_leading_root_sb = np.array([185.,20.81,160.])
top_trailing_root_sb = np.array([215.,20.81,160.])
top_leading_tip_sb = np.array([158.,190.,160.])
top_trailing_tip_sb = np.array([178.,190.,160.])

bot_leading_root_sb = np.array([185.,20.81,110.])
bot_trailing_root_sb = np.array([215.,20.81,110.])
bot_leading_tip_sb = np.array([158.,190.,110.])
bot_trailing_tip_sb = np.array([178.,190.,110.])

top_leading_root_po = np.array([185.,-20.81,160.])
top_trailing_root_po = np.array([215.,-20.81,160.])
top_leading_tip_po = np.array([158.,-190.,160.])
top_trailing_tip_po = np.array([178.,-190.,160.])

bot_leading_root_po = np.array([185.,-20.81,110.])
bot_trailing_root_po = np.array([215.,-20.81,110.])
bot_leading_tip_po = np.array([158.,-190.,110.])
bot_trailing_tip_po = np.array([178.,-190.,110.])

eps1 = 5e-1

num_spars = 3
num_ribs = 5

num_pts_ribs_top_bot = 7
num_side_pts_spars_and_ribs = 5

# Starboard wing structures
spars, ribs = create_uniform_ribs_and_spars(geo, wing_surface_names_sb, bot_leading_root_sb, bot_leading_tip_sb, top_leading_root_sb, top_leading_tip_sb,
        bot_trailing_root_sb, bot_trailing_tip_sb, top_trailing_root_sb, top_trailing_tip_sb,
        np.array([num_spars, num_ribs]), np.array([num_pts_ribs_top_bot, num_side_pts_spars_and_ribs]), epsilon=eps1)
for i in range(len(spars)):
    geo.register_output(spars[i], name = f'sb_spar{i}')
for i in range(len(ribs)):
    geo.register_output(ribs[i], name = f'sb_rib{i}')

# Port wing structures
spars, ribs = create_uniform_ribs_and_spars(geo, wing_surf_names_po, bot_leading_root_po, bot_leading_tip_po, top_leading_root_po, top_leading_tip_po,
        bot_trailing_root_po, bot_trailing_tip_po, top_trailing_root_po, top_trailing_tip_po,
        np.array([num_spars, num_ribs]), np.array([num_pts_ribs_top_bot, num_side_pts_spars_and_ribs]), epsilon=eps1)
for i in range(len(spars)):
    geo.register_output(spars[i], name = f'po_spar{i}')
for i in range(len(ribs)):
    geo.register_output(ribs[i], name = f'po_rib{i}')

# Middle wing structures
spars, ribs = create_uniform_ribs_and_spars(geo, wing_surf_names_po, bot_leading_root_po, bot_leading_root_sb, top_leading_root_po, top_leading_root_sb,
        bot_trailing_root_po, bot_trailing_root_sb, top_trailing_root_po, top_trailing_root_sb,
        np.array([num_spars, num_ribs]), np.array([num_pts_ribs_top_bot, num_side_pts_spars_and_ribs]), epsilon=eps1)
for i in range(len(spars)):
    geo.register_output(spars[i], name = f'mid_spar{i}')
for i in range(len(ribs)):
    geo.register_output(ribs[i], name = f'mid_rib{i}')

'''
Tail structures
'''
#print('Starting tail structures')
tail_surf_names = []
for surf_name in geo.input_bspline_entity_dict.keys():
    if ', Tail' in surf_name:
        tail_surf_names.append(surf_name)
#print(tail_surf_names)

# tail_surface_names = [
#     'Surf_WFWKRQIMCA, Wing, 0, 12',
#     'Surf_WFWKRQIMCA, Wing, 0, 13',
#     'Surf_WFWKRQIMCA, Wing, 0, 14',
#     'Surf_WFWKRQIMCA, Wing, 0, 15',
#     ]
# tail_wing_surface_names = [ 
#     'Surf_WFWKRQIMCA, Wing, 0, 13', 
#     'Surf_WFWKRQIMCA, Wing, 0, 15', 
#     ]
# tail_wing_surface_names = [
#     'Surf_WFWKRQIMCA, Wing, 0, 12',
#     'Surf_WFWKRQIMCA, Wing, 0, 14', 
#     ]
down_direction = np.array([0., 0., -1.])
up_direction = np.array([0., 0., 1.])
num_tail_ribs = 5

# Project the points onto the surface
top_leading_root_sb = np.array([555.,0.,200.])
top_trailing_root_sb = np.array([580.,0.,200.])
top_leading_tip_sb = np.array([565.,90.,200.])
top_trailing_tip_sb = np.array([580.,90.,200.])

bot_leading_root_sb = np.array([555.,0.,120.])
bot_trailing_root_sb = np.array([580.,0.,120.])
bot_leading_tip_sb = np.array([565.,90.,120.])
bot_trailing_tip_sb = np.array([580.,90.,120.])

top_leading_root_po = np.array([555.,0.,200.])
top_trailing_root_po = np.array([580.,0.,200.])
top_leading_tip_po = np.array([565.,-90.,200.])
top_trailing_tip_po = np.array([580.,-90.,200.])

bot_leading_root_po = np.array([555.,0.,120.])
bot_trailing_root_po = np.array([580.,0.,120.])
bot_leading_tip_po = np.array([565.,-90.,120.])
bot_trailing_tip_po = np.array([580.,-90.,120.])

eps1 = 1e-8

num_spars = 2
num_ribs = 4

num_pts_ribs_top_bot = 7
num_side_pts_spars_and_ribs = 4

# spars, ribs = create_uniform_ribs_and_spars(geo, tail_surf_names, bot_leading_root_sb, bot_leading_tip_sb, top_leading_root_sb, top_leading_tip_sb,
#         bot_trailing_root_sb, bot_trailing_tip_sb, top_trailing_root_sb, top_trailing_tip_sb,
#         np.array([num_spars, num_ribs]), np.array([num_pts_ribs_top_bot, num_side_pts_spars_and_ribs]), epsilon=eps1)
# for i in range(len(spars)):
#     geo.register_output(spars[i], name = f'tail_sb_spar{i}')
# for i in range(len(ribs)):
#     geo.register_output(ribs[i], name = f'tail_sb_rib{i}')

# spars, ribs = create_uniform_ribs_and_spars(geo, tail_surf_names, bot_leading_root_po, bot_leading_tip_po, top_leading_root_po, top_leading_tip_po,
#         bot_trailing_root_po, bot_trailing_tip_po, top_trailing_root_po, top_trailing_tip_po,
#         np.array([num_spars, num_ribs]), np.array([num_pts_ribs_top_bot, num_side_pts_spars_and_ribs]), epsilon=eps1)
# for i in range(len(spars)):
#     geo.register_output(spars[i], name = f'tail_po_spar{i}')
# for i in range(len(ribs)):
#     geo.register_output(ribs[i], name = f'tail_po_rib{i}')

# # Construct upper and lower curves of primary spar
# primary_spar_top = np.linspace(top_leading_root, top_leading_tip, num_tail_ribs)
# primary_spar_top_curve, test_bot_pri_spar = geo.project_curve(primary_spar_top, projection_targets_names=tail_surf_names, offset=np.array([0., 0., eps1]), projection_direction = down_direction)

# primary_spar_bot = np.linspace(bot_leading_root, bot_leading_tip, num_tail_ribs)
# primary_spar_bot_curve, test_bot_pri_spar = geo.project_curve(primary_spar_bot, projection_targets_names=tail_surf_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)

# # Construct upper and lower curves of rear spar
# rear_spar_top = np.linspace(top_trailing_root, top_trailing_tip, num_tail_ribs)
# rear_spar_top_curve, test_bot_pri_spar = geo.project_curve(rear_spar_top, projection_targets_names=tail_surf_names, offset=np.array([0., 0., eps1]), projection_direction = down_direction)

# rear_spar_bot = np.linspace(bot_trailing_root, bot_trailing_tip, num_tail_ribs)
# rear_spar_bot_curve, test_bot_pri_spar = geo.project_curve(rear_spar_bot, projection_targets_names=tail_surf_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)


# # Construct side curves of primary spar
# num_side_pts_tail_spars_and_tail_ribs = 5

# top_projected_point1 = geo.extract_pointset(primary_spar_top_curve, np.array([0]), np.array([1]))
# top_projected_point2 = geo.extract_pointset(primary_spar_top_curve, np.array([-1]), np.array([1]))
# bot_projected_point1 = geo.extract_pointset(primary_spar_bot_curve, np.array([0]), np.array([1]))
# bot_projected_point2 = geo.extract_pointset(primary_spar_bot_curve, np.array([-1]), np.array([1]))
# primary_spar_side_root_curve = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, np.array([num_side_pts_tail_spars_and_tail_ribs]))
# primary_spar_side_tip_curve = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, np.array([num_side_pts_tail_spars_and_tail_ribs]))
# # Construct primary spar surface
# primary_spar_surface_mesh = geo.perform_2d_transfinite_interpolation(primary_spar_top_curve, primary_spar_bot_curve, 
#     primary_spar_side_root_curve, primary_spar_side_tip_curve)
# geo.register_output(primary_spar_surface_mesh, name = 'primary_tail_spar')


# # Construct side curves of rear spar
# top_projected_point1 = geo.extract_pointset(rear_spar_top_curve, np.array([0]), np.array([1]))
# top_projected_point2 = geo.extract_pointset(rear_spar_top_curve, np.array([-1]), np.array([1]))
# bot_projected_point1 = geo.extract_pointset(rear_spar_bot_curve, np.array([0]), np.array([1]))
# bot_projected_point2 = geo.extract_pointset(rear_spar_bot_curve, np.array([-1]), np.array([1]))
# rear_spar_side_root_curve = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, np.array([num_side_pts_tail_spars_and_tail_ribs]))
# rear_spar_side_tip_curve = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, np.array([num_side_pts_tail_spars_and_tail_ribs]))
# # Construct rear spar surface
# rear_spar_surface_mesh = geo.perform_2d_transfinite_interpolation(rear_spar_top_curve, rear_spar_bot_curve, 
#     rear_spar_side_root_curve, rear_spar_side_tip_curve)
# geo.register_output(rear_spar_surface_mesh, name = 'rear_tail_spar')

# num_pts_tail_ribs_top_bot = 7
# eps2 = 0

# for i in range(num_tail_ribs):
#     # Construct upper and lower curves of tail_ribs
#     rib_top = np.linspace(primary_spar_top[i]-eps2, rear_spar_top[i]+eps2, num_pts_tail_ribs_top_bot)
#     rib_top_curve, test_top_rib = geo.project_curve(rib_top, projection_targets_names=tail_surf_names, offset=np.array([0., 0., eps1]),
#         projection_direction = down_direction, plot=True)
#     rib_bot = np.linspace(primary_spar_bot[i]-eps2, rear_spar_bot[i]+eps2, num_pts_tail_ribs_top_bot)
#     rib_bot_curve, test_bot_rib = geo.project_curve(rib_bot, projection_targets_names=tail_surf_names, offset=np.array([0., 0., -eps1]),
#         projection_direction = up_direction, plot=True)

#     # Construct side curves of tail_ribs
#     num_side_pts_tail_spars_and_tail_ribs = 10
#     top_leading_root_pt = geo.extract_pointset(rib_top_curve, np.array([0]), np.array([1]))
#     bot_leading_root_pt = geo.extract_pointset(rib_bot_curve, np.array([0]), np.array([1]))
#     rib_side_root_curve = geo.perform_linear_interpolation(top_leading_root_pt, bot_leading_root_pt, np.array([num_side_pts_tail_spars_and_tail_ribs]))
#     top_leading_tip_pt = geo.extract_pointset(rib_top_curve, np.array([-1]), np.array([1]))
#     bot_leading_tip_pt = geo.extract_pointset(rib_bot_curve, np.array([-1]), np.array([1]))
#     rib_side_tip_curve = geo.perform_linear_interpolation(top_leading_tip_pt, bot_leading_tip_pt, np.array([num_side_pts_tail_spars_and_tail_ribs]))
#     # Construct rib surface
#     rib_surface_mesh = geo.perform_2d_transfinite_interpolation(rib_top_curve, rib_bot_curve, 
#         rib_side_root_curve, rib_side_tip_curve)
#     geo.register_output(rib_surface_mesh, name = f'tail_rib{i}')


'''
Fuselage structures
'''
#print('Starting fuselage structures')



# Concatenates vertically all the linear matrices 
geo.assemble()
# Evaluate the physical coornadites of points to be fitted
points_to_be_fitted = geo.evaluate()
time_end = time.time()
runtime = time_end - time_start
print('runtime', runtime)

# print(points_to_be_fitted)

geo.fit_bspline_entities(points_to_be_fitted)
file_name = 'eVTOL_full_structure.igs'
geo.write_iges(path_name + file_name)

#corner_projecting = []
#corner_projected = []

#top_bot = []
#side = []
cps = []
bspline_fitted_cps = []
surface_points = []
#corner_projecting.append(Points(np.vstack(([180.,60.,145.],[210.,60.,145.],[180.,60.,100.],[210.,60.,100.])), r=15, c='crimson').legend('Corner projecting'))
#corner_projected.append(Points(np.vstack((test_top_curve[0,:],test_top_curve[-1,:],test_bot_curve[0,:],test_bot_curve[-1,:])), r=15, c='darkgoldenrod').legend('Corner projected'))
#top_bot.append(Points(test_top_curve, r=15, c='brown').legend('Top curve'))
#top_bot.append(Points(test_bot_curve, r=15, c='brown').legend('Bot curve'))
#side.append(Points(np.linspace(test_top_curve[0,:],test_bot_curve[0,:],num_pts2), r=15, c='chartreuse').legend('Side curve'))
#side.append(Points(np.linspace(test_top_curve[-1,:],test_bot_curve[-1,:],num_pts2), r=15, c='chartreuse').legend('Side curve'))
TFI = Points(points_to_be_fitted, r=10, c='slategray')
# for target in wing_surface_names:
for target in tail_surf_names:
#for target in geo.input_bspline_entity_dict.keys():
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
    surface_points.append(Points(pts, r=5, c='seagreen',alpha=0.3).legend('Surface points'))
    cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))          
for surf in geo.output_bspline_entity_dict.values():
    # if surf.name == 'primary_spar':
    #     bspline_fitted_cps.append(Points(surf.control_points, r=12, c='black').legend('Fitted bspline'))
    # if  surf.name == 'rear_spar':
    #     bspline_fitted_cps.append(Points(surf.control_points, r=12, c='yellow').legend('Fitted bspline'))
    if  'rib' in surf.name:
        bspline_fitted_cps.append(Points(surf.control_points, r=9, c='red').legend('Fitted bspline'))
    # if  surf.name == 'middle_spar':
    #     bspline_fitted_cps.append(Points(surf.control_points, r=12, c='blue').legend('Fitted bspline'))
#vp_test.show(cps, 'Control points of surface to be projected', at=0, viewup="z")#, lb1
#vp_test.show(cps, corner_projecting, 'Surface + projecting points', at=1, viewup="z")#, lb2
#vp_test.show(cps, corner_projected, 'Surface + projected points', at=2, viewup="z")
#vp_test.show(cps, top_bot, 'Surface + projected curves', at=3, viewup="z")
#vp_test.show(cps, top_bot, side, 'Surface + projected curves + interpolated curves', at=4, viewup="z")
vp_test = Plotter(N=2, axes=1)
vp_test.show(cps, surface_points, TFI, 'Surface + transfinite interpolated points', at=0, viewup="z")
vp_test.show(cps, surface_points, bspline_fitted_cps, 'Surface + control points of fitted b-spline surface', at=1, viewup="z", interactive=True)
tem = Plotter()
tem.show(interactive=False)

# fig = plt.figure()
# ax1 = fig.gca(projection='3d')
# for surf in geo.input_bspline_entity_dict.values():
#     if surf.name == 'Surf_WFWKRQIMCA, Wing, 0, 12':
#         ax1.scatter(surf.control_points[:,0],surf.control_points[:,1],surf.control_points[:,2],marker = 'x')
#     if surf.name == 'Surf_WFWKRQIMCA, Wing, 0, 13':
#         ax1.scatter(surf.control_points[:,0],surf.control_points[:,1],surf.control_points[:,2],marker = 'x')
#     if surf.name == 'Surf_WFWKRQIMCA, Wing, 0, 14':
#         ax1.scatter(surf.control_points[:,0],surf.control_points[:,1],surf.control_points[:,2],marker = 'x')
#     if surf.name == 'Surf_WFWKRQIMCA, Wing, 0, 15':
#         ax1.scatter(surf.control_points[:,0],surf.control_points[:,1],surf.control_points[:,2],marker = 'x')
# plt.show()