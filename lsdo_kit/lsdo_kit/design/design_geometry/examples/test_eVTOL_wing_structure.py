import numpy as np
import pyiges

from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
from lsdo_kit.old_files.mesh import Mesh

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
structures_mesh = Mesh(name='structures_mesh')

wing_surface_names = [
    'Wing, 0, 12',
    'Wing, 0, 13',
    'Wing, 0, 14',
    'Wing, 0, 15',
    ]
top_wing_surface_names = [ 
    'Wing, 0, 13', 
    'Wing, 0, 15', 
    ]
bot_wing_surface_names = [
    'Wing, 0, 12',
    'Wing, 0, 14', 
    ]
down_direction = np.array([0., 0., -1.])
up_direction = np.array([0., 0., 1.])
num_ribs = 5

# Project the points onto the surface
top_leading_root = np.array([185.,40.,160.])
top_trailing_root = np.array([215.,40.,160.])
top_leading_tip = np.array([153.5,190.,160.])
top_trailing_tip = np.array([171.5,190.,160.])

bot_leading_root = np.array([185.,40.,110.])
bot_trailing_root = np.array([215.,40.,110.])
bot_leading_tip = np.array([153.5,190.,110.])
bot_trailing_tip = np.array([171.5,190.,110.])

eps1 = 5e-1
# Construct upper and lower curves of primary spar
primary_spar_top = np.linspace(top_leading_root + np.array([0.,-20.,0.]), top_leading_tip+np.array([0.,24.,0.]), num_ribs)
primary_spar_top_curve, primary_spar_top_cood = geo.project_points(primary_spar_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), projection_direction = down_direction,plot = True)
# temp = Plotter()
# temp.show(interactive=True)
# exit()
primary_spar_bot = np.linspace(bot_leading_root + np.array([0.,-20.,0.]), bot_leading_tip+np.array([0.,24.,0.]), num_ribs)
primary_spar_bot_curve, primary_spar_bot_curve_cood = geo.project_points(primary_spar_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)
#print('test11')
# Construct upper and lower curves of rear spar
rear_spar_top = np.linspace(top_trailing_root + np.array([0.,-20.,0.]), top_trailing_tip + np.array([0.,22.,0.]), num_ribs)
rear_spar_top_curve, rear_spar_top_curve_cood = geo.project_points(rear_spar_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), projection_direction = down_direction)
rear_spar_bot = np.linspace(bot_trailing_root + np.array([0.,-20.,0.]), bot_trailing_tip + np.array([0.,22.,0.]), num_ribs)
rear_spar_bot_curve, rear_spar_bot_curve_cood = geo.project_points(rear_spar_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)
#print('test22')
num_side_pts_spars_and_ribs = 5
# Construct side curves of primary spar
top_projected_point1 = geo.extract_pointset(primary_spar_top_curve, np.array([0]), np.array([1]))
top_projected_point2 = geo.extract_pointset(primary_spar_top_curve, np.array([-1]), np.array([1]))
bot_projected_point1 = geo.extract_pointset(primary_spar_bot_curve, np.array([0]), np.array([1]))
bot_projected_point2 = geo.extract_pointset(primary_spar_bot_curve, np.array([-1]), np.array([1]))
primary_spar_side_root_curve = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, np.array([num_side_pts_spars_and_ribs]))
primary_spar_side_tip_curve = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, np.array([num_side_pts_spars_and_ribs]))
# Construct primary spar surface
primary_spar_surface_mesh = geo.perform_2d_transfinite_interpolation(primary_spar_top_curve, primary_spar_bot_curve, 
    primary_spar_side_root_curve, primary_spar_side_tip_curve)
# geo.register_output(primary_spar_surface_mesh, name = 'primary_spar')
structures_mesh.add_pointset(primary_spar_surface_mesh, name='primary_spar')

# Construct side curves of rear spar
top_projected_point1 = geo.extract_pointset(rear_spar_top_curve, np.array([0]), np.array([1]))
top_projected_point2 = geo.extract_pointset(rear_spar_top_curve, np.array([-1]), np.array([1]))
bot_projected_point1 = geo.extract_pointset(rear_spar_bot_curve, np.array([0]), np.array([1]))
bot_projected_point2 = geo.extract_pointset(rear_spar_bot_curve, np.array([-1]), np.array([1]))
rear_spar_side_root_curve = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, np.array([num_side_pts_spars_and_ribs]))
rear_spar_side_tip_curve = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, np.array([num_side_pts_spars_and_ribs]))
# Construct rear spar surface
rear_spar_surface_mesh = geo.perform_2d_transfinite_interpolation(rear_spar_top_curve, rear_spar_bot_curve, 
    rear_spar_side_root_curve, rear_spar_side_tip_curve)
# geo.register_output(rear_spar_surface_mesh, name = 'rear_spar')
structures_mesh.add_pointset(rear_spar_surface_mesh, name = 'rear_spar')


num_pts_ribs_top_bot = 7
eps2 = 7
eps3 = 5e-1
for i in range(num_ribs):
    # Construct upper and lower curves of ribs
    if i ==2:
        primary_spar_top = np.linspace(top_leading_root, top_leading_tip, num_ribs)-eps2 - 2
        rear_spar_top = np.linspace(top_trailing_root, top_trailing_tip, num_ribs)+eps3
        rear_spar_top[:,1] = primary_spar_top[:,1]
        rib_top = np.linspace(primary_spar_top[i] + np.array([0.,1.,0.]), rear_spar_top[i] + np.array([1.5,1.,0.]), num_pts_ribs_top_bot)
        rib_top_curve, rib_top_curve_cood = geo.project_points(rib_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), plot=True, projection_direction = down_direction)
        primary_spar_bot = np.linspace(bot_leading_root, bot_leading_tip, num_ribs)-eps2 - 2
        rear_spar_bot = np.linspace(bot_trailing_root, bot_trailing_tip, num_ribs)+eps3
        rear_spar_bot[:,1] = primary_spar_bot[:,1] 
        rib_bot = np.linspace(primary_spar_bot[i] + np.array([0.,1.,0.]), rear_spar_bot[i] + np.array([1.5,1.,0.]), num_pts_ribs_top_bot)
        rib_bot_curve, rib_bot_curve_cood = geo.project_points(rib_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), plot=True, projection_direction = up_direction)
        #print('test44')
    elif i==3:
        primary_spar_top = np.linspace(top_leading_root, top_leading_tip, num_ribs)-eps2 - 1
        rear_spar_top = np.linspace(top_trailing_root, top_trailing_tip, num_ribs)+eps3
        rear_spar_top[:,1] = primary_spar_top[:,1] 
        rib_top = np.linspace(primary_spar_top[i] + np.array([1.,0.,0.]), rear_spar_top[i] + np.array([4,0.,0.]), num_pts_ribs_top_bot)
        rib_top_curve, rib_top_curve_cood = geo.project_points(rib_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), plot=True, projection_direction = down_direction)
        
        primary_spar_bot = np.linspace(bot_leading_root, bot_leading_tip, num_ribs)-eps2 - 1
        rear_spar_bot = np.linspace(bot_trailing_root, bot_trailing_tip, num_ribs)+eps3
        rear_spar_bot[:,1] = primary_spar_bot[:,1]         
        rib_bot = np.linspace(primary_spar_bot[i] + np.array([1.,0.,0.]), rear_spar_bot[i] + np.array([4,0.,0.]), num_pts_ribs_top_bot)
        rib_bot_curve, rib_bot_curve_cood = geo.project_points(rib_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)
        #print('test55')
    elif i==4:
        primary_spar_top = np.linspace(top_leading_root, top_leading_tip, num_ribs)-eps2 - 1
        rear_spar_top = np.linspace(top_trailing_root, top_trailing_tip, num_ribs)+eps3
        rear_spar_top[:,1] = primary_spar_top[:,1] 
        rib_top = np.linspace(primary_spar_top[i] + np.array([3.,0.,0.]), rear_spar_top[i] + np.array([6.5,0.,0.]), num_pts_ribs_top_bot)
        rib_top_curve, rib_top_curve_cood = geo.project_points(rib_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), plot=True, projection_direction = down_direction)
        
        primary_spar_bot = np.linspace(bot_leading_root, bot_leading_tip, num_ribs)-eps2 - 1
        rear_spar_bot = np.linspace(bot_trailing_root, bot_trailing_tip, num_ribs)+eps3
        rear_spar_bot[:,1] = primary_spar_bot[:,1]         
        rib_bot = np.linspace(primary_spar_bot[i] + np.array([3.,0.,0.]), rear_spar_bot[i] + np.array([6.5,0.,0.]), num_pts_ribs_top_bot)
        rib_bot_curve, rib_bot_curve_cood = geo.project_points(rib_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)
        #print('test66')
    else: 
        primary_spar_top = np.linspace(top_leading_root, top_leading_tip, num_ribs)-eps2 - 2.5
        rear_spar_top = np.linspace(top_trailing_root, top_trailing_tip, num_ribs)+eps3
        rear_spar_top[:,1] = primary_spar_top[:,1]       
        rib_top = np.linspace(primary_spar_top[i], rear_spar_top[i] - np.array([eps3,0.,0.]), num_pts_ribs_top_bot)
        rib_top_curve, rib_top_curve_cood = geo.project_points(rib_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), plot=True, projection_direction = down_direction)

        primary_spar_bot = np.linspace(bot_leading_root, bot_leading_tip, num_ribs)-eps2 - 2.5
        rear_spar_bot = np.linspace(bot_trailing_root, bot_trailing_tip, num_ribs)+eps3
        rear_spar_bot[:,1] = primary_spar_bot[:,1]          
        rib_bot = np.linspace(primary_spar_bot[i], rear_spar_bot[i] - np.array([eps3,0.,0.]), num_pts_ribs_top_bot)
        rib_bot_curve, rib_bot_curve_cood = geo.project_points(rib_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)
        #print('test33')
    # Construct side curves of ribs
    num_side_pts_spars_and_ribs = 10
    top_leading_root_pt = geo.extract_pointset(rib_top_curve, np.array([0]), np.array([1]))
    bot_leading_root_pt = geo.extract_pointset(rib_bot_curve, np.array([0]), np.array([1]))
    rib_side_root_curve = geo.perform_linear_interpolation(top_leading_root_pt, bot_leading_root_pt, np.array([num_side_pts_spars_and_ribs]))
    top_leading_tip_pt = geo.extract_pointset(rib_top_curve, np.array([-1]), np.array([1]))
    bot_leading_tip_pt = geo.extract_pointset(rib_bot_curve, np.array([-1]), np.array([1]))
    rib_side_tip_curve = geo.perform_linear_interpolation(top_leading_tip_pt, bot_leading_tip_pt, np.array([num_side_pts_spars_and_ribs]))
    # Construct rib surface
    rib_surface_mesh = geo.perform_2d_transfinite_interpolation(rib_top_curve, rib_bot_curve, 
        rib_side_root_curve, rib_side_tip_curve)
    # geo.register_output(rib_surface_mesh, name = f'rib{i}')
    structures_mesh.add_pointset(rib_surface_mesh, name = f'rib{i}')
    eps2 -= 0.5
    
# Construct upper and lower curves of middle spar
primary_spar_top = np.linspace(top_leading_root + np.array([0.,-20.,0.]), top_leading_tip + np.array([0.,23.,0.]), num_ribs)
middle_spar_top = primary_spar_top
middle_spar_top[:,0] = 0.5 * (primary_spar_top[:,0] + rear_spar_top[:,0])
middle_spar_top_curve, middle_spar_top_curve_cood = geo.project_points(middle_spar_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), projection_direction = down_direction)
primary_spar_bot = np.linspace(bot_leading_root + np.array([0.,-20.,0.]), bot_leading_tip + np.array([0.,23.,0.]), num_ribs)
middle_spar_bot = primary_spar_bot
middle_spar_bot[:,0] = 0.5 * (primary_spar_bot[:,0] + rear_spar_bot[:,0])
middle_spar_bot_curve, middle_spar_bot_curve_cood = geo.project_points(middle_spar_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)
# Construct side curves of middle spar
top_projected_point1 = geo.extract_pointset(middle_spar_top_curve, np.array([0]), np.array([1]))
top_projected_point2 = geo.extract_pointset(middle_spar_top_curve, np.array([-1]), np.array([1]))
bot_projected_point1 = geo.extract_pointset(middle_spar_bot_curve, np.array([0]), np.array([1]))
bot_projected_point2 = geo.extract_pointset(middle_spar_bot_curve, np.array([-1]), np.array([1]))
middle_spar_side_root_curve = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, np.array([num_side_pts_spars_and_ribs]))
middle_spar_side_tip_curve = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, np.array([num_side_pts_spars_and_ribs]))
# Construct middle spar surface
middle_spar_surface_mesh = geo.perform_2d_transfinite_interpolation(middle_spar_top_curve, middle_spar_bot_curve, 
    middle_spar_side_root_curve, middle_spar_side_tip_curve)

# geo.register_output(middle_spar_surface_mesh, name = 'middle_spar')
structures_mesh.add_pointset(middle_spar_surface_mesh, name = 'middle_spar')



# Concatenates vertically all the linear matrices 
geo.assemble()
# Evaluate the physical coornadites of points to be fitted
points_to_be_fitted = geo.evaluate()
#print(points_to_be_fitted)
#print(points_to_be_fitted.shape)
time_end = time.time()
runtime = time_end - time_start
print('runtime', runtime)

geo.fit_bspline_entities(structures_mesh.pointset_list)

path_name = 'CAD/'#lsdo_geo/DesignGeometry/
file_name = 'eVTOL_wing_structure.igs'
geo.write_iges(path_name + file_name, plot = True)

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
for target in wing_surface_names:
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
for surf in geo.output_bspline_entity_dict.values():
    if surf.name == 'primary_spar':
        bspline_fitted_cps.append(Points(surf.control_points, r=12, c='black').legend('Fitted bspline'))
    if  surf.name == 'rear_spar':
        bspline_fitted_cps.append(Points(surf.control_points, r=12, c='yellow').legend('Fitted bspline'))
    if  'rib' in surf.name:
        bspline_fitted_cps.append(Points(surf.control_points, r=9, c='red').legend('Fitted bspline'))
    if  surf.name == 'middle_spar':
        bspline_fitted_cps.append(Points(surf.control_points, r=12, c='blue').legend('Fitted bspline'))
#vp_test.show(cps, 'Control points of surface to be projected', at=0, viewup="z")#, lb1
#vp_test.show(cps, corner_projecting, 'Surface + projecting points', at=1, viewup="z")#, lb2
#vp_test.show(cps, corner_projected, 'Surface + projected points', at=2, viewup="z")
#vp_test.show(cps, top_bot, 'Surface + projected curves', at=3, viewup="z")
#vp_test.show(cps, top_bot, side, 'Surface + projected curves + interpolated curves', at=4, viewup="z")
vp_test = Plotter(N=2, axes=1)
vp_test.show(cps, surface_points, TFI, 'Surface + transfinite interpolated points', at=0, viewup="z")
vp_test.show(cps, surface_points, bspline_fitted_cps, 'Surface + control points of fitted b-spline surface', at=1, viewup="z", interactive=True)
temp = Plotter()
temp.show(interactive=True)

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