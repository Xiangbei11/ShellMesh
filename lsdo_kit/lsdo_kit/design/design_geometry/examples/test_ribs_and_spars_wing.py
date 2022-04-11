import numpy as np
import pyiges

from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 

import matplotlib.pyplot as plt
import time

from vedo import Points, Plotter, LegendBox

# from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
import scipy.sparse as sps

time_start = time.time()

''' Spars and ribs creation script '''

path_name = 'CAD/'
file_name = 'eVTOL.stp'
geo = DesignGeometry(path_name + file_name, plot=True)

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
top_leading_root = np.array([185.,40.,150.])
top_trailing_root = np.array([210.,40.,140.])
top_leading_tip = np.array([158.,190.,150.])
top_trailing_tip = np.array([180.,190.,145.])

bot_leading_root = np.array([185.,40.,105.])
bot_trailing_root = np.array([210.,40.,105.])
bot_leading_tip = np.array([158.,190.,105.])
bot_trailing_tip = np.array([180.,190.,105.])

eps1 = 1e-1
# Construct upper and lower curves of primary spar
primary_spar_top = np.linspace(top_leading_root+np.array([0.,-1.,0.]), top_leading_tip, num_ribs)
primary_spar_top_curve, test_bot_pri_spar = geo.project_points(primary_spar_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), plot=True, projection_direction = down_direction)
primary_spar_top = np.linspace(top_leading_root, top_leading_tip, num_ribs)

primary_spar_bot = np.linspace(bot_leading_root+np.array([0.,-1.,0.]), bot_leading_tip, num_ribs)
primary_spar_bot_curve, test_bot_pri_spar = geo.project_points(primary_spar_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), plot=True, projection_direction = up_direction)
primary_spar_bot = np.linspace(bot_leading_root, bot_leading_tip, num_ribs)

# Construct upper and lower curves of rear spar
rear_spar_top = np.linspace(top_trailing_root, top_trailing_tip+np.array([0.,1.,0.]), num_ribs)
rear_spar_top_curve, test_bot_pri_spar = geo.project_points(rear_spar_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), plot=True, projection_direction = down_direction)
rear_spar_top = np.linspace(top_trailing_root, top_trailing_tip, num_ribs)

rear_spar_bot = np.linspace(bot_trailing_root, bot_trailing_tip+np.array([0.,1.,0.]), num_ribs)
rear_spar_bot_curve, test_bot_pri_spar = geo.project_points(rear_spar_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), plot=True, projection_direction = up_direction)
rear_spar_bot = np.linspace(bot_trailing_root, bot_trailing_tip, num_ribs)

# Construct side curves of primary spar
num_side_pts_spars_and_ribs = 5

top_projected_point1 = geo.extract_pointset(primary_spar_top_curve, np.array([0]), np.array([1]))
top_projected_point2 = geo.extract_pointset(primary_spar_top_curve, np.array([-1]), np.array([1]))
bot_projected_point1 = geo.extract_pointset(primary_spar_bot_curve, np.array([0]), np.array([1]))
bot_projected_point2 = geo.extract_pointset(primary_spar_bot_curve, np.array([-1]), np.array([1]))
primary_spar_side_root_curve = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, np.array([num_side_pts_spars_and_ribs]))
primary_spar_side_tip_curve = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, np.array([num_side_pts_spars_and_ribs]))
# Construct primary spar surface
primary_spar_surface_mesh = geo.perform_2d_transfinite_interpolation(primary_spar_top_curve, primary_spar_bot_curve, 
    primary_spar_side_root_curve, primary_spar_side_tip_curve)
geo.register_output(primary_spar_surface_mesh, name = 'primary_spar')


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
geo.register_output(rear_spar_surface_mesh, name = 'rear_spar')

num_pts_ribs_top_bot = 7
eps2 = 1
for i in range(num_ribs):
    # Construct upper and lower curves of ribs
    rib_top = np.linspace(primary_spar_top[i]-eps2, rear_spar_top[i]+eps2, num_pts_ribs_top_bot)
    rib_top_curve, test_top_rib = geo.project_points(rib_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), plot=True, projection_direction = down_direction)
    
    rib_bot = np.linspace(primary_spar_bot[i]-eps2, rear_spar_bot[i]+eps2, num_pts_ribs_top_bot)
    rib_bot_curve, test_bot_rib = geo.project_points(rib_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), plot=True, projection_direction = up_direction)


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
    geo.register_output(rib_surface_mesh, name = f'rib{i}')

# Concatenates vertically all the linear matrices 
geo.assemble()
# Evaluate the physical coornadites of points to be fitted
points_to_be_fitted = geo.evaluate()
#print(points_to_be_fitted)
#print(points_to_be_fitted.shape)
time_end = time.time()
runtime = time_end - time_start
print('runtime', runtime)

# Fit Bspline Curve/Surface with least square method (Volume//TODO)
# For now,  after the user call fit_bspline_entities, 
# it will add a new bspline entity object to the output_bspline_entity_dict 
# wihch also has all the input bspline entity object
geo.fit_bspline_entities(points_to_be_fitted)

path_name = 'CAD/'
file_name = 'eVTOL_ribs_and_spars.igs'
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
    cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))
    bspline_entity = geo.input_bspline_entity_dict[target]

    num_points_u = 50   # TODO might want to pass these in as input
    num_points_v = 50

    u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u), np.ones(num_points_v)).flatten()
    v_vec = np.einsum('i,j->ij', np.ones(num_points_u), np.linspace(0., 1., num_points_v)).flatten()

    pts = bspline_entity.evaluate_points(u_vec, v_vec)

    surface_points.append(Points(pts, r=5, c='seagreen',alpha=0.3).legend('Surface points'))
    cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))          
for surf in geo.output_bspline_entity_dict.values():
    if surf.name == 'primary_spar':
        bspline_fitted_cps.append(Points(surf.control_points, r=12, c='black', alpha=0.5).legend('Fitted bspline'))
    if  surf.name == 'rear_spar':
        bspline_fitted_cps.append(Points(surf.control_points, r=12, c='yellow', alpha=0.5).legend('Fitted bspline'))
    if  'rib' in surf.name:
        bspline_fitted_cps.append(Points(surf.control_points, r=9, c='red').legend('Fitted bspline'))
#vp_test.show(cps, 'Control points of surface to be projected', at=0, viewup="z")#, lb1
#vp_test.show(cps, corner_projecting, 'Surface + projecting points', at=1, viewup="z")#, lb2
#vp_test.show(cps, corner_projected, 'Surface + projected points', at=2, viewup="z")
#vp_test.show(cps, top_bot, 'Surface + projected curves', at=3, viewup="z")
#vp_test.show(cps, top_bot, side, 'Surface + projected curves + interpolated curves', at=4, viewup="z")
vp_test = Plotter(N=2, axes=1)
vp_test.show(cps, surface_points, TFI, 'Surface + transfinite interpolated points', at=0, viewup="z")
vp_test.show(cps, surface_points, bspline_fitted_cps, 'Surface + control points of fitted b-spline surface', at=1, viewup="z", interactive=True)
exit()       

