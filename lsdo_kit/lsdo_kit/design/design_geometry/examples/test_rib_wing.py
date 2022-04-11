import numpy as np
import pyiges
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 

from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
import scipy.sparse as sps

import matplotlib.pyplot as plt

from vedo import Points, Plotter, LegendBox

# for debugging in vscode only
import os
os.chdir("../lsdo_geo/lsdo_kit/design_geometry/examples")

''' Rib creation script - demonstrates how to import a step file and process of creating a component rib'''

# importing geometry and creating geometry object
path_name = 'CAD/'
file_name = 'eVTOL.stp'
geo = DesignGeometry(path_name + file_name)

# specify what part of the 3D object the user wants to work with
wing_surface_names = [
    'Wing, 0, 12',
    'Wing, 0, 13',
    ]
top_wing_surface_names = [
    'Wing, 0, 13',
    ]
bot_wing_surface_names = [
    'Wing, 0, 12',
    ]
down_direction = np.array([0., 0., -1.])
up_direction = np.array([0., 0., 1.])

# Project top and bottom curve onto the surface, user will specify number of points they want to project onto the surface
# and specify the endpoints of the set of points - output is a pointset object that represents the top curve of the rib and physical coordinates of that curve

num_pts1 = 5

rib_top_curve, test_top_curve = geo.project_points(np.linspace(np.array([180.,63.,145.]),np.array([210.,63.,145.]), num_pts1),
    projection_targets_names=top_wing_surface_names, offset=np.array([0., 0., 0.5]), plot=True, projection_direction = down_direction)#top_wing_surface_names, 

rib_bot_curve, test_bot_curve = geo.project_points(np.linspace(np.array([180.,63.,105.]),np.array([210.,63.,105.]), num_pts1),
    projection_targets_names = bot_wing_surface_names, offset=np.array([0., 0., -0.5]), plot=True, projection_direction = up_direction)#

# Get the first and last indices of the projected point set - output is a pointset object that represents what indices the user wanted

top_projected_point1 = geo.extract_pointset(rib_top_curve, np.array([0]), np.array([1]))
top_projected_point2 = geo.extract_pointset(rib_top_curve, np.array([-1]), np.array([1]))
bot_projected_point1 = geo.extract_pointset(rib_bot_curve, np.array([0]), np.array([1]))
bot_projected_point2 = geo.extract_pointset(rib_bot_curve, np.array([-1]), np.array([1]))

# Create rib side curves through linear interpolation - output is a pointset object of linear interpolated points that represent the side curves

num_pts2 = [15]
rib_side_curve1 = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, num_pts2)
rib_side_curve2 = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, num_pts2)

# Create rib surface mesh through transfinite interpolation - output is a pointset object that represents the rib surface

rib_surface_mesh = geo.perform_2d_transfinite_interpolation(rib_top_curve, rib_bot_curve,  rib_side_curve1, rib_side_curve2)


geo.register_output(rib_surface_mesh, name = 'Rib')
geo.assemble()                                                  # Concatenates vertically all the linear matrices 
geo.evaluate()
points_to_be_fitted = geo.evaluate()                            # Evaluate the physical coornadites of points to be fitted
geo.fit_bspline_entities([rib_surface_mesh])
path_name = 'CAD/'
file_name = 'eVTOL_rib.igs'
geo.write_iges(path_name + file_name)

# plotting

# fig = plt.figure(figsize=(4,4))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(geo.total_cntrl_pts_vector[:,0],
# geo.total_cntrl_pts_vector[:,1],
# geo.total_cntrl_pts_vector[:,2]) # plot the point (2,3,4) on the figure
# plt.show()

corner_projecting = []
corner_projected = []
vp_test = Plotter(N=7, axes=1)
top_bot = []
side = []
cps = []
surface_points = []
corner_projecting.append(Points(np.vstack(([180.,60.,145.],[210.,60.,145.],[180.,60.,100.],[210.,60.,100.])), r=15, c='crimson').legend('Corner projecting'))
corner_projected.append(Points(np.vstack((test_top_curve[0,:],test_top_curve[-1,:],test_bot_curve[0,:],test_bot_curve[-1,:])), r=15, c='darkgoldenrod').legend('Corner projected'))
top_bot.append(Points(test_top_curve, r=15, c='brown').legend('Top curve'))
top_bot.append(Points(test_bot_curve, r=15, c='brown').legend('Bot curve'))
side.append(Points(np.linspace(test_top_curve[0,:],test_bot_curve[0,:],num_pts2[0]), r=15, c='chartreuse').legend('Side curve'))
side.append(Points(np.linspace(test_top_curve[-1,:],test_bot_curve[-1,:],num_pts2[0]), r=15, c='chartreuse').legend('Side curve'))
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
    if surf.name == 'Rib':
        bspline_fitted_cps = Points(surf.control_points, r=9, c='plum').legend('Fitted bspline')
vp_test.show(cps, surface_points, 'Control points of surface to be projected', at=0, viewup="z")#, lb1
vp_test.show(cps, surface_points, corner_projecting, 'Surface + projecting points', at=1, viewup="z")#, lb2
vp_test.show(cps, surface_points, corner_projected, 'Surface + projected points', at=2, viewup="z")
vp_test.show(cps, surface_points, top_bot, 'Surface + projected curves', at=3, viewup="z")
vp_test.show(cps, surface_points,  top_bot, side, 'Surface + projected curves + interpolated curves', at=4, viewup="z")
vp_test.show(cps, surface_points, top_bot, side, TFI, 'Surface + transfinite interpolated points', at=5, viewup="z")
vp_test.show(cps, surface_points, bspline_fitted_cps, 'Surface + control points of fitted b-spline surface', at=6, viewup="z", interactive=True)
       


