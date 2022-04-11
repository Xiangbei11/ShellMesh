import numpy as np
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 

import matplotlib.pyplot as plt


#import matplotlib.pyplot as plt

''' Rib creation script with test wing '''

path_name = 'CAD/'
file_name = 'test_wing.stp'
geo = DesignGeometry(path_name + file_name)



# fig = plt.figure()
# ax1 = fig.gca(projection='3d')
# #ax1.scatter(geo.total_cntrl_pts_vector[:,0], geo.total_cntrl_pts_vector[:,1], geo.total_cntrl_pts_vector[:,2], s= 10, marker = 'o', alpha=0.3)
# for surf in geo.input_bspline_entity_dict.values():
#     ax1.scatter(surf.control_points[:,0],surf.control_points[:,1],surf.control_points[:,2], marker = 'o',s = 10, alpha =0.1)
# ax1.scatter(0.04,0.5,0.03)
# ax1.scatter(0.2 ,0.5,0.03)
# ax1.scatter(0.04,0.5,-0.02)
# ax1.scatter(0.2 ,0.5,-0.02)
top_curve_pts = np.linspace(np.array([0.04,0.5,0.03]),np.array([0.2 ,0.5,0.03]), 10)
bot_curve_pts = np.linspace(np.array([0.04,0.5,-0.02]),np.array([0.2 ,0.5,-0.02]), 10)
# ax1.scatter(top_curve_pts[:,0], top_curve_pts[:,1], top_curve_pts[:,2])
# ax1.scatter(bot_curve_pts[:,0], bot_curve_pts[:,1], bot_curve_pts[:,2])
# plt.show()

# Project top and bottom curve onto the surface
num_pts1 = 10
direction = np.array([0., 0., -1.])
rib_top_curve, pts_top_curve = geo.project_points(np.linspace(np.array([0.04,0.5,0.03]),np.array([0.2 ,0.5,0.03]), num_pts1), projection_direction=direction)
direction = np.array([0., 0., 1.])
rib_bot_curve, pts_bot_curve = geo.project_points(np.linspace(np.array([0.04,0.5,-0.02]),np.array([0.2 ,0.5,-0.02]), num_pts1), projection_direction=direction)

# # Project the points onto the surface 
# # Note: Not necessary now that we have the extract method which is used below
# down_direction = np.array([0., 0., -1.])
# top_projected_point1,_ = geo.project_point([0.04,0.5,0.03], projection_direction=down_direction)
# top_projected_point2,_ = geo.project_point([0.2 ,0.5,0.03], projection_direction=down_direction)
# up_direction = np.array([0., 0., 1.])
# bot_projected_point1,_ = geo.project_point([0.04,0.5,-0.02], projection_direction=up_direction)
# bot_projected_point2,_ = geo.project_point([0.2 ,0.5,-0.02], projection_direction=up_direction)

top_projected_point1 = geo.extract_pointset(rib_top_curve, np.array([0]), np.array([1]))
top_projected_point2 = geo.extract_pointset(rib_top_curve, np.array([-1]), np.array([1]))
bot_projected_point1 = geo.extract_pointset(rib_bot_curve, np.array([0]), np.array([1]))
bot_projected_point2 = geo.extract_pointset(rib_bot_curve, np.array([-1]), np.array([1]))

# Create rib side curves
num_pts2 = 15
rib_side_curve1 = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, num_pts2)
rib_side_curve2 = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, num_pts2)

# Define rib surface mesh
rib_surface_mesh = geo.perform_2d_transfinite_interpolation(rib_top_curve, rib_bot_curve,  rib_side_curve1, rib_side_curve2)

# Register output
geo.register_output(rib_surface_mesh, name = 'Rib')
# Concatenates vertically all the linear matrices 
geo.assemble()
# Evaluate the physical coornadites of points to be fitted
points_to_be_fitted = geo.evaluate()

# Fit Bspline Curve/Surface with least square method (Volume//TODO)
# For now,  after the user call fit_bspline_entities, 
# it will add a new bspline entity object to the output_bspline_entity_dict 
geo.fit_bspline_entities(points_to_be_fitted)

path_name = 'lsdo_geo/DesignGeometry/CAD/'
file_name = 'test_wing_rib.igs'
geo.write_iges(path_name + file_name, plot = True)
