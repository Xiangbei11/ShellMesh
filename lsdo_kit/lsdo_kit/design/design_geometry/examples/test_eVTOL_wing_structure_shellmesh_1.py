import numpy as np
import pyiges

from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
from lsdo_kit.old_files.mesh import Mesh

import matplotlib.pyplot as plt
import time

import vedo
from vedo import Points, Plotter, LegendBox

from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
import scipy.sparse as sps

time_start = time.time()

''' Spars and ribs creation script '''
path_name = 'CAD/'
file_name = 'eVTOL_wing.stp' #_wing_with_tip #_wing_structure_PENGoLINS
geo = DesignGeometry(path_name + file_name)

# path_name = 'CAD/'
# file_name = 'eVTOL_wing_with_tip.igs'
# geo.write_iges(path_name + file_name, plot = True)
# exit()

structures_mesh = Mesh(name='structures_mesh')
wing_surface_names = [
    'Surface, 0',
    'Surface, 1',
    'Surface, 2',
    'Surface, 3',
    ]
top_wing_surface_names = [ 
    'Surface, 1', 
    'Surface, 3', 
    ]
bot_wing_surface_names = [
    'Surface, 0',
    'Surface, 2', 
    ]
down_direction = np.array([0., 0., -1.])
up_direction = np.array([0., 0., 1.])

# Project the points onto the surface
x_leading_root = 190
x_leading_tip = 155
y_leading_root = 20
y_leading_tip = 210

x_trailing_root = 220
x_trailing_tip = 175
y_trailing_root = y_leading_root
y_trailing_tip = y_leading_tip

top_leading_root = np.array([x_leading_root,y_leading_root,160.]) 
top_trailing_root = np.array([x_trailing_root,y_trailing_root,160.])
top_leading_tip = np.array([x_leading_tip,y_leading_tip,160.]) 
top_trailing_tip = np.array([x_trailing_tip,y_trailing_tip,160.]) 

bot_leading_root = np.array([x_leading_root,y_leading_root,110.]) 
bot_trailing_root = np.array([x_trailing_root,y_trailing_root,110.]) 
bot_leading_tip = np.array([x_leading_tip,y_leading_tip,110.]) 
bot_trailing_tip = np.array([x_trailing_tip,y_trailing_tip,110.]) 

eps1 = 0
num_pts_spars_top_bot = 10
num_side_pts_spars = 7

# Construct upper and lower curves of primary spar
primary_spar_top = np.linspace(top_leading_root, top_leading_tip, num_pts_spars_top_bot)
primary_spar_top_curve, primary_spar_top_cood = geo.project_points(primary_spar_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), projection_direction = down_direction)
primary_spar_bot = np.linspace(bot_leading_root, bot_leading_tip, num_pts_spars_top_bot)
primary_spar_bot_curve, primary_spar_bot_curve_cood = geo.project_points(primary_spar_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)
# Construct side curves of primary spar
top_projected_point1 = geo.extract_pointset(primary_spar_top_curve, np.array([0]), np.array([1]))
top_projected_point2 = geo.extract_pointset(primary_spar_top_curve, np.array([-1]), np.array([1]))
bot_projected_point1 = geo.extract_pointset(primary_spar_bot_curve, np.array([0]), np.array([1]))
bot_projected_point2 = geo.extract_pointset(primary_spar_bot_curve, np.array([-1]), np.array([1]))
primary_spar_side_root_curve = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, np.array([num_side_pts_spars]))
primary_spar_side_tip_curve = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, np.array([num_side_pts_spars]))
# Construct primary spar surface
primary_spar_surface_mesh = geo.perform_2d_transfinite_interpolation(primary_spar_top_curve, primary_spar_bot_curve, 
    primary_spar_side_root_curve, primary_spar_side_tip_curve)
# geo.register_output(primary_spar_surface_mesh, name = 'primary_spar')
structures_mesh.add_pointset(primary_spar_surface_mesh, name='primary_spar')

# Construct upper and lower curves of rear spar
rear_spar_top = np.linspace(top_trailing_root, top_trailing_tip, num_pts_spars_top_bot)
rear_spar_top_curve, rear_spar_top_curve_cood = geo.project_points(rear_spar_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), projection_direction = down_direction)
rear_spar_bot = np.linspace(bot_trailing_root, bot_trailing_tip, num_pts_spars_top_bot)
rear_spar_bot_curve, rear_spar_bot_curve_cood = geo.project_points(rear_spar_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)
# Construct side curves of rear spar
top_projected_point1 = geo.extract_pointset(rear_spar_top_curve, np.array([0]), np.array([1]))
top_projected_point2 = geo.extract_pointset(rear_spar_top_curve, np.array([-1]), np.array([1]))
bot_projected_point1 = geo.extract_pointset(rear_spar_bot_curve, np.array([0]), np.array([1]))
bot_projected_point2 = geo.extract_pointset(rear_spar_bot_curve, np.array([-1]), np.array([1]))
rear_spar_side_root_curve = geo.perform_linear_interpolation(top_projected_point1, bot_projected_point1, np.array([num_side_pts_spars]))
rear_spar_side_tip_curve = geo.perform_linear_interpolation(top_projected_point2, bot_projected_point2, np.array([num_side_pts_spars]))
# Construct rear spar surface
rear_spar_surface_mesh = geo.perform_2d_transfinite_interpolation(rear_spar_top_curve, rear_spar_bot_curve, 
    rear_spar_side_root_curve, rear_spar_side_tip_curve)
# geo.register_output(rear_spar_surface_mesh, name = 'rear_spar')
structures_mesh.add_pointset(rear_spar_surface_mesh, name = 'rear_spar')


num_ribs = 12
num_pts_ribs_top_bot = 7
num_side_pts_ribs = 7
num_pts = 14
rib_top_leading_root = np.linspace(top_leading_root, top_leading_tip, num_pts)[1]
rib_top_leading_tip = np.linspace(top_leading_root, top_leading_tip, num_pts)[-2]
rib_top_trailing_root = np.linspace(top_trailing_root, top_trailing_tip, num_pts)[1]
rib_top_trailing_tip = np.linspace(top_trailing_root, top_trailing_tip, num_pts)[-2]
rib_bot_leading_root = np.linspace(bot_leading_root, bot_leading_tip, num_pts)[1]
rib_bot_leading_tip = np.linspace(bot_leading_root, bot_leading_tip, num_pts)[-2]
rib_bot_trailing_root = np.linspace(bot_trailing_root, bot_trailing_tip, num_pts)[1]
rib_bot_trailing_tip = np.linspace(bot_trailing_root, bot_trailing_tip, num_pts)[-2]

primary_spar_top = np.linspace(rib_top_leading_root, rib_top_leading_tip, num_ribs)
rear_spar_top = np.linspace(rib_top_trailing_root, rib_top_trailing_tip, num_ribs)
rear_spar_top[:,1] = primary_spar_top[:,1] 
primary_spar_bot = np.linspace(rib_bot_leading_root, rib_bot_leading_tip, num_ribs)
rear_spar_bot = np.linspace(rib_bot_trailing_root, rib_bot_trailing_tip, num_ribs)
rear_spar_bot[:,1] = primary_spar_bot[:,1]

for i in range(num_ribs):
    # Construct upper and lower curves of ribs      
    rib_top = np.linspace(primary_spar_top[i], rear_spar_top[i], num_pts_ribs_top_bot)
    rib_top_curve, rib_top_curve_cood = geo.project_points(rib_top, top_wing_surface_names, offset=np.array([0., 0., eps1]), projection_direction = down_direction)
         
    rib_bot = np.linspace(primary_spar_bot[i], rear_spar_bot[i], num_pts_ribs_top_bot)
    rib_bot_curve, rib_bot_curve_cood = geo.project_points(rib_bot, bot_wing_surface_names, offset=np.array([0., 0., -eps1]), projection_direction = up_direction)

    # Construct side curves of ribs
    top_leading_root_pt = geo.extract_pointset(rib_top_curve, np.array([0]), np.array([1]))
    bot_leading_root_pt = geo.extract_pointset(rib_bot_curve, np.array([0]), np.array([1]))
    rib_side_root_curve = geo.perform_linear_interpolation(top_leading_root_pt, bot_leading_root_pt, np.array([num_side_pts_ribs]))
    top_leading_tip_pt = geo.extract_pointset(rib_top_curve, np.array([-1]), np.array([1]))
    bot_leading_tip_pt = geo.extract_pointset(rib_bot_curve, np.array([-1]), np.array([1]))
    rib_side_tip_curve = geo.perform_linear_interpolation(top_leading_tip_pt, bot_leading_tip_pt, np.array([num_side_pts_ribs]))
    
    # Construct rib surface
    rib_surface_mesh = geo.perform_2d_transfinite_interpolation(rib_top_curve, rib_bot_curve, 
        rib_side_root_curve, rib_side_tip_curve)
    structures_mesh.add_pointset(rib_surface_mesh, name = f'rib{i}')
    
# Concatenates vertically all the linear matrices 
geo.assemble()
# Evaluate the physical coornadites of points to be fitted
points_to_be_fitted = geo.evaluate()
time_end = time.time()
runtime = time_end - time_start

members_ctrl_pointset_list = geo.fit_bspline_ctrl_pointsets(structures_mesh.pointset_list)

geo.fit_bspline_entities(structures_mesh.pointset_list)
print('runtime', runtime)

# path_name = 'CAD/'
# file_name = 'shellmesh_test_eVTOL_wing_structure_2.igs'
# geo.write_iges(path_name + file_name, plot = False)

'''plot'''
cpts = []
bspline_fitted_cps = []
spts = []

for target in wing_surface_names:
    bspline_entity = geo.input_bspline_entity_dict[target]
    order_u = bspline_entity.order_u
    order_v = bspline_entity.order_v
    num_control_points_u = bspline_entity.shape[0]
    num_control_points_v = bspline_entity.shape[1]
    num_points_u = 50   
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
    spts.append(Points(pts, r=5, c='seagreen',alpha=0.3).legend('Surface points'))
    cpts.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Whole Structure control points'))          

for surf in geo.output_bspline_entity_dict.values():
    if surf.name == 'primary_spar':
        bspline_fitted_cps.append(Points(surf.control_points, r=12, c='black').legend('Fitted bspline'))
    if  surf.name == 'rear_spar':
        bspline_fitted_cps.append(Points(surf.control_points, r=12, c='yellow').legend('Fitted bspline'))
    if  'rib' in surf.name:
        bspline_fitted_cps.append(Points(surf.control_points, r=9, c='red').legend('Fitted bspline'))

vp_test = Plotter(axes=1)
#vp_test.show(cpts, spts, bspline_fitted_cps, 'Wing + control points of fitted b-spline surface', viewup="z", interactive=False)
