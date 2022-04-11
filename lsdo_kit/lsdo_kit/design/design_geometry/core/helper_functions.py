'''
Temporary helper functions file
'''

import numpy as np

def create_rib_or_spar(geo, projection_target_names, point00, point10, point01, point11, shape, epsilon):
    down_direction = np.array([0., 0., -1.])
    up_direction = np.array([0., 0., 1.])

    # Project to find top and bottom curves
    top_curve_points_to_be_projected = np.linspace(point01, point11, shape[0])
    top_curve, top_curve_eval_pts = geo.project_curve(top_curve_points_to_be_projected, 
        projection_targets_names=projection_target_names, offset=np.array([0., 0., epsilon]), projection_direction = down_direction, plot=False)

    bot_curve_points_to_be_projected = np.linspace(point00, point10, shape[0])
    bot_curve, bot_curve_eval_pts = geo.project_curve(bot_curve_points_to_be_projected,
        projection_targets_names=projection_target_names, offset=np.array([0., 0., -epsilon]), projection_direction = up_direction, plot=False)

    # Construct side curves of primary spar
    top_left_corner = geo.extract_pointset(top_curve, np.array([0]), np.array([1]))
    top_right_corner = geo.extract_pointset(top_curve, np.array([-1]), np.array([1]))
    bot_left_corner = geo.extract_pointset(bot_curve, np.array([0]), np.array([1]))
    bot_right_corner = geo.extract_pointset(bot_curve, np.array([-1]), np.array([1]))
    left_curve = geo.perform_linear_interpolation(top_left_corner, bot_left_corner, np.array([shape[1]]))
    right_curve = geo.perform_linear_interpolation(top_right_corner, bot_right_corner, np.array([shape[1]]))
    # Construct primary spar surface
    surf = geo.perform_2d_transfinite_interpolation(top_curve, bot_curve, left_curve, right_curve)
    return surf

def create_uniform_ribs_and_spars(geo, projection_target_names, point000, point010, point001, point011,
        point100, point110, point101, point111, overall_structure_shape, member_structure_shape, epsilon):
    
    # TODO allow the user to generate a mesh and return an array of these that can be registered by the user
    lower_spar_root_locations = np.linspace(point000, point100, overall_structure_shape[0])
    lower_spar_tip_locations = np.linspace(point010, point110, overall_structure_shape[0])
    upper_spar_root_locations = np.linspace(point001, point101, overall_structure_shape[0])
    upper_spar_tip_locations = np.linspace(point011, point111, overall_structure_shape[0])
    spars = []
    for i in range(overall_structure_shape[0]):
        print('spar', i)
        spars.append(create_rib_or_spar(geo, projection_target_names, lower_spar_root_locations[i], lower_spar_tip_locations[i],
            upper_spar_root_locations[i], upper_spar_tip_locations[i], member_structure_shape, epsilon))

    lower_rib_leading_edge_locations = np.linspace(point000, point010, overall_structure_shape[1])
    lower_rib_trailing_edge_locations = np.linspace(point100, point110, overall_structure_shape[1])
    upper_rib_leading_edge_locations = np.linspace(point001, point011, overall_structure_shape[1])
    upper_rib_trailing_edge_locations = np.linspace(point101, point111, overall_structure_shape[1])
    ribs = []
    for i in range(overall_structure_shape[1]):
        print('rib', i)
        ribs.append(create_rib_or_spar(geo, projection_target_names, lower_rib_leading_edge_locations[i], lower_rib_trailing_edge_locations[i],
            upper_rib_leading_edge_locations[i], upper_rib_trailing_edge_locations[i], member_structure_shape, epsilon))

    return spars, ribs
