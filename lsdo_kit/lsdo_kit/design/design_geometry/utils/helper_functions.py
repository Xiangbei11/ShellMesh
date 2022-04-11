'''
Temporary helper functions file
'''

import numpy as np

def create_rib_or_spar(geo, projection_target_names, point00, point10, point01, point11, shape, epsilon):
    down_direction = np.array([0., 0., -1.])
    up_direction = np.array([0., 0., 1.])

    # Project to find top and bottom curves
    top_curve_points_to_be_projected = np.linspace(point01, point11, shape[0])
    top_curve, top_curve_cood = geo.project_points(top_curve_points_to_be_projected, 
        projection_targets_names=projection_target_names, offset=np.array([0., 0., epsilon]), projection_direction = down_direction, plot=False)

    bot_curve_points_to_be_projected = np.linspace(point00, point10, shape[0])
    bot_curve, bot_curve_cood = geo.project_points(bot_curve_points_to_be_projected,
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

# def create_rib_from_chord(geo, projection_target_names, point00, point01, point10, point11,shape,loc_front_back=[0.5,0,1]):
def create_element_from_chord(geo, top_target_names, bot_target_names, point00, point01, point10, point11,shape,loc_front_back=[0.5,0,1]):
    #shape = [num_points_along_chord, num_points_along_span,num_points_height)
    #rib_locations must be np.array with ndmin = 2
    #loc_front_back = [location,front,back]
    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])
    
    '''Project points'''
    wing_lead_root, wing_lead_root_coord = geo.project_points(point00, projection_targets_names = top_target_names, projection_direction = down_direction)
    wing_lead_tip, wing_lead_tip_coord = geo.project_points(point01, projection_targets_names = top_target_names, projection_direction = down_direction)
    wing_trail_tip, wing_trail_tip_coord = geo.project_points(point10, projection_targets_names = top_target_names, projection_direction = down_direction)
    wing_trail_root, wing_trail_root_coord = geo.project_points(point11, projection_targets_names = top_target_names, projection_direction = down_direction)

    
    # #Create Chord Lines First
    # if all(shape[:2] != 1):
    #     print('Dimensions of pointset make ambiguous geometry. shape should have 1 in one of the first two arguments.')
    #     return
    # elif shape[0] == 1:#Spar
    #     output_parameters0 = np.linspace(loc_front_back[1],loc_front_back[2],shape[1])
    #     root_curve = geo.perform_linear_interpolation(wing_lead_root,wing_trail_root,shape[0],output_parameters=loc_front_back[0])
    #     tip_curve = geo.perform_linear_interpolation(wing_lead_tip,wing_trail_tip,shape[0],output_parameters=loc_front_back[0])
    #     span_line = geo.perform_linear_interpolation(root_curve,tip_curve, shape[0:2], output_parameters=output_parameters0)
    # elif shape[1] == 1:#Rib
    #     output_parameters0 = np.linspace(loc_front_back[1],loc_front_back[2],shape[0])
    #     root_curve = geo.perform_linear_interpolation(wing_lead_root,wing_trail_root,shape[0],output_parameters=output_parameters0)
    #     tip_curve = geo.perform_linear_interpolation(wing_lead_tip,wing_trail_tip,shape[0],output_parameters=output_parameters0)
    #     span_line = geo.perform_linear_interpolation(root_curve,tip_curve, shape[0:2], output_parameters=loc_front_back[0])

    #Change Which Edges to Start With
    if all(shape[:2] != 1):
        print('Dimensions of pointset make ambiguous geometry. shape should have 1 in one of the first two arguments.')
        return
    elif shape[0] == 1:#Spar
        output_parameters0 = np.linspace(loc_front_back[1],loc_front_back[2],shape[1])
        side_curve1 = geo.perform_linear_interpolation(wing_lead_root,wing_lead_tip,[shape[1]],output_parameters=output_parameters0)
        side_curve2= geo.perform_linear_interpolation(wing_trail_root,wing_trail_tip,[shape[1]],output_parameters=output_parameters0)
        output_parameters1 = loc_front_back[0]*np.ones([shape[1],1])
        span_line = geo.perform_linear_interpolation(side_curve1,side_curve2,[shape[1],shape[0]],output_parameters=output_parameters1)#yuck
    elif shape[1] == 1:#Rib
        # print('test1')
        output_parameters0 = np.linspace(loc_front_back[1],loc_front_back[2],shape[0])
        side_curve1 = geo.perform_linear_interpolation(wing_lead_root,wing_trail_root,[shape[0]],output_parameters=output_parameters0)
        side_curve2 = geo.perform_linear_interpolation(wing_lead_tip,wing_trail_tip,[shape[0]],output_parameters=output_parameters0)
        output_parameters1 = loc_front_back[0]*np.ones([shape[0],1])    
        span_line = geo.perform_linear_interpolation(side_curve1,side_curve2,shape[:2],output_parameters=output_parameters1)
    
    geo.assemble(span_line)
    span_line_eval = geo.evaluate(span_line)

    '''Project onto Wing Surfaces'''
    # print('test2')
    span_shift_up = span_line_eval
    span_shift_up[2] = span_line_eval[2]+12

    span_shift_down = span_line_eval
    span_shift_down[2] = span_line_eval[2]-12

    # top_points, top_plot = geo.project_surface(span_line_eval, projection_target_names, up_direction)
    # bot_points, bot_plot = geo.project_surface(span_line_eval, projection_target_names, down_direction)
    top_points, top_plot = geo.project_points(span_shift_up, top_target_names, down_direction)
    bot_points, bot_plot = geo.project_points(span_shift_down, bot_target_names, up_direction)
    # top_points, top_plot = geo.project_surface(span_line, top_target_names, up_direction)
    # bot_points, bot_plot = geo.project_surface(span_line, bot_target_names, down_direction)

    '''Interpolate between wing surfaces'''
    geo.assemble(top_points)
    tp = geo.evaluate(top_points)
    # print(tp)
    # print('test4')
    geo.assemble(bot_points)
    bp = geo.evaluate(bot_points)
    # print(bp)
    geometry = geo.perform_linear_interpolation(top_points, bot_points, shape)
    # print('test5')
    return geometry
    

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
