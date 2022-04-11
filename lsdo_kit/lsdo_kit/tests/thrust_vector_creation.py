import numpy as np
from lsdo_kit.design.design_geometry.core.component import Component

def generate_thrust_vector(geo, nacelle_comp):


    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    _nacelle_length = nacelle_comp.x_max - nacelle_comp.x_min
    _quarter_nacelle_length = (1/4) * _nacelle_length
    _tip_indices = np.argwhere(nacelle_comp.embedded_entities_control_points == nacelle_comp.x_min)
    _tip_index = _tip_indices[:,0]

    _tip_centroid = np.array([0., 0., 0.])

    for _tip_index in _tip_indices[:,0]:
        thrust_tip_control_point = nacelle_comp.embedded_entities_control_points[_tip_index, :]
        _tip_centroid = _tip_centroid + thrust_tip_control_point

    _tip_centroid = (_tip_centroid / len(_tip_indices[:,0]) )
    # print('_tip_centroid: ', _tip_centroid)

    _nacelle_z_max = nacelle_comp.z_max
    _nacelle_z_min = nacelle_comp.z_min

    thrust_tip_control_point = _tip_centroid

    thrust_top_origin_control_point = np.array([0., 0., 0.])
    thrust_top_origin_control_point[0] = thrust_tip_control_point[0] + _quarter_nacelle_length
    thrust_top_origin_control_point[1] = thrust_tip_control_point[1] 
    thrust_top_origin_control_point[2] = _nacelle_z_max

    thrust_bot_origin_control_point = np.array([0., 0., 0.])
    thrust_bot_origin_control_point[0] = thrust_tip_control_point[0] + _quarter_nacelle_length
    thrust_bot_origin_control_point[1] = thrust_tip_control_point[1]
    thrust_bot_origin_control_point[2] = _nacelle_z_min

    thrust_top_origin, thrust_top_origin_coord = geo.project_points(thrust_top_origin_control_point, projection_direction = down_direction, projection_targets_names=[nacelle_comp.name])
    thrust_bot_origin, thrust_bot_origin_coord = geo.project_points(thrust_bot_origin_control_point, projection_direction = up_direction, projection_targets_names=[nacelle_comp.name])
    thrust_origin = geo.perform_linear_interpolation(thrust_top_origin, thrust_bot_origin ,[1] ,output_parameters = np.array([0.5]))

    thrust_tip, thrust_tip_coord  = geo.project_points(thrust_tip_control_point, projection_direction = down_direction, projection_targets_names=[nacelle_comp.name])
    
    thrust_vector = geo.subtract_pointsets(thrust_tip, thrust_origin)

    geo.assemble()
    geo.evaluate()

    # print('_quarter_nacelle_length: ', _quarter_nacelle_length)
    # print('thrust_tip_control_point: ', thrust_tip_control_point)
    # print('thrust_top_origin_control_point[0] :', thrust_top_origin_control_point[0] + _quarter_nacelle_length)
    # print('Tip Point: ', thrust_tip_control_point)
    # print('thrust_top_origin_coord: ', thrust_top_origin_coord)
    # print('thrust_bot_origin_coord', thrust_bot_origin_coord)
    # print('Thrust Origin: ', thrust_origin.physical_coordinates)
    

    thrust_mag = np.linalg.norm(thrust_vector.physical_coordinates)
    thrust_vector_direction = geo.divide_pointset_by_scalar(thrust_vector, thrust_mag)

    # print('thrust_vector_direction: ', thrust_vector_direction)

    return [thrust_origin, thrust_vector_direction]
