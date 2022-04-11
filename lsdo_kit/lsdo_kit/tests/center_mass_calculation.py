import numpy as np

def generate_CoM(geo, comp, CoMcoord, ip=0.5):

    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    _top_coord = np.array([0.0, 0.0, 0.0])
    _bot_coord = np.array([0.0, 0.0, 0.0])

    _top_coord[:] = CoMcoord
    _bot_coord[:] = CoMcoord

    # print('CoMcoord: ', CoMcoord)
    # print('_top_coord: ', _top_coord)
    # print('_bot_coord: ', _bot_coord)

    _top_coord[2] = comp.z_max + 1.0
    _bot_coord[2] = comp.z_min - 1.0

    # print('_top_coord new: ', _top_coord)
    # print('_bot_coord new: ', _bot_coord)

    _top_pt, _top_pt_coord = geo.project_points(_top_coord, projection_direction=down_direction, projection_targets_names=[comp.name])
    _bot_pt, _bot_pt_coord = geo.project_points(_bot_coord, projection_direction=up_direction, projection_targets_names=[comp.name])

    CoM = geo.perform_linear_interpolation(_top_pt, _bot_pt ,[1] ,output_parameters = np.array([ip]))

    return _top_pt_coord, _bot_pt_coord, CoM

