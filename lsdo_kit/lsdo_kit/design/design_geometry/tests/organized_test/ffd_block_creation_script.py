import csdl
import numpy as np
from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD

# from lsdo_kit.tests.organized_test.test_inner_opt import geo

def generate_ffd_blocks(geo):
    ''' Creating the FFD block and objects '''

    ''' Created based on the coordinates from rect_wing.stp '''
    nxp = 5
    nyp = 7
    nzp = 5

    pt000 = np.array([0.0, -9500., -2000.])
    pt001 = np.array([0.0, -9500., 2000.])

    pt010 = np.array([0.0, 9500., -2000.])
    pt011 = np.array([0.0, 9500., 2000.])

    pt100 = np.array([4000.0, -9500., -2000.])
    pt101 = np.array([4000.0, -9500., 2000.])
    
    pt110 = np.array([4000.0, 9500., -2000.])
    pt111 = np.array([4000.0, 9500., 2000.])

    control_points = np.zeros((2,2,2,3))

    control_points[0,0,0,:] = pt000
    control_points[0,0,1,:] = pt001

    control_points[0,1,0,:] = pt010
    control_points[0,1,1,:] = pt011

    control_points[1,0,0,:] = pt100
    control_points[1,0,1,:] = pt101

    control_points[1,1,0,:] = pt110
    control_points[1,1,1,:] = pt111

    ffd_control_points = create_ffd(control_points/1000, nxp, nyp, nzp)

    ffd_block_1 = FFD('test_ffd1', ffd_control_points, embedded_entities_pointers=geo.input_bspline_entity_dict.values())
    ffd_block_2 = FFD('test_ffd2', ffd_control_points, embedded_entities_pointers=geo.input_bspline_entity_dict.values())

    # # Create ffd_blocks list, and the csdl variables needed
    # ffd_block_1.add_shape_parameter(property_name = 'rot_x', parameter_names=['linear'], order=[2], num_cp=[3], dv=[False])
    ffd_block_1.add_shape_parameter(property_name = 'scale_y', parameter_name='linear', order=1, num_cp=2, dv=False, val=np.ones(2))
    # ffd_block_1.add_shape_parameter(property_name = 'scale_y', parameter_name='quadratic', order=3, num_cp=4, dv=False)
    
    # ffd_block_1.add_shape_parameter('linear', 'scale_y', order=1, num_cp=2, dv=False, val=np.ones(2))
    # ffd_block_1.add_shape_parameter('quadratic', 'scale_y', order=3, num_cp=4, dv=False, val=np.ones(4)*2)

    # ffd_block_1.add_shape_parameter(property_name = 'trans_x', parameter_names=['linear', 'quadratic'], order=[2, 3], num_cp=[3, 4], dv=[True, False], val=[np.arange(3), None])

    # ffd_block_2.add_shape_parameter(property_name ='rot_x', parameter_names=['linear', 'quadratic'], order=[2,3], num_cp=[3,5], dv=[False, False], val=[np.arange(3), np.arange(5)])

    # ffd_blocks = [ffd_block_1, ffd_block_2]
    ffd_blocks = [ffd_block_1]
    return ffd_blocks


    ''' Points associated with test_wing.stp. Divide the points by 100/2.9 '''
    # point000 = np.array([20., -557., -6.])
    # point001 = np.array([20., -557., 10.])

    # point100 = np.array([43., -557., -6.])
    # point101 = np.array([43., -557., 10.])

    # point010 = np.array([-20., 0., -6.])
    # point011 = np.array([-20., 0., 10.])

    # point110 = np.array([98., 0, -6.])
    # point111 = np.array([98., 0, 10.])

    # point020 = np.array([20., 557., -6.])
    # point021 = np.array([20., 557., 10.])

    # point120 = np.array([43., 557., -6.])
    # point121 = np.array([43., 557., 10.])

    # control_points = np.zeros((2,3,2,3))

    # control_points[0,0,0,:] = point000
    # control_points[0,0,1,:] = point001

    # control_points[0,1,0,:] = point010
    # control_points[0,1,1,:] = point011

    # control_points[1,0,0,:] = point100
    # control_points[1,0,1,:] = point101

    # control_points[1,1,0,:] = point110
    # control_points[1,1,1,:] = point111

    # control_points[0,2,0,:] = point020
    # control_points[0,2,1,:] = point021

    # control_points[1,2,0,:] = point120
    # control_points[1,2,1,:] = point121