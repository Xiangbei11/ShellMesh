from lsdo_kit.design.design_geometry.inner_optimization_model import InnerOptimizationModel
import csdl
import numpy as np
from numpy.testing import assert_allclose
from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD
from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
from lsdo_kit.old_files.mesh import Mesh
from lsdo_kit.design.design_geometry.core.geometric_outputs import GeometricOuputs
from lsdo_kit.design.design_geometry.core.component import Component
import pytest
from csdl_om import Simulator
import matplotlib.pyplot as plt
from vedo import Points, Plotter, LegendBox
# import os
# os.chdir("../lsdo_geo/lsdo_kit/design_geometry/tests/inner_opt_model_tests/geo_outputs_tests")

'''
Angle between two pointsets unit test:

Let A and B be two pointsets. For our purposes, we can define them as vectors in R3. This test will compute the angle between 
the A and B will begin by dividing the dot product between A and B with the 2-norm of A times 2-norm of B. Taking the arccos
of this quantity will give the angle between the two vectors.

* run "pytest angle_test.py" in the terminal to test this script

'''

def test_angle():
    ''' import file and create geometry object/ data structure '''
    path_name = '../../../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    ''' Creating the FFD block and objects '''
    from lsdo_kit.design.design_geometry.tests.organized_test.ffd_block_creation_script import generate_ffd_blocks
    ffd_blocks = generate_ffd_blocks(geo)

    ''' Creating the chamber line mesh for mesh eval '''
    from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes
    meshes, camber_surface_mesh = generate_meshes(geo)

    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    left_lead_point = np.array([0.0, -9000., 2000.])/1000
    left_trail_point = np.array([4000.0, -9000.0, 2000.])/1000
    right_lead_point = np.array([0.0, -9000.0, 2000.])/1000
    right_trail_point = np.array([4000.0, 9000.0, 2000.])/1000

    '''Project user-specified points onto the imported geometry'''
    wing_lead_left, wing_lead_left_coord = geo.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_left ,wing_trail_left_coord = geo.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_lead_right ,wing_lead_right_coord = geo.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_right ,wing_trail_right_coord = geo.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

    ''' Call the Geometric Outputs python class to create the dict for InnerOpt -> this dictionary will specify what calculations the user wants'''
    geo_outputs = GeometricOuputs(geo=geo)
    geo_outputs.compute_angle('test_angle', wing_lead_right, wing_trail_right)

    ''' Simulate system '''
    print('Geo dict: ', geo_outputs.geometric_outputs_dict)                                                             # prints dictionary of computations
    sim = Simulator(InnerOptimizationModel(ffd_blocks = ffd_blocks, geometry=geo, geometric_outputs=geo_outputs))       # feed model into simulator class
    sim.run()
    # geo.evaluate()

    print('')
    print('TEST ANGLE CALCULATIONS')
    print('')
    print('=========================================')
    print('INMODEL point1: ', sim.prob['geo_outputs.points1'])
    print('INMODEL point2: ', sim.prob['geo_outputs.points2'])
    
    ''' Calculating numerators, denominators, and norms for angle calculation '''

    print('wing_lead_left: ', wing_lead_right.physical_coordinates)
    print('wing_trail_right: ', wing_trail_right.physical_coordinates)
    p1 = wing_lead_right.physical_coordinates; p2 = wing_trail_right.physical_coordinates
    p1 = np.reshape(p1, (3,)) ; p2 = np.reshape(p2, (3,)) 
    print('post reshape lead left: ', p1)
    print('post reshape trail right: ', p2)
    pt1unit = np.linalg.norm(p1)
    pt2unit = np.linalg.norm(p2)
    num_dot = np.dot(p1, p2)
    den_dot = pt1unit * pt2unit
    ptangle = np.arccos(num_dot / den_dot)

    print('=========================================')
    print('true numerator: ', num_dot)
    print('model numerator: ', sim.prob['geo_outputs.num_dot'])
    print('=========================================')
    print('true denominator: ', den_dot)
    print('model denominator: ', sim.prob['geo_outputs.den_prod'])
    print('=========================================')
    print('true angle: ', ptangle)
    print('model angle: ', sim.prob['geo_outputs.test_angle'])
    print('=========================================')

    assert assert_allclose(sim.prob['geo_outputs.test_angle'], ptangle, atol=.1)==None

test_angle()