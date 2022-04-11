''' import models and submodels'''
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
from csdl_om import Simulator
import matplotlib.pyplot as plt
from vedo import Points, Plotter, LegendBox
import pytest

# for debugging in Visual Studio Code only
# import os
# os.chdir("../lsdo_geo/lsdo_kit/design_geometry/tests/inner_opt_model_tests/geo_outputs_tests")

def test_displacement():
    ''' import file and create geometry object/ data structure '''
    path_name = '../../../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    ''' Creating the FFD block and objects '''
    from lsdo_kit.design.design_geometry.tests.organized_test.ffd_block_creation_script import generate_ffd_blocks
    ffd_blocks = generate_ffd_blocks(geo)

    ''' Add wing component to project points onto'''
    wing_comp = Component(stp_entity_names=['RectWing'], name='wing')
    geo.add_component(wing_comp)

    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    left_lead_point = np.array([0.0, -9000., 2000.])/1000
    left_trail_point = np.array([4000.0, -9000.0, 2000.])/1000
    right_lead_point = np.array([0.0, -9000.0, 2000.])/1000
    right_trail_point = np.array([4000.0, 9000.0, 2000.])/1000

    '''Project points'''
    wing_lead_left, wing_lead_left_coord = geo.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_left ,wing_trail_left_coord = geo.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_lead_right ,wing_lead_right_coord = geo.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_right ,wing_trail_right_coord = geo.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

    ''' Call the Geometric Outputs python class to create the dict for InnerOpt'''
    geo_outputs = GeometricOuputs(geo=geo)
    geo_outputs.compute_displacement('test_displacement', wing_lead_left, wing_trail_left)

    print('Geo dict: ', geo_outputs.geometric_outputs_dict)
    sim = Simulator(InnerOptimizationModel(ffd_blocks = ffd_blocks, geometry=geo, geometric_outputs=geo_outputs))
    sim.run()
    geo.evaluate()

    print('TEST DISPLACEMENT CALCULATIONS')
    print('wing_lead_left: ', wing_lead_left.physical_coordinates)
    print('wing_trail_left: ', wing_trail_left.physical_coordinates)
    print('true val displacement: ',        wing_lead_left.physical_coordinates - wing_trail_left.physical_coordinates)
    print('displacement: ',                 sim.prob['geo_outputs.test_displacement'])    
    print(' ')

    assert assert_allclose(sim.prob['geo_outputs.test_displacement'], wing_lead_left.physical_coordinates - wing_trail_left.physical_coordinates,atol=.1) == None

test_displacement()