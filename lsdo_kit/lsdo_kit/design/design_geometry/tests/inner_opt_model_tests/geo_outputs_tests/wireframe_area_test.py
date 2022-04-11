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
from lsdo_kit.design.design_geometry.wireframe_area_functions import *
import pytest
from csdl_om import Simulator
import matplotlib.pyplot as plt
from vedo import Points, Plotter, LegendBox
import os
# os.chdir("../lsdo_geo/lsdo_kit/design_geometry/tests/inner_opt_model_tests/geo_outputs_tests")

'''
Wireframe Area Test Script:

Calculates the area of a mesh through taking cross products of each mesh element.

* run "pytest wireframe_area_test.py" in the terminal to test this script

'''

def test_wireframe_area():
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
    geo_outputs.compute_wireframe_area('test_wireframe_area', camber_surface_mesh)

    ''' Simulate system '''
    print('Geo dict: ', geo_outputs.geometric_outputs_dict)                                                             # prints dictionary of computations
    sim = Simulator(InnerOptimizationModel(ffd_blocks = ffd_blocks, geometry=geo, geometric_outputs=geo_outputs))       # feed model into simulator class
    sim.run()           # run simulation
    geo.evaluate()      # evaluate pointsets -> enable pointsets.physical_coordinates to be called.       

    ''' comparing and testing csdl simulation with manual calculation of wireframe area '''
    print("TEST WIREFRAME AREA CALCULATIONS")
    print("Wireframe Area from simulation: ", sim.prob['geo_outputs.test_wireframe_area'])
    mesh_coordinates = mesh_coordinates_list(camber_surface_mesh.shape[0], camber_surface_mesh.shape[1], camber_surface_mesh.physical_coordinates)
    area_sum = get_area(mesh_coordinates)
    print("Manual Area Calculation: ", area_sum)
    assert assert_allclose(sim.prob['geo_outputs.test_wireframe_area'], area_sum, atol=.1)==None

test_wireframe_area()