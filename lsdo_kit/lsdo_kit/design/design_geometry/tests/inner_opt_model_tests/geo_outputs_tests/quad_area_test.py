from lsdo_kit.design.design_geometry.inner_optimization_model import InnerOptimizationModel

# import os
# os.chdir("../lsdo_geo/lsdo_kit/design_geometry/tests/inner_opt_model_tests/geo_outputs_tests")

'''
Quadrilaterial Area Computation Unit Test: 

This objective of this test is to see if the quadrilateral area computation is correct. The area is computed through taking the 
cross product of two vectors (represented as pointsets) and taking the 2-norm.

'''
if __name__ == '__main__':
    import csdl
    import numpy as np
    from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    from lsdo_kit.old_files.mesh import Mesh
    from lsdo_kit.design.design_geometry.core.geometric_outputs import GeometricOuputs
    from lsdo_kit.design.design_geometry.core.component import Component

    from csdl_om import Simulator

    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    ''' import file and create geometry object/ data structure '''
    path_name = '../../../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)


    ''' Creating the FFD block and objects '''
    from lsdo_kit.design.design_geometry.tests.organized_test.ffd_block_creation_script import generate_ffd_blocks
    ffd_blocks = generate_ffd_blocks(geo)

    ''' Creating the camber line mesh for mesh eval '''
    from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes
    meshes, camber_surface_mesh = generate_meshes(geo)
    
    ''' Visualize geometry control points '''
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(geo.total_cntrl_pts_vector[:,0],
    geo.total_cntrl_pts_vector[:,1],
    geo.total_cntrl_pts_vector[:,2]) # plot the point (2,3,4) on the figure
    plt.show()

    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    ''' User specified points to project down onto the geometry '''
    left_lead_point = np.array([0.0, -9000., 2000.])/1000
    left_trail_point = np.array([4000.0, -9000.0, 2000.])/1000
    right_lead_point = np.array([0.0, -9000.0, 2000.])/1000
    right_trail_point = np.array([4000.0, 9000.0, 2000.])/1000

    '''Project user-specified points onto the imported geometry -> outputs are the projected points onto the geometry'''
    wing_lead_left, wing_lead_left_coord = geo.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_left ,wing_trail_left_coord = geo.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_lead_right ,wing_lead_right_coord = geo.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_right ,wing_trail_right_coord = geo.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

    ''' Call the Geometric Outputs python class to create the dict for InnerOpt -> this dictionary will specify what calculations the user wants'''
    geo_outputs = GeometricOuputs(geo=geo)
    geo_outputs.compute_quad_area('test_quad_area', wing_lead_left,wing_trail_left,wing_trail_right)

    ''' Simulate system '''
    print('Geo dict: ', geo_outputs.geometric_outputs_dict)                                                             # prints dictionary of computations
    sim = Simulator(InnerOptimizationModel(ffd_blocks = ffd_blocks, geometry=geo, geometric_outputs=geo_outputs))       # feed model into simulator class
    sim.run()

    ''' Compare results from CSDL versus manual calculation '''
    print("TEST QUAD AREA CALCULATIONS")
    print("Quad Area from simulation: ", sim.prob['geo_outputs.test_quad_area'])
    quad_area_test_vector_1 = wing_lead_left_coord-wing_trail_left_coord
    quad_area_test_vector_2 = wing_trail_right_coord - wing_trail_left_coord
    calculated_area = np.linalg.norm(np.cross(quad_area_test_vector_1,quad_area_test_vector_2))
    print("Quad Area from manual calculation: ", calculated_area)