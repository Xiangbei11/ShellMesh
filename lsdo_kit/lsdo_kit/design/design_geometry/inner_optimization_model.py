import csdl
import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

from lsdo_kit.design.design_geometry.ffd_model.ffd_model import FFDModel 
from lsdo_kit.design.design_geometry.geometric_outputs_model import GeometricOutputsModel 
from lsdo_kit.design.design_geometry.evaluate_pointsets import EvaluatePointsets
from lsdo_kit.design.design_geometry.wireframe_area_functions import *
# from lsdo_kit.design.design_geometry.geometric_outputs_model import GeometricOutputsModel 

# import os
# os.chdir("../lsdo_kit/lsdo_kit/design/design_geometry")

class InnerOptimizationModel(csdl.Model):

    def initialize(self):
        self.parameters.declare('ffd_blocks')
        self.parameters.declare('geometry')
        self.parameters.declare('geometric_outputs')
        # self.parameters.declare('meshes')

        
    def define(self):
        ffd_blocks = self.parameters['ffd_blocks']
        geometry = self.parameters['geometry']
        geometric_outputs = self.parameters['geometric_outputs']
        # meshes = self.parameters['meshes']

        ffd_model = FFDModel(ffd_blocks=ffd_blocks, geometry=geometry)
        self.add(ffd_model, name='ffd_model', promotes=[])

        # app_points = self.declare_variable('app_points')
        # self.connect(f'ffd_model.application.total_geometry_control_points', f'app_points') 

        eval_pointsets = EvaluatePointsets(geometry=geometry)
        self.add(eval_pointsets, name='eval_pointsets', promotes=[])
        self.connect(f'ffd_model.application.total_geometry_control_points', f'eval_pointsets.control_points')
        
        # geometric_outputs_model = GeometricOutputsModel(geometry=geometry, geometric_outputs=geometric_outputs)
        # self.add(geometric_outputs_model, name='geo_outputs', promotes=[])
        # self.connect('eval_pointsets.points', 'geo_outputs.points')
        # Then call geometric outputs model 


        ''' 
            Still need to discuss how geometric constraints will be handled
        '''


if __name__ == "__main__":
    import csdl
    import numpy as np
    from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    from lsdo_kit.design.design_geometry.core.mesh import Mesh
    from lsdo_kit.design.design_geometry.core.geometric_outputs import GeometricOuputs
    from lsdo_kit.design.design_geometry.design_geometry_primitive import Component

    from csdl_om import Simulator

    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    ''' import file and create geometry object/ data structure '''
    path_name = '../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    # ''' Creating the FFD block and objects '''
    from lsdo_kit.design.design_geometry.tests.organized_test.ffd_block_creation_script import generate_ffd_blocks
    ffd_blocks = generate_ffd_blocks(geo)

    ''' Creating the chamber line mesh for mesh eval '''
    from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes
    meshes, camber_surface_mesh = generate_meshes(geo)

    # ''' Call the Geometric Outputs python class to create the dict for InnerOpt'''
    # wing_comp = Component(name='wing', stp_entity_names=['RectWing'])
    # # c._add_GeometryPrimitve_to_vehicle()
    # geo.add_component(wing_comp)

    # geo.assemble()
    # geo.evaluate()
    # fig = plt.figure(figsize=(4,4))
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(geo.total_cntrl_pts_vector[:,0],
    # geo.total_cntrl_pts_vector[:,1],
    # geo.total_cntrl_pts_vector[:,2])
    # plt.show()

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

    geo_outputs = GeometricOuputs(geo=geo)
    # geo_outputs.compute_displacement('test_displacement', wing_lead_left, wing_trail_left)
    # geo_outputs.compute_magnitude('test_mag', wing_lead_left, wing_trail_left)
    # geo_outputs.compute_angle('test_angle', wing_lead_right, wing_trail_right)
    geo_outputs.compute_wireframe_area('test_wireframe_area', camber_surface_mesh)
    # geo_outputs.compute_quad_area('test_quad_area',wing_lead_left,wing_trail_left,wing_trail_right)
    # geo_outputs.compute_elliptical_area('test_elliptical_area', 2,3)

    print('Geo dict: ', geo_outputs.geometric_outputs_dict)
    sim = Simulator(InnerOptimizationModel(ffd_blocks = ffd_blocks, geometry=geo, geometric_outputs=geo_outputs))
    sim.run()
    geo.evaluate()

    print("TEST WIREFRAME AREA CALCULATIONS")
    print("Wireframe Area from simulation: ", sim.prob['geo_outputs.test_wireframe_area'])
    mesh_coordinates = mesh_coordinates_list(camber_surface_mesh.shape[0], camber_surface_mesh.shape[1], camber_surface_mesh.physical_coordinates)
    area_sum = get_area(mesh_coordinates)
    print("Manual Area Calculation: ", area_sum)

    # print("TEST QUAD AREA CALCULATIONS")
    # print("Quad Area from simulation: ", sim.prob['geo_outputs.test_quad_area'])
    # quad_area_test_vector_1 = wing_lead_left_coord-wing_trail_left_coord
    # quad_area_test_vector_2 = wing_trail_right_coord - wing_trail_left_coord
    # calculated_area = np.linalg.norm(np.cross(quad_area_test_vector_1,quad_area_test_vector_2))
    # print("Quad Area from manual calculation: ", calculated_area)

    # print('TEST DISPLACEMENT CALCULATIONS')
    # print('wing_lead_left: ', wing_lead_left.physical_coordinates)
    # print('wing_trail_left: ', wing_trail_left.physical_coordinates)
    # print('true val displacement: ', wing_lead_left.physical_coordinates - wing_trail_left.physical_coordinates)
    # print('displacement: ', sim.prob['geo_outputs.test_displacement'])    
    # print(' ')

    # print('TEST MAGNITUDE CALCULATIONS')
    # print('wing_lead_left: ', wing_lead_left.physical_coordinates)
    # print('wing_trail_left: ', wing_trail_left.physical_coordinates)
    # temp =  wing_lead_left.physical_coordinates - wing_trail_left.physical_coordinates    
    # print('True Value magnitude: ', np.linalg.norm(temp))
    # print('Model magnitude: ', sim.prob['geo_outputs.test_mag'])    
    # print(' ')

    # print('TEST ANGLE CALCULATIONS')
    # print('INMODEL point1: ', sim.prob['geo_outputs.points1'])
    # print('INMODEL point2: ', sim.prob['geo_outputs.points2'])
    
    # print('wing_lead_left: ', wing_lead_right.physical_coordinates)
    # print('wing_trail_right: ', wing_trail_right.physical_coordinates)
    # p1 = wing_lead_right.physical_coordinates; p2 = wing_trail_right.physical_coordinates
    # p1 = np.reshape(p1, (3,)) ; p2 = np.reshape(p2, (3,)) 
    # print('post reshape lead left: ', p1)
    # print('post reshape trail right: ', p2)
    # pt1unit = np.linalg.norm(p1)
    # pt2unit = np.linalg.norm(p2)

    # num_dot = np.dot(p1, p2)
    # den_dot = pt1unit * pt2unit
    # ptangle = np.arccos(num_dot / den_dot)

    # print('true numerator: ', num_dot)
    # print('model numerator: ', sim.prob['geo_outputs.num_dot'])
    # print('true denominator: ', den_dot)
    # print('model denominator: ', sim.prob['geo_outputs.den_prod'])
    # print('true angle: ', ptangle)
    # print('model angle: ', sim.prob['geo_outputs.test_angle'])

    # sim.prob.model.list_inputs(prom_name=True)
    # sim.prob.model.list_outputs(prom_name=True, print_arrays=True)
    # sim.visualize_implementation()


    # up_direction = np.array([0., 0., 1.])
    # down_direction = np.array([0., 0., -1.])

    # # TODO Definitely want to automate this for all types of components
    # #  could look into using normal vectors to distinguish between up and down
    # wing_surface_names = [
    #     'RectWing, 0, 3', 'RectWing, 1, 8', 
    #     'RectWing, 0, 2', 'RectWing, 1, 9', 
    #     ]
    # top_wing_surface_names = [
    #     'RectWing, 0, 3', 
    #     'RectWing, 1, 9', 
    #     ]
    # bot_wing_surface_names = [
    #     'RectWing, 0, 2',
    #     'RectWing, 1, 8', 
    #     ]

    # # TODO Make helper function for camber surface
    # ''' Points to be projected'''
    # lead_port_point = np.array([0., -9., 3.])
    # lead_sb_point = np.array([0., 9., 3.])
    # trail_sb_point = np.array([4., 9., 3.])
    # trail_port_point = np.array([4., -9., 3.])

    # '''Project points'''
    # wing_lead_port, wing_lead_port_coord = geo.project_points(lead_port_point, projection_targets_names=top_wing_surface_names, projection_direction = down_direction)
    # wing_lead_sb ,wing_lead_sb_coord = geo.project_points(lead_sb_point, projection_targets_names=top_wing_surface_names, projection_direction = down_direction)
    # wing_trail_sb ,wing_trail_sb_coord = geo.project_points(trail_sb_point, projection_targets_names=top_wing_surface_names, projection_direction = down_direction)
    # wing_trail_port ,wing_trail_port_coord = geo.project_points(trail_port_point, projection_targets_names=top_wing_surface_names,projection_direction = down_direction)


