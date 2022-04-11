import csdl
import numpy as np

from lsdo_kit.geometry.section_properties_model import SectionPropertiesModel 
from lsdo_kit.geometry.ffd_block_update_model import FfdBlockUpdateModel
from lsdo_kit.geometry.ffd_application_model import FFDapplicationModel
from lsdo_kit.geometry.mesh_evaluation_model import MeshEvaluationModel

from lsdo_kit.geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD

from csdl_om import Simulator

class TopFFDModel(csdl.Model):

    def initialize(self):
        self.parameters.declare('ffd_blocks')
        self.parameters.declare('geometry')
        self.parameters.declare('meshes')

        
    def define(self):
        ffd_blocks = self.parameters['ffd_blocks']
        geometry = self.parameters['geometry']
        meshes = self.parameters['meshes']

        section_properties = SectionPropertiesModel(ffd_blocks=ffd_blocks)
        self.add(section_properties, name='section_properties', promotes=[])

        block_update = FfdBlockUpdateModel(ffd_blocks=ffd_blocks)
        self.add(block_update, name='update', promotes=[])

        application_model = FFDapplicationModel(ffd_blocks=ffd_blocks, geo_control_points_dict=geometry.output_geo_control_points_dict)
        self.add(application_model, name='application', promotes=[])

        mesh_eval = MeshEvaluationModel(geometry=geometry, meshes=meshes)
        self.add(mesh_eval, name='mesh', promotes=[])

        for ffd_block in ffd_blocks:
            self.connect(f'update.{ffd_block.name}_updated_ffd_cp', f'application.{ffd_block.name}_updated_ffd_control_points')
            for property_name in ffd_block.properties_list:
                self.connect(f'section_properties.{ffd_block.name}.{property_name}_sum', f'update.{ffd_block.name}_{property_name}')

        self.connect('application.total_geometry_control_points', 'mesh.geometry_control_points')

if __name__ == "__main__":
    import csdl
    import numpy as np
    from lsdo_kit.geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    from lsdo_kit.old_files.mesh import Mesh


    ''' Camber surface creation script for this case '''
    path_name = '../examples/CAD/'
    file_name = 'test_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    ''' Creating the FFD block and objects '''
    nxp = 5
    nyp = 7
    nzp = 5

    point000 = np.array([10., -557., -6.])
    point001 = np.array([10., -557., 10.])

    point100 = np.array([43., -557., -6.])
    point101 = np.array([43., -557., 10.])

    point010 = np.array([-10., 0., -6.])
    point011 = np.array([-10., 0., 10.])

    point110 = np.array([98., 0, -6.])
    point111 = np.array([98., 0, 10.])

    point020 = np.array([10., 557., -6.])
    point021 = np.array([10., 557., 10.])
    
    point120 = np.array([43., 557., -6.])
    point121 = np.array([43., 557., 10.])

    control_points = np.zeros((2,3,2,3))
    
    control_points[0,0,0,:] = point000
    control_points[0,0,1,:] = point001

    control_points[0,1,0,:] = point010
    control_points[0,1,1,:] = point011
    
    control_points[1,0,0,:] = point100
    control_points[1,0,1,:] = point101
    
    control_points[1,1,0,:] = point110
    control_points[1,1,1,:] = point111

    control_points[0,2,0,:] = point020
    control_points[0,2,1,:] = point021

    control_points[1,2,0,:] = point120
    control_points[1,2,1,:] = point121

    ffd_control_points = create_ffd(control_points/100/2.9, nxp, nyp, nzp)

    ffd_block_1 = FFD('test_ffd1', ffd_control_points, embedded_entities_pointers=geo.input_bspline_entity_dict.values())
    ffd_block_2 = FFD('test_ffd2', ffd_control_points, embedded_entities_pointers=geo.input_bspline_entity_dict.values())

    ''' Creating the chamber line mesh for mesh eval '''

    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    lead_root_point = np.array([0., 0., 2.])
    lead_tip_point = np.array([18.3, 556.26, 2.])/100/2.9
    trail_tip_point = np.array([42.67, 556.26, 2.])/100/2.9
    trail_root_point = np.array([97.53, 0., 2.])/100/2.9

    '''Project points'''
    wing_lead_root, wing_lead_root_coord = geo.project_points(lead_root_point, projection_direction = down_direction)
    wing_lead_tip ,wing_lead_tip_coord = geo.project_points(lead_tip_point, projection_direction = down_direction)
    wing_trail_tip ,wing_trail_tip_coord = geo.project_points(trail_tip_point, projection_direction = down_direction)
    wing_trail_root ,wing_trail_root_coord = geo.project_points(trail_root_point, projection_direction = down_direction)

    num_pts1 = [20]#20
    num_pts2 = [20]#5
    surface_shape1 = np.append(num_pts1,num_pts2)
    lead_edge_curve = geo.perform_linear_interpolation(wing_lead_root, wing_lead_tip, num_pts1)
    trail_edge_curve = geo.perform_linear_interpolation(wing_trail_root, wing_trail_tip, num_pts1)

    chord_surface = geo.perform_linear_interpolation(lead_edge_curve, trail_edge_curve, surface_shape1)
    geo.assemble(chord_surface)
    chord_surface_mesh = geo.evaluate(chord_surface)

    '''Project from chord to wing surfaces'''
    top_points, top_plot = geo.project_points(chord_surface_mesh, projection_direction=up_direction)
    bot_points, bot_plot = geo.project_points(chord_surface_mesh, projection_direction=down_direction)

    '''Create Camber Surface'''
    surface_shape2 = np.append(surface_shape1,1)
    output_parameters = 0.5*np.ones([num_pts1[0]*num_pts2[0],1])
    camber_surface_mesh = geo.perform_linear_interpolation(top_points, bot_points, surface_shape2, output_parameters = output_parameters)

    vlm_mesh = Mesh('vlm_mesh')
    vlm_mesh.add_pointset(camber_surface_mesh, name="camber_surface")   
    
     # geo.register_output(camber_surface_mesh, name = 'camber_surface')

    # # Create ffd_blocks list, and the csdl variables needed
    ffd_block_1.add_shape_parameter(property_name = 'rot_x', parameter_name='linear', order=2, num_cp=3, dv=False)
    ffd_block_1.add_shape_parameter(property_name = 'rot_x', parameter_name='quadratic', order=3, num_cp=4, dv=False)

    # ffd_block_2.add_shape_parameter(property_name ='rot_x', parameter_names=['linear', 'quadratic'], order=[2,3], num_cp=[3,5], dv=[False, False], val=[np.arange(3), np.arange(5)])

    ffd_blocks = [ffd_block_1, ffd_block_2]

    sim = Simulator(TopFFDModel(ffd_blocks=ffd_blocks, geometry=geo, meshes=[vlm_mesh]))
    sim.run()
    sim.prob.model.list_inputs(prom_name=True)
    sim.prob.model.list_outputs(prom_name=True)

    sim.visualize_implementation()

