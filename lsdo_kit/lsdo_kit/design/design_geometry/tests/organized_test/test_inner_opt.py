import csdl
import numpy as np

from lsdo_kit.design.design_geometry.ffd_model.section_properties_model import SectionPropertiesModel 
from lsdo_kit.design.design_geometry.ffd_model.ffd_block_update_model import FFDBlockUpdateModel
from lsdo_kit.design.design_geometry.ffd_model.ffd_application_model import FFDApplicationModel
from lsdo_kit.design.design_geometry.mesh_evaluation_model import MeshEvaluationModel

from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD

import matplotlib.pyplot as plt
from vedo import Points, Plotter, LegendBox

from csdl_om import Simulator

# use only for debugging in vscode
# import os
# os.chdir("../lsdo_geo/lsdo_kit/design_geometry/tests/organized_test")

# Kinda wanna call this DesignOpt, but that's a bit misleading
# Should be an implicit model?
class InnerOpt(csdl.Model):


    def initialize(self):
        self.parameters.declare('ffd_blocks')
        self.parameters.declare('geometry')
        self.parameters.declare('meshes')

        
    def define(self):
        ffd_blocks = self.parameters['ffd_blocks']
        geometry = self.parameters['geometry']
        meshes = self.parameters['meshes']

        bogus_output = self.declare_variable('bogus_output', val=1.0)
        bogus_output_2 = bogus_output + 1.0
        self.register_output('bogus_output_2', bogus_output_2)

        section_properties = SectionPropertiesModel(ffd_blocks=ffd_blocks)
        self.add(section_properties, name='section_properties', promotes=[])

        block_update = FFDBlockUpdateModel(ffd_blocks=ffd_blocks)
        self.add(block_update, name='update', promotes=[])

        application_model = FFDApplicationModel(ffd_blocks=ffd_blocks, geo_control_points_dict=geometry.output_geo_control_points_dict)
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
    from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    from lsdo_kit.old_files.mesh import Mesh


    ''' Camber surface creation script for this case '''
    path_name = '../../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    ''' Plotting to determine surface names '''
    # for surface in geo.input_bspline_entity_dict.values():
    #     vp_init = Plotter()
    #     vps1 = Points(surface.control_points, r=8, c = 'blue')
    #     # vps.append(vps2)
    #     vp_init.show(vps1, f'{surface.name}', axes=1, viewup="z", interactive = True)

    # exit()

    ''' Creating the FFD block and objects '''
    from lsdo_kit.design.design_geometry.tests.organized_test.ffd_block_creation_script import generate_ffd_blocks
    ffd_blocks = generate_ffd_blocks(geo)
    

    ''' Creating the chamber line mesh for mesh eval '''
    from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes
    meshes = generate_meshes(geo)
    
    sim = Simulator(InnerOpt(ffd_blocks=ffd_blocks, geometry=geo, meshes=meshes))
    sim.run()

    # camber_surface = geo.evaluate(camber_surface_mesh)
    # camber_surface = sim.prob['mesh.vlm_mesh_camber_surface_geo']
    # chord_surface = sim.prob['mesh.vlm_mesh_chord_surface']

    # chord_surface = sim.prob['mesh.chord_surface_mesh']

    temp = ffd_blocks[0]
    original_ffd_pts = temp.control_points
    original_ffd_pts = np.reshape(original_ffd_pts, (-1,3))
    updated_ffd_pts = sim.prob['update.test_ffd1_updated_ffd_cp']

    print('updated_ffd_pts: ', updated_ffd_pts.shape)
    updated_ffd_pts = np.reshape(updated_ffd_pts, (-1,3))
    print('updated_ffd_pts_reshape', updated_ffd_pts.shape)
    # geo.fit_bspline_entities(pointset_list=[camber_surface_mesh], output_vec=camber_surface)

    application_cp_tot = sim.prob['application.total_geometry_control_points']
    application_cp_tot = np.reshape(application_cp_tot, (-1,3))
    vp_init = Plotter()
    vps1 = Points(application_cp_tot, r=8, c = 'blue')
    vp_init.show(vps1, 'Application Total DesignGeometry', axes=1, viewup="z", interactive = True)    

    application_cp_geo_cp = sim.prob['application.test_ffd1_geometry_control_points']
    application_cp_geo_cp = np.reshape(application_cp_geo_cp, (-1,3))
    vp_init = Plotter()
    vps1 = Points(application_cp_geo_cp, r=8, c = 'blue')
    vp_init.show(vps1, 'Application Geo Control Points', axes=1, viewup="z", interactive = True)    

    camber_surface = sim.prob['mesh.vlm_mesh_camber_surface']
    vp_init = Plotter()
    vps1 = Points(camber_surface, r=8, c = 'blue')
    # vps2 = Points(chord_surface, r=8, c = 'red')
    vp_init.show(vps1, 'Camber', axes=1, viewup="z", interactive = True)

    # vp_init = Plotter()
    # vps1 = Points(geo.total_cntrl_pts_vector, r=8, c = 'blue')
    # vp_init.show(vps1, 'Entire Surface', axes=1, viewup="z", interactive = True)
'''
    # vp_init = Plotter()
    # vps1 = Points(camber_surface, r=8, c = 'blue')
    # vps2 = Points(chord_surface, r=8, c = 'red')
    # vp_init.show(vps1, vps2, 'Camber', axes=1, viewup="z", interactive = True)
'''
    # vp_init = Plotter()
    # vps1 = Points(updated_ffd_pts, r=8, c = 'red')
    # vps2 = Points(original_ffd_pts, r=8, c = 'black')
    # # vps3 = Points(np.reshape(temp.embedded_points, (-1, 3)), r=8, c = 'orange')
    # vps = [vps1, vps2]
    # vp_init.show(vps, 'FFD Pts', axes=1, viewup="z", interactive = True)

    # sim.prob.model.list_inputs(prom_name=True)
    # sim.prob.model.list_outputs(prom_name=True)

    # sim.visualize_model()

