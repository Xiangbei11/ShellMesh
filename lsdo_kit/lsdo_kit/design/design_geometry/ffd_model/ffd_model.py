import csdl
import numpy as np

from lsdo_kit.design.design_geometry.ffd_model.section_properties_model import SectionPropertiesModel 
from lsdo_kit.design.design_geometry.ffd_model.ffd_block_update_model import FFDBlockUpdateModel
from lsdo_kit.design.design_geometry.ffd_model.ffd_application_model import FFDApplicationModel
# from lsdo_kit.geometry.mesh_evaluation_model import MeshEvaluationModel

# import os
# os.chdir("../lsdo_geo/lsdo_kit/design_geometry/csdl/ffd_csdl")

class FFDModel(csdl.Model):

    def initialize(self):
        # self.parameters.declare('ffd_blocks')
        # self.parameters.declare('geometry')
        # self.parameters.declare('meshes')
        self.parameters.declare('design_geometry_obj')
        
    def define(self):
        # ffd_blocks = self.parameters['ffd_blocks']
        # geometry = self.parameters['geometry']
        # meshes = self.parameters['meshes']
        design_geometry_obj = self.parameters['design_geometry_obj']

        section_properties = SectionPropertiesModel(design_geometry_obj=design_geometry_obj)
        self.add(section_properties, name='section_properties', promotes=[])

        block_update = FFDBlockUpdateModel(design_geometry_obj=design_geometry_obj)
        self.add(block_update, name='update', promotes=[])

        application_model = FFDApplicationModel(design_geometry_obj=design_geometry_obj)
        self.add(application_model, name='application', promotes=[])

        ffd_blocks = list(design_geometry_obj.components_ffd_dict.values())
        
        for ffd_block in ffd_blocks:
            self.connect(f'update.{ffd_block.name}_updated_ffd_cp', f'application.{ffd_block.name}_updated_ffd_control_points')
            for property_name in ffd_block.properties_list:
                self.connect(f'section_properties.{ffd_block.name}.{property_name}_sum', f'update.{ffd_block.name}_{property_name}')



if __name__ == "__main__":
    import csdl
    import numpy as np
    from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    from lsdo_kit.old_files.mesh import Mesh

    from csdl_om import Simulator

    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    path_name = '../../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    ''' Creating the FFD block and objects '''
    from lsdo_kit.design.design_geometry.tests.organized_test.ffd_block_creation_script import generate_ffd_blocks
    ffd_blocks = generate_ffd_blocks(geo)

    sim = Simulator(FFDModel(ffd_blocks=ffd_blocks, geometry=geo))
    sim.run()
    sim.visualize_implementation()

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
