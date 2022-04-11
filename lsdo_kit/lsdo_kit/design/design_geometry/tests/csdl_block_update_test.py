
import csdl
import numpy as np 

from lsdo_kit.geometry.ffd_block_update_model import FfdBlockUpdateModel
from csdl_om import Simulator

from lsdo_kit.geometry.utils.generate_ffd_plot import generate_ffd_plot

# TODO: Determine how to set values for properties that were not specified such that they equate to no
# change in the ffd_block_update model. 

class CsdlBlockUpdateTest(csdl.Model):

    def initialize(self):
        self.parameters.declare('ffd_block_list')

    def define(self):
        ffd_block_list = self.parameters['ffd_block_list']

        for ffd_block in ffd_block_list:
            
            change_trans = np.ones((ffd_block.nxp, )) * 50
            no_change_trans = np.zeros((ffd_block.nxp, ))
            
            change_rot = np.linspace(-1, 1, ffd_block.nxp)
            no_change_rot = np.zeros((ffd_block.nxp, ))

            change_scale = np.ones((ffd_block.nxp, )) * 5
            no_change_scale = np.ones((ffd_block.nxp, ))

            rot_x = self.create_input('{}_rot_x'.format(ffd_block.name), val=no_change_rot)
            rot_y = self.create_input('{}_rot_y'.format(ffd_block.name), val=no_change_rot)
            rot_z = self.create_input('{}_rot_z'.format(ffd_block.name), val=no_change_rot)

            trans_x = self.create_input('{}_trans_x'.format(ffd_block.name), val=no_change_trans)
            trans_y = self.create_input('{}_trans_y'.format(ffd_block.name), val=no_change_trans * 10)
            trans_z = self.create_input('{}_trans_z'.format(ffd_block.name), val=no_change_trans)

            scale_y = self.create_input('{}_scale_y'.format(ffd_block.name), val=change_scale)
            scale_z = self.create_input('{}_scale_z'.format(ffd_block.name), val=no_change_scale)
            
            # shape_var = self.create_input('{}_shape_var'.format(ffd_block.name), val=ffd_block.local_control_points)

            update_model = FfdBlockUpdateModel(ffd_blocks=ffd_block_list)
            self.add(update_model, name=f'{ffd_block.name}', promotes=[])

if __name__ == "__main__":
    import csdl
    import numpy as np
    from lsdo_kit.geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD

    nxp = 4
    nyp = 10
    nzp = 6

    point000 = np.array([170. ,0. ,100.])
    point001 = np.array([170., 0., 170.])

    point010 = np.array([170., 230., 100.])
    point011 = np.array([170., 230., 170.])
    
    point100 = np.array([1000. ,0. ,100.])
    point101 = np.array([1000. ,0. ,170.])

    point110 = np.array([1000. ,230. ,100.])
    point111 = np.array([1000. ,230. ,170.])

    control_points = np.zeros((2,2,2,3))
    
    control_points[0,0,0,:] = point000
    control_points[0,0,1,:] = point001

    control_points[0,1,0,:] = point010
    control_points[0,1,1,:] = point011
    
    control_points[1,0,0,:] = point100
    control_points[1,0,1,:] = point101
    
    control_points[1,1,0,:] = point110
    control_points[1,1,1,:] = point111

    ffd_control_points = create_ffd(control_points, nxp, nyp, nzp)

    ffd_block_1 = FFD('ffd_block_1', ffd_control_points)

    # Create ffd_blocks list, and the csdl variables needed

    ffd_blocks = [ffd_block_1]

    sim = Simulator(CsdlBlockUpdateTest(ffd_block_list = ffd_blocks))
    sim.run()
    sim.visualize_model()
    sim.prob.model.list_inputs(prom_name=True)
    sim.prob.model.list_outputs(prom_name=True)

    # generate_ffd_plot(ffd_block_1.control_points, sim['Update.updated_ffd_control_points'])


    