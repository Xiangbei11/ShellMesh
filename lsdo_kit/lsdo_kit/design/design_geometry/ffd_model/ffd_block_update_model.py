import csdl
from csdl_om import Simulator
import numpy as np

from lsdo_kit.design.design_geometry.utils.calculate_rotation_mat import calculate_rotation_mat

class FFDBlockUpdateModel(csdl.Model):
    '''
    Maps from the section properties to the FFD control points.
    '''

    def initialize(self):
        # self.parameters.declare('ffd_blocks')
        self.parameters.declare('design_geometry_obj')
    def define(self):
 
        design_geometry_obj = self.parameters['design_geometry_obj']
        ffd_blocks = list(design_geometry_obj.components_ffd_dict.values())
      
        for ffd_block in ffd_blocks: 
            # These rotations, translations, and scalings all come from the section properties model.
            # These values below are not matrices, they are arrays. 

            # Add a val argument such that if no connections are made, initial values correspond to no changes
            rot_x = self.declare_variable(f'{ffd_block.name}_rot_x', shape=(ffd_block.nxp,), val=np.zeros((ffd_block.nxp, )))
            rot_y = self.declare_variable(f'{ffd_block.name}_rot_y', shape=(ffd_block.nxp,), val=np.zeros((ffd_block.nxp, )))
            rot_z = self.declare_variable(f'{ffd_block.name}_rot_z', shape=(ffd_block.nxp,), val=np.zeros((ffd_block.nxp, )))

            trans_x = self.declare_variable(f'{ffd_block.name}_trans_x', shape=(ffd_block.nxp,), val=np.zeros((ffd_block.nxp, )))
            trans_y = self.declare_variable(f'{ffd_block.name}_trans_y', shape=(ffd_block.nxp,), val=np.zeros((ffd_block.nxp, )))
            trans_z = self.declare_variable(f'{ffd_block.name}_trans_z', shape=(ffd_block.nxp,), val=np.zeros((ffd_block.nxp, )))

            scale_y = self.declare_variable(f'{ffd_block.name}_scale_y', shape=(ffd_block.nxp,), val=np.ones((ffd_block.nxp, )))
            scale_z = self.declare_variable(f'{ffd_block.name}_scale_z', shape=(ffd_block.nxp,), val=np.ones((ffd_block.nxp, )))
            
            #shape_var = self.declare_variable(f'{ffd_block.name}_shape', shape=ffd_block.control_points.shape, val=ffd_block.local_control_points)
            shape_var = self.declare_variable(f'{ffd_block.name}_shape', shape=ffd_block.control_points.shape, val=np.zeros_like(ffd_block.local_control_points))
            initial_shape = self.create_input(f'{ffd_block.name}_initial_shape', val=ffd_block.local_control_points)

            # ''' Bulk movement '''
            # bulk_trans_x = self.declare_variable(f'{ffd_block.name}_bulk_trans_x', shape=1, val=0)
            # bulk_trans_y = self.declare_variable(f'{ffd_block.name}_bulk_trans_y', shape=1, val=0)
            # bulk_trans_z = self.declare_variable(f'{ffd_block.name}_bulk_trans_z', shape=1, val=0)
            # bulk_rot_x = self.declare_variable(f'{ffd_block.name}_bulk_rot_x', shape=1, val=0)
            # bulk_rot_y = self.declare_variable(f'{ffd_block.name}_bulk_rot_y', shape=1, val=0)
            # bulk_rot_z = self.declare_variable(f'{ffd_block.name}_bulk_rot_z', shape=1, val=0)

            
            # local_ffd_points = self.create_input('local_ffd_points', val=ffd_block.local_control_points)                        
            
            temp = np.eye(3)

            # test_cos = csdl.expand(csdl.cos(rot_x), (ffd_block.nxp, 1, 1), 'i->ijk')

            # These test values below are being used for troubleshooting 
            # test1 = self.create_input(f'{ffd_block.name}_test1', val= -5)
            # test2 = self.create_input(f'{ffd_block.name}_test2', val= -10)
            # test3 = self.create_input(f'{ffd_block.name}_test3', val= 2)
            # test4 = self.create_input(f'{ffd_block.name}_test4', val= 8)

            rot_mat_x = self.create_output(f'{ffd_block.name}_rot_mat_x', shape=(ffd_block.nxp, 3, 3), val=np.repeat(temp[np.newaxis, :, :], ffd_block.nxp, axis=0))
            rot_mat_x[:, 1, 1] = csdl.expand(csdl.cos(rot_x), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_x[:, 1, 2] = csdl.expand(-csdl.sin(rot_x), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_x[:, 2, 1] = csdl.expand(csdl.sin(rot_x), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_x[:, 2, 2] = csdl.expand(csdl.cos(rot_x), (ffd_block.nxp, 1, 1), 'i->ijk')

            rot_mat_y = self.create_output(f'{ffd_block.name}_rot_mat_y', shape=(ffd_block.nxp, 3, 3), val=np.repeat(temp[np.newaxis, :, :], ffd_block.nxp, axis=0))
            rot_mat_y[:, 0, 0] = csdl.expand(csdl.cos(rot_y), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_y[:, 0, 2] = csdl.expand(-csdl.sin(rot_y), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_y[:, 2, 0] = csdl.expand(csdl.sin(rot_y), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_y[:, 2, 2] = csdl.expand(csdl.cos(rot_y), (ffd_block.nxp, 1, 1), 'i->ijk')

            rot_mat_z = self.create_output(f'{ffd_block.name}_rot_mat_z', shape=(ffd_block.nxp, 3, 3), val=np.repeat(temp[np.newaxis, :, :], ffd_block.nxp, axis=0))
            rot_mat_z[:, 0, 0] = csdl.expand(csdl.cos(rot_z), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_z[:, 0, 1] = csdl.expand(csdl.sin(rot_z), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_z[:, 1, 0] = csdl.expand(-csdl.sin(rot_z), (ffd_block.nxp, 1, 1), 'i->ijk')
            rot_mat_z[:, 1, 1] = csdl.expand(csdl.cos(rot_z), (ffd_block.nxp, 1, 1), 'i->ijk')
            
            # Order of application: yaw, pitch, roll 
            T =  rot_mat_x * rot_mat_y * rot_mat_z
            
            scl = self.create_output('{}_scl_mat'.format(ffd_block.name), shape=shape_var.shape, val=np.ones(shape_var.shape))
            scl[:,:,:,1] = csdl.expand(scale_y * ffd_block.initial_scale_y, (ffd_block.nxp, ffd_block.nyp, ffd_block.nzp, 1), 'i->ijkl') 
            scl[:,:,:,2] = csdl.expand(scale_z * ffd_block.initial_scale_z, (ffd_block.nxp, ffd_block.nyp, ffd_block.nzp, 1), 'i->ijkl')

            pos = self.create_output('{}_pos'.format(ffd_block.name), shape=ffd_block.origin.shape)
            pos_x = self.create_input(f'{ffd_block.name}_pos_x', val=ffd_block.origin[:, 0, np.newaxis])
            pos_y = self.create_input(f'{ffd_block.name}_pos_y', val=ffd_block.origin[:, 1, np.newaxis])
            pos_z = self.create_input(f'{ffd_block.name}_pos_z', val=ffd_block.origin[:, 2, np.newaxis])

            pos[:, 0] = pos_x + csdl.expand(trans_x, (ffd_block.nxp, 1), 'i->ij')
            pos[:, 1] = pos_y + csdl.expand(trans_y, (ffd_block.nxp, 1), 'i->ij')
            pos[:, 2] = pos_z + csdl.expand(trans_z, (ffd_block.nxp, 1), 'i->ij')
           
            pos_mat = self.register_output(f'{ffd_block.name}_pos_mat', csdl.expand(pos, (ffd_block.nxp, ffd_block.nyp, ffd_block.nzp, 3), 'ij->iklj'))

            new_shape = shape_var + initial_shape

            # scl_points = shape_var * scl 
            scl_points = new_shape * scl   # Apply new scaling
            # self.register_output('scl_points', scl_points)

            T_transposed = csdl.einsum(T, subscripts='kij->kji')    # Apply new rotation (delta) and original rotation

            # Apply delta translation and original position
            final_design_cp = scl_points + pos_mat
            # print(csdl.einsum(T_transposed, scl_points, subscripts='kji,klmj->klmi').shape)

            # Apply bulk movements
            # cp_final = R_bulk*(final_design_cp - origin) + bulk_trans + origin
            # updated_ffd_cp = TODO

            # self.register_output('pointset_points', csdl.einsum(T_transposed, scl_points, subscripts='kji,klmj->klmi') + pos_mat)
            self.register_output(f'{ffd_block.name}_updated_ffd_cp', final_design_cp)


if __name__ == "__main__":
    import csdl
    import numpy as np
    from lsdo_kit.geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.tests.csdl_block_update_test import CsdlBlockUpdateTest

    nxp = 5
    nyp = 5
    nzp = 5

    point000 = np.array([170. ,0. ,100.])
    point010 = np.array([130., 230., 100.])
    point001 = np.array([170., 0., 170.])
    point011 = np.array([130., 230., 170.])
    
    point100 = np.array([240. ,0. ,100.])
    point101 = np.array([240. ,0. ,170.])
    point110 = np.array([200. ,230. ,100.])
    point111 = np.array([200. ,230. ,170.])

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

    # Create ffd_blocks list, and the csdl variables needed

    ffd_blocks = [ffd_control_points]

    sim = Simulator(CsdlBlockUpdateTest())
    sim.run()
    sim.visualize_model()
