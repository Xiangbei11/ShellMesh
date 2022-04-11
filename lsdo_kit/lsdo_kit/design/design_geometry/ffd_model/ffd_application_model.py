import csdl
import numpy as np 

from csdl_om import Simulator
from vedo import Points, Plotter, LegendBox

# TODO: Make it loop over FFD Blocks 

class FFDApplicationModel(csdl.Model):
    '''
    Maps from the FFD control points to the OML geometry control points.
    '''

    def initialize(self):
        # self.parameters.declare('ffd_blocks')
        # self.parameters.declare('output_geo_control_points_dict')
        self.parameters.declare('design_geometry_obj')
    def define(self):
        design_geometry_obj = self.parameters['design_geometry_obj']
        ffd_blocks = list(design_geometry_obj.components_ffd_dict.values())

        # ffd_blocks = self.parameters['ffd_blocks']
        # probably want collections.OrderedDict
        # output_geo_control_points_dict = self.parameters['output_geo_control_points_dict']
        total_geometry_control_points = np.array([])

        entity_end_index = 0
        for ffd_block in ffd_blocks: 
            
            ffd_map = ffd_block.project_points_FFD()

            ffd_map_var = self.create_input(f'{ffd_block.name}_map', val=ffd_map.todense())

            # The ffd_control_points are received from ffd_block_update_model
            ffd_control_points = self.declare_variable(f'{ffd_block.name}_updated_ffd_control_points', val=ffd_block.control_points)

            geometry_control_points = csdl.matmat(ffd_map_var, csdl.reshape(ffd_control_points, (ffd_block.control_points.shape[0] * ffd_block.control_points.shape[1] * ffd_block.control_points.shape[2], 3)))
            
            entity_starting_index = 0
            for entity in ffd_block.embedded_entities_pointers:
                entity_name = entity.name
                entity_end_index = entity_starting_index + np.cumprod(entity.shape)[-2]
                entity_control_points = geometry_control_points[entity_starting_index:entity_end_index, :]
                entity_starting_index = entity_end_index
                design_geometry_obj.output_geo_control_points_dict[entity_name] = entity_control_points

            self.register_output(f'{ffd_block.name}_geometry_control_points', geometry_control_points)

        # register one variable that is a concatenation of a bunch of other variables (e.g. variables stored in a dict)
        total_geometry_control_points = self.create_output('total_geometry_control_points', shape=(entity_end_index, 3))
        entity_starting_index = 0
        for gcp in design_geometry_obj.output_geo_control_points_dict.values():
            
            total_geometry_control_points[entity_starting_index:entity_starting_index+gcp.shape[0], :] = gcp
            entity_starting_index+=gcp.shape[0]

        # self.print_var(total_geometry_control_points)
        

if __name__ == "__main__":
    import csdl
    import numpy as np
    from lsdo_kit.geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 

    ''' Camber surface creation script for this case '''
    path_name = '../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    # wing_surface_names = [
    # 'Surf_WFWKRQIMCA, Wing, 0, 12', 'Surf_WFWKRQIMCA, Wing, 0, 13', 
    # 'Surf_WFWKRQIMCA, Wing, 0, 14', 'Surf_WFWKRQIMCA, Wing, 0, 15', 
    # ]
    
    # print(geo.input_bspline_entity_dict.values())

    ''' Creating the FFD Blocks and FFD Objects'''
    nxp = 5
    nyp = 5
    nzp = 5

    point000 = np.array([0., -9000., -200.])
    point001 = np.array([0., -9000., 200.])

    point100 = np.array([4000., -9000., -200.])
    point101 = np.array([4000., -9000., 200.])

    point010 = np.array([0., 9000., -200.])
    point011 = np.array([0., 9000., 200.])

    point110 = np.array([4000., 9000., -200.])
    point111 = np.array([4000., 9000., 200.])

    # point020 = np.array([10., 557., -6.])
    # point021 = np.array([10., 557., 10.])
    
    # point120 = np.array([43., 557., -6.])
    # point121 = np.array([43., 557., 10.])

    control_points = np.zeros((2,2,2,3))
    
    control_points[0,0,0,:] = point000
    control_points[0,0,1,:] = point001

    control_points[0,1,0,:] = point010
    control_points[0,1,1,:] = point011
    
    control_points[1,0,0,:] = point100
    control_points[1,0,1,:] = point101
    
    control_points[1,1,0,:] = point110
    control_points[1,1,1,:] = point111

    # control_points[0,2,0,:] = point020
    # control_points[0,2,1,:] = point021

    # control_points[1,2,0,:] = point120
    # control_points[1,2,1,:] = point121

    ffd_control_points = create_ffd(control_points/1000 * 1.2, nxp, nyp, nzp)

    ffd_block = FFD('ffd_block', ffd_control_points, embedded_entities_pointers=geo.input_bspline_entity_dict.values())
    ffd_block.plot(nxp, nyp, nzp)
    ffd_blocks=[ffd_block]

    # sim = Simulator(FFDapplicationModel(ffd_blocks=ffd_blocks))
    # sim.run()
    # sim.prob.model.list_inputs(prom_name=True)
    # sim.prob.model.list_outputs(prom_name=True)

    ''' Plotting results of Application model '''
    # pts_shape = ffd_block.control_points.shape 
    # nxp = pts_shape[0]
    # nyp = pts_shape[1]
    # nzp = pts_shape[2]

    # org_pts_reshape = np.reshape(ffd_block.control_points, (nxp * nyp * nzp, 3))
    # ffd_control_points = np.reshape(ffd_b)
    # # mod_pts_reshape = np.reshape(sim[f'{ffd_block.name}_geometry_control_points'], (2700, 3))

    # vp_init = Plotter()
    # vps = []
    # vps1 = Points(org_pts_reshape, r=8, c = 'blue')
    # # vps2 = Points(mod_pts_reshape, r=8, c='red')
    # vps.append(vps1)
    # # vps.append(vps2)

    # vp_init.show(vps, 'FFD Changes', axes=1, viewup="z", interactive = True)


    # sim.visualize_model()


    # ffd_block_1 = FFD('ffd_block_1', ffd_control_points)

    # # Create ffd_blocks list, and the csdl variables needed

    # ffd_blocks = [ffd_block_1]

    # sim = Simulator(CsdlBlockUpdateTest(ffd_block_list = ffd_blocks))
    # sim.run()

    # generate_ffd_plot(ffd_block_1.control_points, sim['Update.updated_ffd_control_points'])


    