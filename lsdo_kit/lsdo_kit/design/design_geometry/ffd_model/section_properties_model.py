from lsdo_kit.design.design_geometry.bsplines.bspline_curve import BSplineCurve
import csdl
import numpy as np

from csdl.core.input import Input
# property_names_list = ['rot_x', 'rot_y', 'rot_z', 'trans_x', ' trans_y', 'trans_z', 'scale_y', 'scale_z']

'''
wing_component.add_shape_parameter('twist', num=2, order=2)
'''

class SectionPropertiesModel(csdl.Model):

    def initialize(self):
        # self.parameters.declare('ffd_blocks', types=list)
        self.parameters.declare('design_geometry_obj')
    def define(self):
        '''
        TOP LEVEL
        parameter_names = dict()
            - property_names are the keys of this top level dictionary
            - The property_names keys contain a dictionary that has parameter_name = string, and parameter_dict = dict()
            - parameter_dict has two keys order and num_cp
        
        '''

        # ffd_blocks = self.parameters['ffd_blocks']
        design_geometry_obj = self.parameters['design_geometry_obj']
        ffd_blocks = list(design_geometry_obj.components_ffd_dict.values())

        # print('Section Properties FFD Block: ', ffd_blocks[0])
        # print('Section Properties FFD Block: ', ffd_blocks[0].param_dict)
        

        for ffd_block in ffd_blocks:
            num_sections = ffd_block.nxp

            ffd_block_csdl_model = csdl.Model()

            parameter_products = []
            # print(parameter_products)
            for property_name, property_dict in ffd_block.properties_dict.items():

                # property_var = ffd_block_csdl_model.create_input(property_name)
                # parametes_list = ffd_block.properties_dict[property_name]
                # print(property_name)
                # print(property_dict)
                # for parameter_dict in parameters_list['parameters']:
                for parameter_name, parameter_info in property_dict['parameters'].items():
                    # print('parameter info: ', parameter_info['val'])
                    # print('parameter info type: ', type(parameter_info['val']))

                    # TODO: Add a design variable check, need two more if clauses!!

                    if parameter_info['val'] is not None and parameter_info['dv'] is not False:
                        parameter_var = ffd_block_csdl_model.create_input(property_name + '_' + parameter_name, shape=(parameter_info['num'],), val=parameter_info['val'])
                        ffd_block_csdl_model.add_design_variable(property_name + '_' + parameter_name)

                    elif parameter_info['val'] is not None and parameter_info['dv'] is not True:
                        parameter_var = ffd_block_csdl_model.create_input(property_name + '_' + parameter_name, shape=(parameter_info['num'],), val=parameter_info['val'])

                    elif parameter_info['val'] is None:
                        parameter_var = ffd_block_csdl_model.declare_variable(property_name + '_' + parameter_name, shape=(parameter_info['num'],))
                        # print('VAL: ', parameter_var.val)
                        # print('TYPE: ', type(parameter_var.val))
                        
                    # # TODO: Implement way to create initial values of Control Points for the parameter curve. Create default control points for shape 
                    # if parameter_info['num'] != None: 
                    #     control_points = np.zeros(parameter_info['num'])                       
                    
                    bsp_curve = BSplineCurve(name=parameter_name + '_curve', order_u=parameter_info['order'], control_points=parameter_var.val)
                    sp_mtx = bsp_curve.compute_eval_map_points(np.linspace(0., 1., num_sections))

                    # We store the products in a list to sum them later
                    parameter_products = parameter_products + [csdl.matvec(sp_mtx, parameter_var)]
                    # print(f'{parameter_name}: ', len(parameter_products))

                # print(parameter_list)
                if len(parameter_products) == 1:
                    ffd_block_csdl_model.register_output(property_name + '_sum', *parameter_products)
                else:
                    # print('Parameter: ', parameter_products)
                    ffd_block_csdl_model.register_output(property_name + '_sum', csdl.sum(*parameter_products))
            self.add(ffd_block_csdl_model, name=ffd_block.name, promotes=[])
        

if __name__ == "__main__":
    import csdl
    import numpy as np

    from lsdo_kit.geometry.section_properties_model import SectionPropertiesModel 
    from lsdo_kit.geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD

    from csdl_om import Simulator

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

    ffd_control_points1 = create_ffd(control_points, nxp, nyp, nzp)

    point000 = np.array([90. ,0. ,100.])
    point010 = np.array([180., 230., 100.])
    point001 = np.array([90., 0., 170.])
    point011 = np.array([180., 230., 170.])

    point100 = np.array([130. ,0. ,100.])
    point101 = np.array([130. ,0. ,170.])
    point110 = np.array([70. ,230. ,100.])
    point111 = np.array([70. ,230. ,170.])

    control_points = np.zeros((2,2,2,3))

    control_points[0,0,0,:] = point000
    control_points[0,0,1,:] = point001

    control_points[0,1,0,:] = point010
    control_points[0,1,1,:] = point011

    control_points[1,0,0,:] = point100
    control_points[1,0,1,:] = point101

    control_points[1,1,0,:] = point110
    control_points[1,1,1,:] = point111

    ffd_control_points2 = create_ffd(control_points, nxp, nyp, nzp)

    ffd_block_1 = FFD('test_ffd1', ffd_control_points1)
    ffd_block_2 = FFD('test_ffd2', ffd_control_points2)

    ffd_block_1.add_shape_parameter('rot_x', 'linear', 2, 3, False, val=5)
    ffd_block_1.add_shape_parameter('rot_x', 'quadratic', 4,5,False,val=2)
    ffd_block_2.add_shape_parameter('rot_y', 'quad', 4, 5, False)
    ffd_blocks = [ffd_block_1, ffd_block_2]

    # print(ffd_block_test.local_control_points)

    sim = Simulator(SectionPropertiesModel(ffd_blocks=ffd_blocks))
    sim.run()
    sim.prob.model.list_inputs(prom_name=True, print_arrays=True)
    sim.prob.model.list_outputs(prom_name=True, print_arrays=True)    
    # sim.visualize_implementation()


                
