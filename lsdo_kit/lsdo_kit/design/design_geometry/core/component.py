import numpy as np
from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD

class Component(object):
    def __init__(self, name: str, stp_entity_names: list, ffd_obj=None, nxp=10, nyp=10, nzp=10):
        '''
        This class sets up FFD blocks for a specified part of the step file.
        
        '''
        # self.geo = geo
        self.ffd_obj = ffd_obj
        self.embedded_entity_names = []
        self.stp_entity_names = stp_entity_names
        self.name = name
        self.ffd_name = f'{self.name}_ffd'
        self.nxp = nxp
        self.nyp = nyp
        self.nzp = nzp
        self.shape_parameter_list = []

        # if self.ffd_obj is None: 
        #     self._extract_embedded_entity_data(self)
        #     self._auto_generate_ffd()
        
        # else:
        #     self.define_ffd(ffd_obj)
    def add_shape_parameter(self, property_name: str, parameter_name: str, order: int , num_cp: int, dv: str, val = None):
        if self.ffd_obj is None:
            self.shape_parameter_list.append([property_name, parameter_name, order, num_cp, dv, val])
        else:
            self.ffd.add_shape_parameter(self, property_name=property_name, parameter_name=parameter_name, order=order , num_cp=num_cp, dv=dv, val = val)

    def set_thickness(self, thickness):
        self.thickness = thickness

    def define_ffd(self, ffd_obj):

        self.ffd = ffd_obj
    #     ''' 
    #     Allow the user to pass in both the value and name, but is the name is passed in then the name 
    #     takes priority over the value
    #     '''

    #     if mass is not None and mass_name is not None:
    #         Raise Exception('Cannot feed in both mass, and mass_name')

    #     self.mass = mass
    #     self.mass_name = mass_name

    #     if cg is not None and cg_name is not None:
    #         Raise Exception('Cannot feed in both cg, and cg_name')

    #     self.cg = cg
    #     self.cg_name = cg_name


    #     if MoI is not None and MoI_name is not None:
    #         Raise Exception('Cannot feed in both MoI, and MoI_name')

    #     self.MoI = MoI
    #     self.MoI_name = MoI_name

    
    def _auto_generate_ffd(self):
        if self.ffd_obj is None: 
            self._run_auto_rect_ffd()
        
        elif self.ffd_obj is not None:
            self.ffd = self.ffd_obj                                                                                                                                                                                                                                                                                 

        self.ffd.add_embedded_entities(self.embedded_entities_objects)

    def _run_auto_cylinder_ffd(self):
        # Extracts the min/max values of the entities along all axes
        self._extract_min_max_values()

        # Auto creates the ffd block for the component
        self._auto_create_cylinder_ffd_block()

    def _run_auto_rect_ffd(self):
        # Extracts the min/max values of the entities along all axes
        self._extract_min_max_values()

        # Auto creates the ffd block for the component
        self._auto_create_rect_ffd_block()


    def _auto_create_cylinder_ffd_block(self):
        pass

    def _auto_create_rect_ffd_block(self):
        ''' 
        Take the maximum and minimum values for x y z that were computed, and create the ffd block
        '''

        ''' 
        Dont know what these values should be. Maybe we can look at the lengths of the different dimensions 
        and come up with a algorithm that determines how many sections there are based on dimension length.
        '''

        x_buffer = (self.x_max - self.x_min) * 0.01
        y_buffer = (self.y_max - self.y_min) * 0.01
        z_buffer = (self.z_max - self.z_min) * 0.01

        # Arbitrary buffers set for FFD block to make sure points are encapsulated
        mins = np.array([self.x_min - x_buffer, self.y_min - y_buffer, self.z_min - z_buffer])
        maxs = np.array([self.x_max + x_buffer, self.y_max + y_buffer, self.z_max + z_buffer]) 

        pt000 = np.array([mins[0], mins[1], mins[2]])
        pt001 = np.array([mins[0], mins[1], maxs[2]])

        pt010 = np.array([mins[0], maxs[1], mins[2]])
        pt011 = np.array([mins[0], maxs[1], maxs[2]])

        pt100 = np.array([maxs[0], mins[1], mins[2]])
        pt101 = np.array([maxs[0], mins[1], maxs[2]])
        
        pt110 = np.array([maxs[0], maxs[1], mins[2]])
        pt111 = np.array([maxs[0], maxs[1], maxs[2]])

        control_points = np.zeros((2,2,2,3))

        control_points[0,0,0,:] = pt000
        control_points[0,0,1,:] = pt001

        control_points[0,1,0,:] = pt010
        control_points[0,1,1,:] = pt011

        control_points[1,0,0,:] = pt100
        control_points[1,0,1,:] = pt101

        control_points[1,1,0,:] = pt110
        control_points[1,1,1,:] = pt111

        self.ffd_control_points = create_ffd(control_points, self.nxp, self.nyp, self.nzp)
        self.ffd = FFD(self.ffd_name, self.ffd_control_points)

        for param_list in self.shape_parameter_list:
            self.ffd.add_shape_parameter(property_name=param_list[0], parameter_name=param_list[1], order=param_list[2], 
                num_cp=param_list[3], dv=param_list[4], val=param_list[5])

    
        
    def _extract_embedded_entity_data(self, geo):
        ''' 
        This method is meant to extract the control points from the entity specified, and find 
        the min/max values in the x,y,z directions in order to make a rectangular prism FFD block
        around the specified entities.

        THIS METHOD REQUIRES THE USER TO JUST FEED IN THE NAME OF THE ENTITY. FOR EXAMPLE 'RECTWING'
        AND NOTHING MORE. IT IS LEAVING UP TO THE USER TO NAME ALL THE ENTITIES DIFFERENTLY, i.e NOT 
        HAVING MULTIPLE ENTITIES NAMED THE SAME THING. 
        '''

        self.embedded_entities_objects= [] # We want to compile a list of all the objects that this component contains
        self.embedded_entities_control_points = [] # We want to compile a list of all the control points that this component represents
        bspline_keys_list = list(geo.input_bspline_entity_dict.keys())
                
        ''' These nested for loops count the number of entities that have that name '''
        for entity_name in self.stp_entity_names:    
            for key in bspline_keys_list:
                if entity_name in key:
                    self.embedded_entity_names.append(key)
                    self.embedded_entities_objects.append(geo.input_bspline_entity_dict[key])
                    self.embedded_entities_control_points.append(geo.input_bspline_entity_dict[key].control_points)

        # print("OBJECT ATTRIBUTE: ", self.embedded_entities_objects[-1].starting_geometry_index)
        self.embedded_entities_control_points = np.concatenate(self.embedded_entities_control_points, axis=0) # Concatenating all entity control points

   
    def _extract_min_max_values(self):
        
        self.x_min, self.y_min, self.z_min = self.embedded_entities_control_points.min(axis=0)
        self.x_max, self.y_max, self.z_max = self.embedded_entities_control_points.max(axis=0)


if __name__ == "__main__":
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox
    import os
    os.chdir("/Users/myronphan/lsdo_geo/lsdo_kit/design/design_geometry/core")

    ''' Using rect_wing as an example for a Component '''
    path_name = '../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)
    c = Component(stp_entity_names=['RectWing'], name='wing')
    # c._add_component_to_vehicle()
    geo.add_component(c)
    # print(geo.components.keys())

    ''' ------------------- PLOT --------------------- ''' 
    # print(c.embedded_entities_control_points.shape)
    # print(c.ffd.control_points.shape)
    c.ffd.control_points = np.reshape(c.ffd.control_points, (c.nxp * c.nyp * c.nzp, 3))
    c.ffd.plot(10,10,10)

    vp_init = Plotter()
    vps1 = Points(c.embedded_entities_control_points, r=8, c = 'blue')
    vps2 = Points(c.ffd.control_points, r=8, c = 'red')
    vp_init.show([vps1, vps2], 'FFD', axes=1, viewup="z", interactive = True)



