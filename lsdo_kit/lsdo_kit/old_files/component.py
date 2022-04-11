import numpy as np
from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD

class Component(object):
    def __init__(self, geo, stp_entity_names: list , alias: str, ffd_control_points=None):
        self.geo = geo
        self.stp_entity_names = stp_entity_names
        self.alias = alias
        self.ffd_name = f'{self.alias}_ffd'

    def set_properties(self, mass=None, mass_name=None, cg=None, cg_name=None, MoI=None, MoI_name=None):
        ''' 
        Allow the user to pass in both the value and name, but is the name is passed in then the name 
        takes priority over the value
        '''

        if mass is not None and mass_name is not None:
            raise Exception('Cannot feed in both mass, and mass_name')

        self.mass = mass
        self.mass_name = mass_name

        if cg is not None and cg_name is not None:
            raise Exception('Cannot feed in both cg, and cg_name')

        self.cg = cg
        self.cg_name = cg_name


        if MoI is not None and MoI_name is not None:
            raise Exception('Cannot feed in both MoI, and MoI_name')

        self.MoI = MoI
        self.MoI_name = MoI_name

    def _add_component_to_vehicle(self):
        # Computes the embedded entites control points
        self._extract_embedded_entity_control_points(self.geo, self.stp_entity_names)
        
        # Extracts the min/max values of the entities along all axes
        self._extract_min_max_values()

        # Auto creates the ffd block for the component
        self._auto_create_rect_ffd_block()


        # if ffd_control_points is None:
        #     self.define_ffd(f'{self.alias}_ffd', )



    def _define_ffd(self, name, control_points, embedded_entities, local_axes=dict(xprime=None, yprime=None, zprime=None),
                    origin=None, local_control_points=None, shape=None, order_u=4, order_v=4, order_w=4,
                    knots_u=None, knots_v=None, knots_w=None):

        self.ffd = FFD(self.ffd_name, self.ffd_control_points)

    def _auto_create_rect_ffd_block(self):
        ''' 
        Take the maximum and minimum values for x y z that were computed, and create the ffd block
        '''

        ''' 
        Dont know what these values should be. Maybe we can look at the lengths of the different dimensions 
        and come up with a algorithm that determines how many sections there are based on dimension length.
        '''
        self.nxp = 10
        self.nyp = 10
        self.nzp = 10
        buffer = 1

        # Arbitrary buffers set for FFD block to make sure points are encapsulated
        mins = np.array([self.x_min, self.y_min, self.z_min]) - buffer 
        maxs = np.array([self.x_max, self.y_max, self.z_max]) + buffer 

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
        
    def _extract_embedded_entity_control_points(self, geo, stp_entity_names):
        ''' 
        This method is meant to extract the control points from the entity specified, and find 
        the min/max values in the x,y,z directions in order to make a rectangular prism FFD block
        around the specified entities.

        THIS METHOD REQUIRES THE USER TO JUST FEED IN THE NAME OF THE ENTITY. FOR EXAMPLE 'RECTWING'
        AND NOTHING MORE. IT IS LEAVING UP TO THE USER TO NAME ALL THE ENTITIES DIFFERENTLY, i.e NOT 
        HAVING MULTIPLE ENTITIES NAMED THE SAME THING. 
        '''

        self.embedded_entities_control_points = [] # We want to compile a list of all the control points that this component represents
        bspline_keys_list = list(geo.input_bspline_entity_dict.keys())
                
        ''' These nested for loops count the number of entities that have that name '''
        for entity_name in stp_entity_names:    
            for key in bspline_keys_list:
                if entity_name in key:
                    self.embedded_entities_control_points.append(geo.input_bspline_entity_dict[key].control_points)
                    
        self.embedded_entities_control_points = np.concatenate(self.embedded_entities_control_points, axis=0) # Concatenating all entity control points 
   
    def _extract_min_max_values(self):
        
        self.x_min, self.y_min, self.z_min = self.embedded_entities_control_points.min(axis=0)
        self.x_max, self.y_max, self.z_max = self.embedded_entities_control_points.max(axis=0)


if __name__ == "__main__":
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    ''' Using rect_wing as an example for a Component '''
    path_name = '../../../lsdo_geo/lsdo_geo/examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)
    c = Component(geo=geo, stp_entity_names=['RectWing'], alias='wing')
    c._add_component_to_vehicle()


    ''' ------------------- PLOT --------------------- ''' 
    print(c.embedded_entities_control_points.shape)
    print(c.ffd.control_points.shape)
    c.ffd.control_points = np.reshape(c.ffd.control_points, (c.nxp * c.nyp * c.nzp, 3))

    vp_init = Plotter()
    vps1 = Points(c.embedded_entities_control_points, r=8, c = 'blue')
    vps2 = Points(c.ffd.control_points, r=8, c = 'red')
    vp_init.show([vps1, vps2], 'FFD', axes=1, viewup="z", interactive = True)



