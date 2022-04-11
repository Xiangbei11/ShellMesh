import numpy as np
from lsdo_kit.geometry.utils.generate_ffd import create_ffd
from lsdo_kit.design.design_geometry.core.ffd import FFD

class Component(object):
    def __init__(self, name: str, stp_entity_names: list, pointset = None, ffd_control_points=None, nxp=10, nyp=10, nzp=10, mass_density=1, thickness=1):
        # self.geo = geo

        self.design_variables_list = []        
        self.embedded_entity_names = []

        self.name = name     
        self.ffd_control_points = ffd_control_points
        self.stp_entity_names = stp_entity_names
        
        self.ffd_name = f'{self.name}_ffd'
        self.nxp = nxp
        self.nyp = nyp
        self.nzp = nzp

        self.pointset = pointset

        self.attribute_dict = {'mass_density': mass_density, 'thickness' : thickness}

    def add_design_variable(self, attribute):
        if attribute not in self.attribute_dict.keys():
            raise Exception('Invalied attribute!')
        else:
            self.design_variables_list.append(attribute)

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
        if self.ffd_control_points is None: 
            self._run_auto_rect_ffd()
        
        elif self.ffd_control_points is not None:
            self.ffd = FFD(self.ffd_name, self.ffd_control_points)

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
        self.embedded_entities_control_points = np.concatenate(self.embedded_entities_control_points, axis=0) # Concatenating all entity control points

   
    def _extract_min_max_values(self):
        
        self.x_min, self.y_min, self.z_min = self.embedded_entities_control_points.min(axis=0)
        self.x_max, self.y_max, self.z_max = self.embedded_entities_control_points.max(axis=0)


if __name__ == "__main__":
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    ''' Using rect_wing as an example for a Component '''
    path_name = '../geometry/examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)
    c = Component(stp_entity_names=['Fuse'], name='wing') # Creates a 
    f = FeatureComp(pointset=rib_pointset)
    # c._add_component_to_vehicle()


    c = Component()
    c.add_geo_info(stp_enti..)
    c.add_feature()




    geo.add_geometry_primitive(c)
    print(geo.components.keys())

    ''' ------------------- PLOT --------------------- ''' 
    print(c.embedded_entities_control_points.shape)
    print(c.ffd.control_points.shape)
    c.ffd.control_points = np.reshape(c.ffd.control_points, (c.nxp * c.nyp * c.nzp, 3))
    c.ffd.plot(10,10,10)

    vp_init = Plotter()
    vps1 = Points(c.embedded_entities_control_points, r=8, c = 'blue')
    vps2 = Points(c.ffd.control_points, r=8, c = 'red')
    vp_init.show([vps1, vps2], 'FFD', axes=1, viewup="z", interactive = True)



# import numpy as np

# class Feature(object):
#     def __init__(self, name, pointset, mass_density=1, thickness=1):
#         self.name = name 
#         self.pointset = pointset
#         self.design_variables_list = []

#         self.attribute_dict = {'mass_density': mass_density, 'thickness' : thickness}

#     def add_design_variable(self, attribute):
#         if attribute not in self.attribute_dict.keys():
#             raise Exception('Invalied attribute!')
#         else:
#             self.design_variables_list.append(attribute)

#         # self.mass_density = mass_density
#         # self.thickness = thickness


# if __name__ == "__main__":

#     '''
#     Challenges:
#         - We are going to add Features or Pointsets into the Mesh class DONE
#         - Feature just contains the information 
#     '''

#     # vlm_mesh = Mesh('vlm_mesh')
#     # my_feature = Feature('my_feature', mass_density = 5, thickness = 1, geometry = camber_surface_mesh)
#     # vlm_mesh.add_feature(camber_surface_mesh, name="camber_surface")

#     ''' 
#     TODO:
#         1. Create a method or technique to keep track of all features, should be a part of Mesh
#         2. Create the design variable creation in the mass properties model
#         3. Find a general way of computing area given arbitrary shapes
#         4. Clean up old code and test cases. Make tests directory and set them up for pytest. 
#     '''