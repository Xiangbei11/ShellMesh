import numpy as np

class Feature(object):
    
    def __init__(self, name, pointset, parents=None, mass_density=1, thickness=1):
        self.name = name 
        self.pointset = pointset
        self.design_variables_list = []
        self.parents = parents

        # self.attribute_dict = {'mass_density': mass_density, 'thickness' : thickness}

    def add_mass(self, mass):
        pass

    def add_thickness(self, thickness):
        pass

    # def add_design_variable(self, attribute):
    #     if attribute not in self.attribute_dict.keys():
    #         raise Exception('Invalied attribute!')
    #     else:
    #         self.design_variables_list.append(attribute)

        # self.mass_density = mass_density
        # self.thickness = thickness


if __name__ == "__main__":

    '''
    Challenges:
        - We are going to add Features or Pointsets into the Mesh class DONE
        - Feature just contains the information 
    '''

    # vlm_mesh = Mesh('vlm_mesh')
    # my_feature = Feature('my_feature', mass_density = 5, thickness = 1, geometry = camber_surface_mesh)
    # vlm_mesh.add_feature(camber_surface_mesh, name="camber_surface")

    ''' 
    TODO:
        1. Create a method or technique to keep track of all features, should be a part of Mesh
        2. Create the design variable creation in the mass properties model
        3. Find a general way of computing area given arbitrary shapes
        4. Clean up old code and test cases. Make tests directory and set them up for pytest. 


    '''