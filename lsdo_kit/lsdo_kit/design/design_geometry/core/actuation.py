from lsdo_kit.cython.get_open_uniform_py import get_open_uniform
import numpy as np
import matplotlib.pyplot as plt
from vedo import Points, Plotter, LegendBox

from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry
from lsdo_kit.design.design_geometry.core.pointset import PointSet

# from lsdo_kit.geometry.utils.calculate_rotation_mat import calculate_rotation_mat


# class Joint(object):
#     # TODO Generalize to make sense for design joints and kinematic joints as well
#     def __init__(self, origin, axis=None, components=None, affected_volume=None):
#         self.origin = origin
#         self.axis = axis
#         self.components = components

#         self._find_points_in_volume()


#     # TODO need for ailerons, etc.
#     def _find_points_in_volume():
#         # If numpy array, use them to generate BSplineVolume
#         # If FFD, extract the BSplineVolume
#         # If BSplineVolume, great
#         # 
#         # If child_entities = None, set child_entities to the whole geometry
#         # Project the points of the child entities into this volume to find which points are affected
#         # Save the indices of the affected points as the child_indices_vec
#         pass

        
# class Actuation(Joint):
class Actuation:
    # def __init__(self, rot_or_trans, origin, axis=None, components=None, affected_volume=None, actuation_profile=None):
    def __init__(self, name, actuation_profile, origin, pointset2, actuating_components, rot_or_trans='rot', affected_volume=None):
        
        # axis = pointset2-origin # TODO this probably needs to be done with the design_geometry method.
        self.name = name
        self.actuation_profile = actuation_profile
        self.origin = origin
        self.pointset2 = pointset2
        # self.axis=axis 
        self.actuating_components=actuating_components
        self.affected_volume=affected_volume
        self.rot_or_trans=rot_or_trans
        # self.rot_or_trans = rot_or_trans
        # self.actuation_profile = actuation_profile

        # get parent dependance from points that origin is dependent on
        # self.get_dependance()
    
    # Gets indices of control points that it is dependent on
    def get_dependance(self):

        if self.origin is type(np.ndarray):
            self.parent_indices_vec = np.array([])
            return

        if self.origin is not type(PointSet):
            Exception('Please pass in a pointset or a numpy array for the joint origin.')

        # geo.evaluate_absolute_map(self.origin)
        map = self.origin.absolute_map  #TODO will the map be evaluated at this step?
        _, self.parent_indices_vec = map.nonzero()


if __name__ == "__main__":
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry
    from lsdo_kit.design.design_geometry.core.pointset import PointSet
    from lsdo_kit.design.design import Design
    from lsdo_kit.design.design_geometry.core.component import Component
    from lsdo_kit.simulation.simulation import Simulation   

    import numpy as np

    actuation_profile = np.linspace(0,1,100)
    name = 'test_actuation'

    path_name = '../examples/CAD/'
    file_name = 'rect_wing.stp'

    design = Design(path_name + file_name)
    geo = design.design_geometry

    ''' Creating a component and adding to design '''
    c = Component(stp_entity_names=['RectWing'], name='wing')
    design.add_component(c)


    right_lead_point = np.array([0., 9., 3.])
    right_trail_point = np.array([4., 9., 3.])
    left_lead_point = np.array([0., -9., 3.])
    left_trail_point = np.array([4., -9., 3.])

    ''' Defining Points and directions for projections '''
    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    ''' Project points '''
    left_lead, wing_lead_left_coord = design.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    left_trail ,wing_trail_left_coord = design.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

    # print('left lead: ', wing_lead_left_coord)
    # print('left trail: ', wing_trail_left_coord)

    right_lead ,wing_lead_right_coord = design.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    right_trail ,wing_trail_right_coord = design.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

    # print('right lead: ', wing_lead_right_coord)
    # print('right trail: ', wing_trail_right_coord)
    
    output_parameters = 0.25*np.ones([3,1])
    quarter_left = design.perform_linear_interpolation(pointset_start=left_lead, pointset_end=left_trail, shape=(1,), output_parameters=[0.25])
    quarter_right = design.perform_linear_interpolation(pointset_start=right_lead, pointset_end=right_trail, shape=(1,), output_parameters=[0.25])
    
    geo.assemble()

    actuation_obj = Actuation(name='wing_rot', actuation_profile=actuation_profile, origin=quarter_left, pointset2=quarter_right, actuating_components=[c])




    # geo.assemble()
    # geo.evaluate()

    # print('POST INTERPOLATION')
    # print('Left  quarter chord coord: ', quarter_left.physical_coordinates)
    # print('Right quarter chord coord: ', quarter_right.physical_coordinates)
    