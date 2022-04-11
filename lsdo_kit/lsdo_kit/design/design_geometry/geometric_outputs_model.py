from tracemalloc import start
import csdl
from lsdo_kit.design.design_geometry.evaluate_pointsets import EvaluatePointsets
from lsdo_kit.design.design_geometry.wireframe_area_functions import *
import numpy as np

class GeometricOutputsModel(csdl.Model):

    def initialize(self):
        self.parameters.declare('design_geometry_obj')

    def define(self):
        design_geometry_obj = self.parameters['design_geometry_obj']
        geometric_outputs = design_geometry_obj.geometric_outputs

        geometric_outputs_dict = geometric_outputs.geometric_outputs_dict
        points = self.declare_variable('points', shape=(design_geometry_obj.eval_map.shape[0],3))

        print('Points: ', points)
        print('Points shape: ', points.shape)
        for operation_name, operation_dict in geometric_outputs_dict.items():
            if 'displacement' in operation_dict.keys():
                pointset = list(operation_dict.values())
                ind = pointset[0].output_starting_ind
                pt_shape = np.cumprod(pointset[0].shape)[-2]
                print('IND: ', ind)
                print('IND + shape: ', ind+pt_shape)
                output_pointset = points[ind : ind+pt_shape, :]
                self.register_output(operation_name, output_pointset)

            elif 'magnitude' in operation_dict.keys():
                pointset = list(operation_dict.values())
                ind = pointset[0].output_starting_ind
                pt_shape = np.cumprod(pointset[0].shape)[-2]
                print('ind: ', ind)
                print('Pointset value: ', pointset[0].physical_coordinates)
                print('Original Pointset Shape: ', pointset[0].shape)
                print('pt_shape: ', pt_shape)
                output_pointset = csdl.pnorm(points[ind : ind+pt_shape, :])
                self.register_output(operation_name, output_pointset)

            elif 'angle' in operation_dict.keys():
                pointsets = list(operation_dict.values())
                pointsets = pointsets.pop()
                # print('pointsets in angle: ', pointsets)
                design_geometry_obj.evaluate()
                print('INMODEL pointset1: ', pointsets[0].physical_coordinates)
                print('INMODEL pointset2: ', pointsets[1].physical_coordinates)
                pointset1 = pointsets[0]
                pointset2 = pointsets[1]
                # print('pointset1 name: ', pointset1.name)
                # print('pointset2 name: ', pointset2.name)
                ind1 = pointset1.output_starting_ind
                ind2 = pointset2.output_starting_ind
                pt_shape1 = np.cumprod(pointset1.shape)[-2]
                pt_shape2 = np.cumprod(pointset2.shape)[-2]
                
                num_dot = csdl.dot(points[ind1 : ind1 + pt_shape1, :], points[ind2 : ind2 + pt_shape2, :], axis=1)
                den_prod = csdl.pnorm(points[ind1 : ind1 + pt_shape1, :]) * csdl.pnorm(points[ind2 : ind2 + pt_shape2, :])                
                output_pointset = num_dot / den_prod
                angle = csdl.arccos(output_pointset)

                self.register_output('points1', points[ind1 : ind1 + pt_shape1, :])
                self.register_output('points2', points[ind2 : ind2 + pt_shape2, :])
                self.register_output('num_dot', num_dot)
                self.register_output('den_prod', den_prod)
                self.register_output('division', output_pointset)
                self.register_output(operation_name, angle)

            elif 'quad_area' in operation_dict.keys():
                '''
                Objective: Compute cross product of two vectors and 2-norm to find area between specified points.

                Output: Area between the two vectors/pointsets.

                '''
                pointsets = list(operation_dict.values())
                pointset1 = pointsets[0][0]
                pointset2 = pointsets[0][1] 
                    
                ind1 = int(pointset1.output_starting_ind)
                ind2 = int(pointset2.output_starting_ind)

                pt_shape1 = np.cumprod(pointset1.shape)[-2]
                pt_shape2 = np.cumprod(pointset2.shape)[-2]

                area = csdl.cross(points[ind1 : ind1 + pt_shape1, :], points[ind2 : ind2 + pt_shape2, :],axis = 1)
                area = csdl.pnorm(area)

                self.register_output(operation_name, area)
                     
            elif 'wireframe' in operation_dict.keys():
                pointsets = list(operation_dict.values())
                wireframe_shape = pointsets[0][0].shape                                                         # shape of mesh
                n = wireframe_shape[0]                                                                          # number of rows in pointset subset
                m = wireframe_shape[1]                                                                          # number of columns in pointset subset
                ind1 = pointsets[0][0].output_starting_ind                                                      # starting index of this pointset subset in parent pointset
                pointset_vector = design_geometry_obj.eval_map.dot(design_geometry_obj.total_cntrl_pts_vector)[ind1:ind1+(n*m),:]     # extract pointset subset from parent subset
                mesh_coordinates = mesh_coordinates_list(n,m,pointset_vector)                                   # extract mesh vertices
                up_vectors, side_vectors = get_vectors(mesh_coordinates)
                area_sum = self.declare_variable('area_sum', val = 0)
                cross_list1 = np.cross(up_vectors, side_vectors)

                for i in range(cross_list1.shape[0]):
                    value = self.declare_variable(f'value{i}', val = cross_list1[i])
                    area_sum = area_sum + csdl.pnorm(value)
                
                # up_vectors = pointsets[0]
                # side_vectors = pointsets[1]
                # ind1 = up_vectors.output_starting_ind
                # ind2 = side_vectors.output_starting_ind
                # pt_shape1 = np.cumprod(up_vectors.shape)[-2]
                # pt_shape2 = np.cumprod(side_vectors.shape)[-2]

                # area_vec = csdl.cross(points[ind1:ind1+pt_shape1], points[ind2:ind2+pt_shape], axis=2)
                # area_mag = csdl.pnorm(area_vec)
                # area_sum = csdl.sum(area_mag)
                # area_sum = 0
                # i = 1
                # for num in pt_shape1:
                #     # area = csdl.cross(up_vectors, side_vectors, axis=2)
                #     area = csdl.cross(up_vectors[ind1 : ind1 + i, :], side_vectors[ind2 : ind2 + i, :])
                #     area = csdl.pnorm(area)
                #     ind1 = ind1 + i; ind2 = ind2 + i
                #     area_sum = area_sum + area
                self.register_output('#', points[0,:])   
                self.register_output(operation_name, area_sum) 

            elif 'ellipse_area' in operation_dict.keys():
                pointsets = list(operation_dict.values())
                pointset1 = pointsets[0]
                pointset2 = pointsets[1]

                ind1 = pointset1.output_starting_ind
                ind2 = pointset2.output_starting_ind

                pt_shape1 = np.cumprod(pointset1.shape)[-2]
                pt_shape2 = np.cumprod(pointset2.shape)[-2]

                semi_minor = self.create_input(operation_name + 'minor', points[ind1 : ind1 + pt_shape1, :])
                semi_major = self.create_input(operation_name + 'major', points[ind2 : ind2 + pt_shape2, :])
                semi_minor_axis_length = csdl.pnorm(semi_minor)
                semi_major_axis_length = csdl.pnorm(semi_major)
                ellipse_area = np.pi * semi_major_axis_length/2 * semi_minor_axis_length/2
                self.register_output(operation_name, ellipse_area)

            # elif 'AR' in operation_dict.keys():


if __name__ == "__main__":
    import csdl
    import numpy as np
    from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    from lsdo_kit.old_files.mesh import Mesh
    from lsdo_kit.design.design_geometry.core.component import Component
    from lsdo_kit.design.design_geometry.core.geometric_outputs import GeometricOuputs

    from csdl_om import Simulator

    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    path_name = '../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    top_wing_surface_names = [
        'RectWing, 0, 3', 
        'RectWing, 1, 9', 
        ]

    bot_wing_surface_names = [
        'RectWing, 0, 2',
        'RectWing, 1, 8', 
        ]

    wing_comp = Component(stp_entity_names=['RectWing'], name='wing')
    # c._add_GeometryPrimitve_to_vehicle()
    geo.add_component(wing_comp)

    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    left_lead_point = np.array([0., -9000., 2000.])/1000
    left_trail_point = np.array([4000.0, -9000.0, 2000.])/1000
    right_lead_point = np.array([0.0, 9000.0, 2000.])/1000
    right_trail_point = np.array([4000.0, 9000.0, 2000.])/1000

    '''Project points'''
    wing_lead_left, wing_lead_left_coord = geo.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_left ,wing_trail_left_coord = geo.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_lead_right ,wing_lead_right_coord = geo.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_right ,wing_trail_right_coord = geo.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

    geo_outputs = GeometricOuputs(geo=geo)
    geo_outputs.compute_magnitude('test_mag', wing_lead_left, wing_trail_left)

    sim = Simulator(GeometricOutputsModel(geometric_outputs=geo_outputs))
    sim.run()
    sim.visualize_implementation()

    # GeometricOutputsModel(geometric_outputs=geo_outputs)
    # print(geo_outputs.geometric_outputs_dict)