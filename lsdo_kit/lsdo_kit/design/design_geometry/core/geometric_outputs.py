import numpy as np

#TODO: 1. Add subtraction capability to pointsets; 2. Add a custom compute method  ; 3. Compute Aspect Ratio method using two of the previous methods 


class GeometricOuputs(object):
    def __init__(self, geo):
        self.geo = geo
        self.geometric_outputs_dict = {}

    def compute_displacement(self, name, pointset1, pointset2):
        '''
        Objective: Compute displacement between two pointsets.

        Input: Two pointsets

        Output: Displacement between the two pointsets.
        '''
        # Does the order of subtraction matter? ptst1 - ptst2 or ptst2 - ptst1
        pointset3 = self.geo.subtract_pointsets(pointset1, pointset2)
        self.geometric_outputs_dict[name] = {'displacement' : pointset3}

    def compute_magnitude(self, name, pointset1, pointset2):
        '''
        Objective: Compute magnitude of distance between two pointsets.

        Input: Two pointsets

        Outout: Magntiude between two pointsets.
        '''
        pointset3 = self.geo.subtract_pointsets(pointset1, pointset2)
        self.geometric_outputs_dict[name] = {'magnitude' : pointset3} 

    def compute_angle(self, name, pointset1, pointset2):
        '''
        Objective: Compute angle between two pointsets (vectors).

        Input: Two pointsets.

        Output: Angle between the two pointsets.
        '''
        self.geometric_outputs_dict[name] = {'angle' : [pointset1, pointset2]}

    def compute_wireframe_area(self, name, wireframe):
      '''
      Objective: Compute up and side vectors to compute cross products.

      Input: Pointset

      Output: Up and side vectors for cross product computation
      '''
      # wireframe_shape = wireframe.shape 
      # print('Wireframe shape: ', wireframe.shape)
      # top_points = self.geo.extract_pointset(wireframe, point_indices=np.array(wireframe_shape[:,1:]), shape= np.array([wireframe_shape[0], wireframe_shape[1]-1]))
      # bottom_points = self.geo.extract_pointset(wireframe, point_indices=np.array(wireframe_shape[:,:-1]), shape= np.array([wireframe_shape[0], wireframe_shape[1]-1]))
      # right_points = self.geo.extract_pointset(wireframe, point_indices=np.array(wireframe_shape[1:,:]), shape= np.array([wireframe_shape[0]-1, wireframe_shape[1]]))
      # left_points = self.geo.extract_pointset(wireframe, point_indices=np.array(wireframe_shape[:-1,:]), shape= np.array([wireframe_shape[0]-1, wireframe_shape[1]]))

      # up_vectors = self.geo.subtract_pointsets(top_points,bottom_points)
      # side_vectors = self.geo.subtract_pointsets(right_points,left_points)
      # self.geometric_outputs_dict[name] = {'wireframe' : [up_vectors, side_vectors]}
      self.geometric_outputs_dict[name] = {'wireframe' : [wireframe]}

    def compute_quad_area(self, name, pointset00, pointset01, pointset11):
        '''
        Objective: Compute area through cross product.

        Input: Three pointsets to form two vectors.

        Output: Two vectors to take cross product of.
        '''

        pointset_bisect_1 = self.geo.subtract_pointsets(pointset00, pointset01)
        pointset_bisect_2 = self.geo.subtract_pointsets(pointset01, pointset11)
        
        self.geometric_outputs_dict[name] = {'quad_area' : [pointset_bisect_1, pointset_bisect_2]}

    def compute_elliptical_area(self, name, semi_minor, semi_major):
        '''
        Objective: Compute area of an ellipse.

        Input:

        Output: Area of ellipse.
        '''
        self.geometric_outputs_dict[name] = {'ellipse_area' : [semi_minor, semi_major]}

    def compute_aspect_ratio(self, name, wireframe, pointset1, pointset2):
        self.compute_wireframe_area(name + 'AR_area', wireframe)
        self.compute_magnitude(name + 'AR_span', pointset1, pointset2)
        self.geometric_outputs_dict[name] = {'AR' : [name + 'AR_area', name + 'AR_span']}



if __name__ == "__main__":
    import csdl
    import numpy as np
    from lsdo_kit.design.design_geometry.utils.generate_ffd import create_ffd
    from lsdo_kit.design.design_geometry.core.ffd import FFD
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    from lsdo_kit.old_files.mesh import Mesh
    from lsdo_kit.design.design_geometry.core.component import Component
    from csdl_om import Simulator

    import matplotlib.pyplot as plt
    from vedo import Points, Plotter, LegendBox

    path_name = '../../examples/CAD/'
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

    wing_comp = Component(name='wing', stp_entity_names=['RectWing'])
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


    num_pts1 = [20]#20
    num_pts2 = [20]#5
    surface_shape1 = np.append(num_pts1,num_pts2)
    lead_edge_curve = geo.perform_linear_interpolation(wing_lead_left, wing_lead_right, num_pts1)
    trail_edge_curve = geo.perform_linear_interpolation(wing_trail_left, wing_trail_right, num_pts1)


    chord_surface = geo.perform_linear_interpolation(lead_edge_curve, trail_edge_curve, surface_shape1)
    geo.assemble(chord_surface)
    chord_surface_points = geo.evaluate(chord_surface)

    '''Project from chord to wing surfaces'''
    top_points, top_plot = geo.project_points(chord_surface_points, projection_direction=up_direction, projection_targets_names=top_wing_surface_names)
    bot_points, bot_plot = geo.project_points(chord_surface_points, projection_direction=down_direction, projection_targets_names=bot_wing_surface_names)

    '''Create Camber Surface'''
    surface_shape2 = np.append(surface_shape1,1)
    output_parameters = 0.5*np.ones([num_pts1[0]*num_pts2[0],1])
    camber_surface_mesh = geo.perform_linear_interpolation(top_points, bot_points, surface_shape2, output_parameters = output_parameters)

    vlm_mesh = Mesh('vlm_mesh')
    vlm_mesh.add_pointset(camber_surface_mesh, name="camber_surface")
    # vlm_mesh.add_pointset(chord_surface, name='chord_surface')

    # vlm_mesh = Mesh('vlm_mesh')
    # my_feature = Feature('my_feature', mass_density = 5, thickness = 1, geometry = camber_surface_mesh)
    # vlm_mesh.add_feature(camber_surface_mesh, name="camber_surface")

    chord_surface_mesh = Mesh('chord_surface_mesh')
    chord_surface_mesh.add_pointset(chord_surface, name='chord_surface')

    meshes = [vlm_mesh, chord_surface_mesh]

    # num_pts1 = [20]#20
    # num_pts2 = [20]#5
    # surface_shape1 = np.append(num_pts1,num_pts2)
    # lead_edge_curve = geo.perform_linear_interpolation(wing_lead_left, wing_lead_right, num_pts1)
    # trail_edge_curve = geo.perform_linear_interpolation(wing_trail_left, wing_trail_right, num_pts1)

    geo_outputs = GeometricOuputs(geo=geo)
    geo_outputs.compute_magnitude('test_mag', wing_lead_left, wing_trail_left)
    # geo_outputs.compute_wireframe_area('test_wire', camber_surface_mesh)
    print(geo_outputs.geometric_outputs_dict)

    # subtraction = geo.subtract_pointsets(wing_lead_left, wing_trail_left)
    # print('lead point shape:', wing_lead_left.shape)
    # print('trail point shape:', wing_trail_left.shape)
    
    
    # geo.assemble()
    # geo.evaluate()

    # print('lead point', wing_lead_left.physical_coordinates)
    # print('trail point', wing_trail_left.physical_coordinates)
    # print('subtraction:', subtraction.physical_coordinates)
    # print('sanity: ', wing_lead_left.physical_coordinates - wing_trail_left.physical_coordinates)