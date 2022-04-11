from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry
import numpy as np
# from lsdo_kit.design.design_geometry.design_model import DesignModel

class Design(object):
    def __init__(self, geometry_stp_file):
        self.component_dict = {}
        self.feature_dict = {}
        self.design_geometry = DesignGeometry(geometry_stp_file)
        
    def add_component(self, component):
        self.component_dict[component.name] = component
        self.design_geometry.add_component(component)

    def add_feature(self, feature):
        self.feature_dict[feature.name] = feature

    def add_geometric_outputs(self, geo_outputs):
        self.design_geometry.add_geometric_outputs(geo_outputs=geo_outputs)

    def read_file(self, file_name):
        self.design_geometry.read_file(file_name=file_name)

    def read_openvsp_stp(self, file_name):
        self.design_geometry.read_openvsp_stp(file_name=file_name)

    def read_iges(self,file_name):
        self.design_geometry.read_iges(file_name=file_name)

    def write_step(self, file_name, plot=False):
        self.design_geometry.write_step(file_name=file_name, plot=plot)

    def write_iges(self, file_name, plot = False):
        self.design_geometry.write_iges(file_name=file_name, plot=plot)

    def generate_bspline_list(self, names_list):
        self.design_geometry.generate_bspline_list(names_list=names_list)

    def project_points(self, points_to_be_projected, projection_targets_names=[], projection_direction=np.array([0., 0., 0.]), offset=np.array([0., 0., 0.]), plot=False):
        return self.design_geometry.project_points(points_to_be_projected=points_to_be_projected, projection_targets_names=projection_targets_names, projection_direction=projection_direction, offset=offset, plot=plot)

    def perform_linear_combination(self, parent_pointset_list, coefficients_matrix, shape, offset=np.array([0., 0., 0.])):
        return self.design_geometry.perform_linear_combination(parent_pointset_list=parent_pointset_list, coefficients_matrix=coefficients_matrix, shape=shape, offset=offset)

    def _define_linear_combination(self, parent_pointset_list, relative_map, shape, offset=np.array([0., 0., 0.])):
        return self.design_geometry._define_linear_combination(parent_pointset_list=parent_pointset_list, relative_map=relative_map, shape=shape, offset=offset)

    def perform_linear_interpolation(self, pointset_start, pointset_end, shape, output_parameters=np.array([0]), offset=np.array([0., 0., 0.])):       
        return self.design_geometry.perform_linear_interpolation(pointset_start=pointset_start, pointset_end=pointset_end, shape=shape, output_parameters=output_parameters, offset=offset)

    def perform_bilinear_interpolation(self, point_00, point_10, point_01, point_11, shape, output_parameters=np.array([0., 0., 0.]), offset=np.array([0., 0., 0.])):
        return self.design_geometry.perform_bilinear_interpolation(point_00=point_00, point_10=point_10, point_01=point_01, point_11=point_11, shape=shape, output_parameters=output_parameters, offset=offset)

    def perform_trilinear_interpolation(self, point_000, point_100, point_010, point_110, point_001, point_101, point_011, point_111, 
            shape, output_parameters=np.array([0., 0., 0.]), offset=np.array([0., 0., 0.])):  

        return self.design_geometry.perform_trilinear_interpolation(point_000=point_000, point_100=point_100, point_010=point_010, point_110=point_110, point_001=point_001, point_101=point_101, point_011=point_011, point_111=point_111,
                shape=shape, output_parameters=output_parameters, offset=offset)

    def perform_2d_transfinite_interpolation(self, u_curve0, u_curve1, v_curve0, v_curve1, output_parameters=np.array([0., 0., 0.]), offset=np.array([0., 0., 0.])):
        return self.design_geometry.perform_2d_transfinite_interpolation(u_curve0=u_curve0, u_curve1=u_curve1, v_curve0=v_curve0, v_curve1=v_curve1, output_parameters=output_parameters, offset=offset)

    def perform_3d_transfinite_interpolation(self, output_parameters=None, offset=np.array([0., 0., 0.])):
        return self.design_geometry.perform_3d_transfinite_interpolation(output_parameters=output_parameters, offset=offset)

    def extract_pointset(self, parent_pointset, point_indices, shape, offset=np.array([0., 0., 0.])):
        return self.design_geometry.extract_pointset(parent_pointset=parent_pointset, point_indices=point_indices, shape=shape, offset=offset)

    def add_pointsets(self, pointset1, pointset2, offset=np.array([0., 0., 0.])):
        return self.design_geometry.add_pointsets(pointset1=pointset1, pointset2=pointset2, offset=offset)

    def subtract_pointsets(self, pointset1, pointset2, offset=np.array([0., 0., 0.])):
        return self.design_geometry.subtract_pointsets(pointset1=pointset1, pointset2=pointset2, offset=offset)

    def register_output(self, pointset, name=None, mesh_names=[], bspline_mesh_names=[]):
        self.design_geometry.register_output(pointset=pointset, name=name, mesh_names=mesh_names, bspline_mesh_names=bspline_mesh_names)

    def assemble_all_absolute_maps(self):
        self.design_geometry.assemble_all_absolute_maps()

    def assemble_my_absolute_map(self):
        self.design_geometry.assemble_my_absolute_map()

    def assemble_absolute_map(self, point_set):
        self.design_geometry.assemble_absolute_map(point_set=point_set)

    def assemble(self, pointset = None):
        self.design_geometry.assemble(pointset=pointset)

    def evaluate(self, pointset=None):
        self.design_geometry.evaluate(pointset=pointset)

    def fit_bspline_entities(self, pointset_list, output_vec=None):
        self.design_geometry.fit_bspline_entities(pointset_list=pointset_list, output_vec=output_vec)

    def remove_multiplicity(self, bspline_entity):
        self.design_geometry.remove_multiplicity(bspline_entity=bspline_entity)


if __name__ == "__main__": 
    path_name = 'examples/CAD/'
    file_name = 'rect_wing.stp'

    design = Design(path_name + file_name)

    from lsdo_kit.design.design_geometry.tests.organized_test.ffd_block_creation_script import generate_ffd_blocks
    ffd_blocks = generate_ffd_blocks(design.design_geometry)

    ''' Creating the chamber line mesh for mesh eval '''
    from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes
    meshes, camber_surface_mesh = generate_meshes(design.design_geometry)

    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    left_lead_point = np.array([0.0, -9000., 2000.])/1000
    left_trail_point = np.array([4000.0, -9000.0, 2000.])/1000
    right_lead_point = np.array([0.0, -9000.0, 2000.])/1000
    right_trail_point = np.array([4000.0, 9000.0, 2000.])/1000

    '''Project points'''
    wing_lead_left, wing_lead_left_coord = design.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_left ,wing_trail_left_coord = design.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_lead_right ,wing_lead_right_coord = design.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
    wing_trail_right ,wing_trail_right_coord = design.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])





