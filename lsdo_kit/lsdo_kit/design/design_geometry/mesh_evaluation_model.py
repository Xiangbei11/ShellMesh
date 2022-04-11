import csdl
from csdl.operations.print_var import print_var
from csdl_om import Simulator
import numpy as np

class MeshEvaluationModel(csdl.Model):
    '''
    Maps from the OML control points to the mesh coordinates.
    '''

    def initialize(self):
        design = self.parameters.declare('design')
        meshes = self.parameters.declare('meshes')
        nt = self.parameters.declare('nt')
        # DesignGeometry assemble generates the absolute map that maps the geometry control points to the OML

    def define(self):
        # print('STARTING MESH EVALUATION')
        design = self.parameters['design']
        meshes = self.parameters['meshes']
        nt = self.parameters['nt']

        design_geometry_obj = design.design_geometry
        design_geometry_obj.assemble()
        eval_map = self.create_input('eval_map', val=design_geometry_obj.eval_map.todense())
        offset_map = self.create_input('offset_map', val=design_geometry_obj.offset_eval_map)

        control_points = self.declare_variable('control_points', shape=(nt,) + design.design_geometry.total_cntrl_pts_vector.shape)

        # print('Offset_map: ', design_geometry_obj.offset_eval_map.shape
        
        final_points = self.create_output('final_points', shape=(nt,) + (eval_map.shape[0], 3))
        # print('FINAL POINTS SHAPE: ', final_points.shape)
        # print('FINAL POINTS INDEXED SHAPE: ', final_points[0,:,0].shape)
        for t in range(nt):
            # print('TIME IN MESH EVAL: ', t)
            temp = csdl.reshape(control_points[t,:,:], new_shape=design.design_geometry.total_cntrl_pts_vector.shape)
            temp = csdl.matmat(eval_map, temp)
            temp = csdl.reshape(temp, new_shape=(1,) + temp.shape)
            final_points[t,:,:] = temp
        # final_points = csdl.matmat(eval_map, control_points)

        # print('Control Points Vector Shape: ', design_geometry_obj.total_cntrl_pts_vector.shape)
        # print('Eval Map Shape: ', design_geometry_obj.eval_map.shape)
        
        # self.register_output('points', final_points)

        # final_points = self.declare_variable('final_points', shape=(nt,) + design.design_geometry.total_cntrl_pts_vector.shape)

        # for t in range(nt):
        for mesh in meshes:
            # mesh_points = self.create_output(f'{mesh.name}', shape=(mesh.num_mesh_points,3))
            # mesh_starting_ind = 0
            # mesh_ending_ind = 0
            for pointset in mesh.pointset_list:
                test = self.create_output(f'{mesh.name}_{pointset.name}', shape=(nt,) + pointset.shape)
                pointset_starting_ind = pointset.output_starting_ind
                # print('TYPE: ', type(pointset_starting_ind))

                pointset_ending_ind = pointset_starting_ind + np.cumprod(pointset.shape)[-2]
                pointset_points = final_points[:, int(pointset_starting_ind):int(pointset_ending_ind), :]
                pointset_points_reshaped = csdl.reshape(pointset_points, new_shape=(nt,) + pointset.shape)
                # print('MESH SHAPE: ', pointset_points_reshaped.shape)
                test[:,:,:,:] = pointset_points_reshaped
                # self.register_output(f'{mesh.name}_{pointset.name}', pointset_points_reshaped)
                    
                    # mesh_ending_ind += (pointset_ending_ind - pointset_starting_ind)
                    # mesh_points[int(mesh_starting_ind):int(mesh_ending_ind), :] = pointset_points[t:]
                    # mesh_starting_ind = mesh_ending_ind

            # for feature in mesh.feature_list:
            #     for pointset in feature.pointset_list:
            #         pointset_starting_ind = pointset.output_starting_ind
            #         pointset_ending_ind = pointset_starting_ind + np.cumprod(pointset.shape)[-2]
            #         pointset_points = final_points[pointset_starting_ind:pointset_ending_ind, :]
            #         self.register_output(f'{mesh.name}_{feature.name}_geo', pointset_points)

            #         mesh_ending_ind += (pointset_ending_ind - pointset_starting_ind)
            #         mesh_points[mesh_starting_ind:mesh_ending_ind, :] = pointset_points
            #         mesh_starting_ind = mesh_ending_ind


if __name__ == "__main__":

    import numpy as np
    from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
    from lsdo_kit.old_files.mesh import Mesh

    import matplotlib.pyplot as plt

    from vedo import Points, Plotter, LegendBox

    from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
    import scipy.sparse as sps

    ''' Camber surface creation script '''
    path_name = '../examples/CAD/'
    file_name = 'rect_wing.stp'
    geo = DesignGeometry(path_name + file_name)

    print(geo.input_bspline_entity_dict.keys())
    ''' 
    Plotting the surfaces 
    '''

    for surface in geo.input_bspline_entity_dict.values():
        vp_init = Plotter()
        vps1 = Points(surface.control_points, r=8, c = 'blue')
        # vps.append(vps2)
        vp_init.show(vps1, f'{surface.name}', axes=1, viewup="z", interactive = True)


    up_direction = np.array([0., 0., 1.])
    down_direction = np.array([0., 0., -1.])

    # TODO Definitely want to automate this for all types of components
    #  could look into using normal vectors to distinguish between up and down
    wing_surface_names = [
        'RectWing, 0, 3', 'RectWing, 1, 8', 
        'RectWing, 0, 2', 'RectWing, 1, 9', 
        ]
    top_wing_surface_names = [
        'RectWing, 0, 3', 
        'RectWing, 1, 9', 
        ]
    bot_wing_surface_names = [
        'RectWing, 0, 2',
        'RectWing, 1, 8', 
        ]

    # TODO Make helper function for camber surface
    ''' Points to be projected'''
    lead_port_point = np.array([0., -9., 3.])
    lead_sb_point = np.array([0., 9., 3.])
    trail_sb_point = np.array([4., 9., 3.])
    trail_port_point = np.array([4., -9., 3.])

    '''Project points'''
    wing_lead_port, wing_lead_port_coord = geo.project_points(lead_port_point, projection_targets_names=top_wing_surface_names, projection_direction = down_direction)
    wing_lead_sb ,wing_lead_sb_coord = geo.project_points(lead_sb_point, projection_targets_names=top_wing_surface_names, projection_direction = down_direction)
    wing_trail_sb ,wing_trail_sb_coord = geo.project_points(trail_sb_point, projection_targets_names=top_wing_surface_names, projection_direction = down_direction)
    wing_trail_port ,wing_trail_port_coord = geo.project_points(trail_port_point, projection_targets_names=top_wing_surface_names,projection_direction = down_direction)

    '''Create chord surface using linear interpolation of lead/trail curves'''
    num_pts1 = [20]#20
    num_pts2 = [20]#5
    surface_shape1 = np.append(num_pts1,num_pts2)
    lead_edge_curve = geo.perform_linear_interpolation(wing_lead_port, wing_lead_sb, num_pts1)
    trail_edge_curve = geo.perform_linear_interpolation(wing_trail_port, wing_trail_sb, num_pts1)

    chord_surface = geo.perform_linear_interpolation(lead_edge_curve, trail_edge_curve, surface_shape1)

    '''Find absolute map of chord surface'''
    geo.assemble(chord_surface)
    chord_surface_mesh = geo.evaluate(chord_surface)

    '''Translate chord surface up and down'''
    chord_surface_trans_up = np.zeros([np.shape(chord_surface_mesh)[0],3])
    chord_surface_trans_up[:,:2] = chord_surface_mesh[:,:2]
    chord_surface_trans_up[:,2]= chord_surface_mesh[:,2]+8

    chord_surface_trans_down = np.zeros([np.shape(chord_surface_mesh)[0],3])
    chord_surface_trans_down[:,:2] = chord_surface_mesh[:,:2]
    chord_surface_trans_down[:,2]= chord_surface_mesh[:,2]-8

    '''Project from translated surfaces to wing surfaces'''
    top_points, top_plot = geo.project_points(chord_surface_trans_up, projection_targets_names=top_wing_surface_names, projection_direction=down_direction)
    bot_points, bot_plot = geo.project_points(chord_surface_trans_down, projection_targets_names=bot_wing_surface_names, projection_direction=up_direction)
    
    '''Create Camber Surface'''
    surface_shape2 = np.append(surface_shape1,1)
    output_parameters = 0.5*np.ones([400,1])
    camber_surface_mesh = geo.perform_linear_interpolation(top_points, bot_points, surface_shape2, output_parameters = output_parameters)

    vlm_mesh = Mesh('vlm_mesh')
    vlm_mesh.add_pointset(camber_surface_mesh, name="camber_surface")
    
    sim = Simulator(MeshEvaluationModel(geometry=geo, meshes=[vlm_mesh]))
    sim.run()
    sim.prob.model.list_inputs(prom_name=True)
    sim.prob.model.list_outputs(prom_name=True)
    
    # geo.assemble(camber_surface_mesh)

    '''Evaluate the physical coornadites of points to be fitted'''
    # camber_surface = geo.evaluate(camber_surface_mesh)
    camber_surface = sim.prob['vlm_mesh']
    geo.fit_bspline_entities(pointset_list=[camber_surface_mesh], output_vec=camber_surface)
    geo.write_iges('../examples/CAD/rect_wing_camber.igs')
    print(camber_surface)
    geo.evaluate()

    vp_init = Plotter()
    vps1 = Points(camber_surface, r=8, c = 'blue')
    # vps.append(vps2)
    vp_init.show(vps1, 'Camber', axes=1, viewup="z", interactive = True)

    '''Plot'''
    corner_projecting = []
    corner_projected = []
    vp_test = Plotter(N=1, axes=1)
    top_bot = []
    side = []
    cps = []
    surface_points = []
    chord_points = Points(chord_surface_mesh, r = 7, c = 'black').legend('Chord_Surface')
    #up_chord_points = Points(chord_surface_trans_up, r = 5, c = 'green').legend('Chord Translated Up')
    #down_chord_points = Points(chord_surface_trans_down, r = 5, c = 'green').legend('Chord Translated Down')
    top_surf_points = Points(top_points.physical_coordinates, r = 7, c = 'orange').legend('Top Surface Projection')
    camber_points = Points(camber_surface, r = 7, c = 'blue').legend('Camber_Surface')
    projected_surfaces = Points(np.vstack((top_plot, bot_plot)), r = 15, c = 'pink').legend('Top and Bottom')
    corner_projecting.append(Points(np.vstack((lead_port_point,lead_sb_point,trail_sb_point ,trail_port_point)),
        r=15, c='crimson').legend('Corner projecting'))
    corner_projected.append(Points(np.vstack((wing_lead_port_coord,wing_lead_sb_coord,wing_trail_sb_coord,wing_trail_port_coord)),
        r=15, c='darkgoldenrod').legend('Corner projected'))
    for target in wing_surface_names:
        cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))
        bspline_entity = geo.input_bspline_entity_dict[target]
        order_u = bspline_entity.order_u
        order_v = bspline_entity.order_v
        num_control_points_u = bspline_entity.shape[0]
        num_control_points_v = bspline_entity.shape[1]
        num_points_u = 50   # TODO might want to pass these in as input
        num_points_v = 50

        nnz = num_points_u * num_points_v * order_u * order_v
        data = np.zeros(nnz)
        row_indices = np.zeros(nnz, np.int32)
        col_indices = np.zeros(nnz, np.int32)

        knot_vector_u = bspline_entity.knots_u
        knot_vector_v = bspline_entity.knots_u
        u_vec = np.einsum('i,j->ij', np.linspace(0., 1., num_points_u), np.ones(num_points_v)).flatten()
        v_vec = np.einsum('i,j->ij', np.ones(num_points_u), np.linspace(0., 1., num_points_v)).flatten()

        get_basis_surface_matrix(
            order_u, num_control_points_u, 0, u_vec, knot_vector_u,
            order_v, num_control_points_v, 0, v_vec, knot_vector_v,
            num_points_u * num_points_v, data, row_indices, col_indices,
        )

        basis0 = sps.csc_matrix(
            (data, (row_indices, col_indices)),
            shape=(num_points_u * num_points_v, num_control_points_u * num_control_points_v),
        )
        pts = basis0.dot(bspline_entity.control_points) 
        surface_points.append(Points(pts, r=5, c='seagreen',alpha=0.3).legend('Surface points'))
        cps.append(Points(geo.input_bspline_entity_dict[target].control_points, r=8, c='cyan',alpha=0.3).legend('Control points'))        
    for surf in geo.output_bspline_entity_dict.values():
        if surf.name == 'Rib':
            bspline_fitted_cps = Points(surf.control_points, r=9, c='plum').legend('Fitted bspline')


    # vp_test.show(cps, surface_points, corner_projecting, corner_projected, 'Control points of surface to be projected', at=0, viewup="z", interactive = True)#, lb1
    vp_test.show(cps, surface_points, corner_projecting, corner_projected, chord_points, camber_points, 'Control points of surface to be projected', at=0, viewup="z", interactive = True)#, lb1

