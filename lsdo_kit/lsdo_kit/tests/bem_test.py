import numpy as np
from lsdo_kit.problem import Problem
from lsdo_kit.design.design import Design
from lsdo_kit.simulation.cruise_simulation import Simulation, CruiseSimulation
from lsdo_kit.design.design_geometry.core.component import Component
from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes
from lsdo_kit.design.design_geometry.core.geometric_outputs import GeometricOuputs
from lsdo_kit.design.design_geometry.core.actuation import Actuation
from lsdo_kit.simulation.mesh.vlm_mesh import VLMMesh
from lsdo_kit.simulation.solver.vlm_solver import VlmSolver
from lsdo_kit.simulation.solver.bem_solver import BemSolver
from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes
from vedo import Points, Plotter, LegendBox
import matplotlib.pyplot as plt
import mplcursors

# Modify files that accept geometry to accept design instead and then unpack geometry

''' Creating File path for Rect Wing DesignGeometry'''
path_name = '../design/design_geometry/examples/CAD/'
file_name = 'rotor.stp'

''' Creating instances of problem and design class '''
problem = Problem()
design = Design(path_name + file_name)
geo = design.design_geometry

# for surface in geo.input_bspline_entity_dict.values():
#     vp_init = Plotter()
#     vps1 = Points(surface.control_points, r=8, c = 'blue')
#     # vps.append(vps2)
#     vp_init.show(vps1, f'{surface.name}', axes=1, viewup="z", interactive = True)

cruise = CruiseSimulation('cruise_sim')
cruise.nt = 2
cruise.add_solver_info('velocity', np.array([[50, 0, 1],[100,5,1]]))
cruise.add_solver_info('RPM', [2000, 3600])
cruise.add_solver_info('altitude', 100)

rotor_comp = Component(stp_entity_names=['Front_Rotor'], name='rotor')
rotor_comp.add_shape_parameter(property_name = 'scale_y', parameter_name='linear', order=1, num_cp=2, dv=False, val=np.ones(2))
design.add_component(rotor_comp)

up_direction = np.array([0., 0., 1.])
down_direction = np.array([0., 0., -1.])

root_lead = np.array([10.0792, -21.2592, 9.])
root_trail = np.array([10.6, -21.9329, 8.25])
tip_lead = np.array([10.0925, -18., 8.7])
tip_trail = np.array([10.33, -18, 8.])

''' Project points '''
root_lead_ps, root_lead_coord = design.project_points(root_lead, projection_direction = down_direction, projection_targets_names=["rotor"])
root_trail_ps , root_trail_coord = design.project_points(root_trail, projection_direction = down_direction, projection_targets_names=["rotor"])

tip_lead_ps , tip_lead_coord = design.project_points(tip_lead, projection_direction = down_direction, projection_targets_names=["rotor"])
tip_trail_ps , tip_trail_coord = design.project_points(tip_trail, projection_direction = up_direction, projection_targets_names=["rotor"])

root_length = geo.subtract_pointsets(root_lead_ps, root_trail_ps)
tip_length = geo.subtract_pointsets(tip_lead_ps, tip_trail_ps)

''' Create rotor camber surface '''
num_pts2 = [20]#20
num_pts1 = [5]#5
surface_shape1 = np.append(num_pts1,num_pts2)

root_chord = geo.perform_linear_interpolation(root_lead_ps, root_trail_ps, num_pts1)
tip_chord = geo.perform_linear_interpolation(tip_lead_ps, tip_trail_ps, num_pts1)

rotor_chord_surface = geo.perform_linear_interpolation(root_chord, tip_chord, surface_shape1)
geo.assemble(rotor_chord_surface)
chord_surface_points = geo.evaluate(rotor_chord_surface)

top_points, top_plot = geo.project_points(chord_surface_points, projection_targets_names=["Front_Rotor, 0, 2"])
bot_points, bot_plot = geo.project_points(chord_surface_points, projection_targets_names=["Front_Rotor, 0, 3"])

'''Create Camber Surface'''
surface_shape2 = np.append(surface_shape1,1)
output_parameters = 0.5*np.ones([num_pts1[0]*num_pts2[0],1])
camber_rotor = geo.perform_linear_interpolation(top_points, bot_points, surface_shape2, output_parameters = output_parameters)
camber_rotor.name = 'camber_surface_mesh'

geo.assemble()
geo.evaluate()

vps = []
vp_init = Plotter()
vps1 = Points(geo.input_bspline_entity_dict["Front_Rotor, 0, 2"].control_points, r=8, c = 'blue')
vps2 = Points(geo.input_bspline_entity_dict["Front_Rotor, 0, 3"].control_points, r=8, c = 'blue')
vps3 = Points(camber_rotor.physical_coordinates, r=8, c = 'red')
vps.append(vps1)
vps.append(vps2)
vps.append(vps3)
vp_init.show(vps, 'Projection', axes=1, viewup="z", interactive = True)


''' Calculate root length, '''
bem_solver = BemSolver(
    name = 'BEM',   
    airfoil = 'NACA_0012', 
    num_blades = 5, 
    num_radial = 25, 
    num_tangential = 25,
    )


cruise.add_solver(bem_solver)

''' Add Design to the Problem class and call the ProblemModel '''
geo.assemble()
geo.evaluate()
problem.set_design(design)
problem.add_simulation(cruise)

problem.assemble()
problem.run()
# problem.sim.check_partials()

# vps = []
# vp_init = Plotter()
# vps1 = Points(geo.total_cntrl_pts_vector, r=8, c = 'blue')
# vps2 = Points(root_lead_coord , r=15, c = 'red')
# vps3 = Points(root_trail_coord , r=15, c = 'red')
# vps4 = Points(tip_lead_coord , r=15, c = 'red')
# vps5 = Points(tip_trail_coord , r=15, c = 'red')

# vps.append(vps1)
# vps.append(vps2)
# vps.append(vps3)
# vps.append(vps4)
# vps.append(vps5)
# vp_init.addHoverLegend()
# vp_init.show(vps, 'Projection', axes=1, viewup="z", interactive = True)


# print('geo total control points shape', geo.total_cntrl_pts_vector.shape)

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# xs = geo.total_cntrl_pts_vector[:,0]
# ys = geo.total_cntrl_pts_vector[:,1]
# zs = geo.total_cntrl_pts_vector[:,2]
# lines = ax.scatter(xs, ys, zs)
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# mplcursors.cursor(hover=True)
# fig.set_size_inches(18.5, 18.5)
# plt.show()


# problem.set_design(design)
# problem.add_simulation(cruise)
# # problem.add_simulation(test_sim)

# problem.assemble()
# problem.run()





# velocity = np.array([[50, 0, 1], [50, 0, 1], [50, 0, 1]])
# test_sim = CruiseSimulation('test_sim')
# test_sim.add_solver_info('velocity', velocity)
# test_sim.nt = 3
# print('DesignGeometry Object in Design: ', design.design_geometry)

# ''' Creating a component and adding to design '''
# c = Component(stp_entity_names=['RectWing'], name='wing')
# c.add_shape_parameter(property_name = 'scale_y', parameter_name='linear', order=1, num_cp=2, dv=False, val=np.ones(2))
# design.add_component(c)

# print('Design Component Dictionary: ', design.component_dict)
# # print('DesignGeometry Component Dictionary: ', design.design_geometry.components)
# # print('DesignGeometry Component FFD Dictionary: ', design.design_geometry.components_ffd_dict)

# ''' Creating the chamber line mesh for mesh eval '''
# lead_edge_curve, trail_edge_curve, camber_surface_mesh = generate_meshes(geo)
# # print('camber surface mesh shape: ', camber_surface_mesh.shape)
# ''' Exporting mesh to txt file for Jiayao '''
# # geo.assemble()
# # geo.evaluate()
# # a_file = open("points.txt", "w")

# # for point in camber_surface_mesh.physical_coordinates:
# #     np.savetxt(a_file, point)

# # a_file.close()

# ''' Defining Points and directions for projections '''
# up_direction = np.array([0., 0., 1.])
# down_direction = np.array([0., 0., -1.])

# left_lead_point = np.array([0.0, -9000., 2000.])/1000
# left_trail_point = np.array([4000.0, -9000.0, 2000.])/1000
# right_lead_point = np.array([0.0, 9000.0, 2000.])/1000
# right_trail_point = np.array([4000.0, 9000.0, 2000.])/1000

# ''' Project points '''
# wing_lead_left, wing_lead_left_coord = design.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
# wing_trail_left ,wing_trail_left_coord = design.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

# wing_lead_right ,wing_lead_right_coord = design.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
# wing_trail_right ,wing_trail_right_coord = design.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

# ''' Creating a pointing vector across the wing for rotation '''
# output_parameters = np.array([0.25])
# quarter_left = design.perform_linear_interpolation(pointset_start=wing_lead_left, pointset_end=wing_trail_left, shape=(1,), output_parameters=output_parameters)
# quarter_right = design.perform_linear_interpolation(pointset_start=wing_lead_right, pointset_end=wing_trail_right, shape=(1,), output_parameters=output_parameters)

# ''' Creating an actuation object and adding it to the simulation '''
# # actuation_profile = np.linspace(0,1,100)
# actuation_profile = np.array([0])
# name = 'test_actuation'
# actuation_obj = Actuation(name='wing_rot', actuation_profile=actuation_profile, origin=quarter_left, pointset2=quarter_right, actuating_components=[c])
# cruise.add_actuations(design=design, actuation_list=[actuation_obj])

# # actuation_profile = np.linspace(0, np.pi/8, 3)
# # # actuation_profile = np.zeros((3,))
# # name = 'multiple_timesteps'
# # actuation_obj = Actuation(name='multi_wing_rot', actuation_profile=actuation_profile, origin=quarter_left, pointset2=quarter_right, actuating_components=[c])
# # test_sim.add_actuations(design=design, actuation_list=[actuation_obj])


# ''' Setting up Geometric Outputs and adding them to the design class '''
# # geo.assemble()
# geo_outputs = GeometricOuputs(geo=geo)
# geo_outputs.compute_displacement('test_displacement', wing_lead_left, wing_trail_left)
# geo_outputs.compute_magnitude('test_mag', wing_lead_left, wing_trail_left)
# # print('Geometric Outputs Dictionary: ', geo_outputs.geometric_outputs_dict)

# design.add_geometric_outputs(geo_outputs)
# # print('DesignGeometry Geometric Outputs: ', design.design_geometry.geometric_outputs)

# ''' Setting up the VLM mesh and creating inputs necessary for the VLM Solver '''
# vlm_mesh = VLMMesh('vlm_mesh', [camber_surface_mesh])
# vlm_solver = VlmSolver('vlm_solver', mesh=vlm_mesh)
# cruise.add_solver(vlm_solver)
# # test_sim.add_solver(vlm_solver)

# ''' Add Design to the Problem class and call the ProblemModel '''
# geo.assemble()
# geo.evaluate()

# problem.set_design(design)
# problem.add_simulation(cruise)
# # problem.add_simulation(test_sim)

# problem.assemble()
# problem.run()


# ''' Printing of VLM results '''

# for simulation_name, simulation_obj in problem.simulation_dict.items():
#     solver_name = list(simulation_obj.solver_dict.keys())[0]
#     for t in range(simulation_obj.nt):
        
#         print(f'{simulation_name} LIFT t={t}:', problem.sim[f'simulation_models.{simulation_name}.{solver_name}.vlm_model_{t}.L'])
#     print()

# for simulation_name, simulation_obj in problem.simulation_dict.items():
#     solver_name = list(simulation_obj.solver_dict.keys())[0]
#     for t in range(simulation_obj.nt):
        
#         print(f'{simulation_name} DRAG t={t}:', problem.sim[f'simulation_models.{simulation_name}.{solver_name}.vlm_model_{t}.D'])
#     print()

# for simulation_name, simulation_obj in problem.simulation_dict.items():
#     solver_name = list(simulation_obj.solver_dict.keys())[0]
#     for t in range(simulation_obj.nt):
        
#         print(f'{simulation_name} C_L t={t}:', problem.sim[f'simulation_models.{simulation_name}.{solver_name}.vlm_model_{t}.C_L'])
#     print()

# for simulation_name, simulation_obj in problem.simulation_dict.items():
#     solver_name = list(simulation_obj.solver_dict.keys())[0]
#     for t in range(simulation_obj.nt):
        
#         print(f'{simulation_name} C_D t={t}:', problem.sim[f'simulation_models.{simulation_name}.{solver_name}.vlm_model_{t}.C_D'])
#     print()
