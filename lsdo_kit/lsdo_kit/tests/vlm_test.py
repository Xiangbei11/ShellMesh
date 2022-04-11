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
from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes
from vedo import Points, Plotter, LegendBox

import os
os.chdir("../lsdo_geo/lsdo_kit/tests")

# Modify files that accept geometry to accept design instead and then unpack geometry

''' Creating File path for Rect Wing DesignGeometry'''
path_name = '../design/design_geometry/examples/CAD/'
file_name = 'rect_wing.stp'

''' Creating instances of problem and design class '''
problem = Problem()
design = Design(path_name + file_name)
geo = design.design_geometry

# vps = []
# vp_init = Plotter()
# vps1 = Points(geo.total_cntrl_pts_vector, r=8, c = 'blue')
# vps.append(vps1)
# vp_init.show(vps, 'Projection', axes=1, viewup="z", interactive = True)

cruise = CruiseSimulation('cruise_sim')
cruise.add_solver_info('velocity', np.array([50, 0, 1]))

nt = 2
velocity = np.tile([50, 0, 1], (nt,1))
test_sim = CruiseSimulation('test_sim')
test_sim.add_solver_info('velocity', velocity)
test_sim.nt = nt
# print('DesignGeometry Object in Design: ', design.design_geometry)

''' Creating a component and adding to design '''
c = Component(stp_entity_names=['RectWing'], name='wing')
c.add_shape_parameter(property_name = 'scale_y', parameter_name='linear', order=1, num_cp=2, dv=False, val=np.ones(2))
design.add_component(c)

# print('Design Component Dictionary: ', design.component_dict)
# print('DesignGeometry Component Dictionary: ', design.design_geometry.components)
# print('DesignGeometry Component FFD Dictionary: ', design.design_geometry.components_ffd_dict)

''' Creating the chamber line mesh for mesh eval '''
lead_edge_curve, trail_edge_curve, camber_surface_mesh = generate_meshes(geo)

geo.assemble()
geo.evaluate()

# vps = []
# vp_init = Plotter()
# vps3 = Points(camber_surface_mesh.physical_coordinates, r=8, c = 'red')
# vps.append(vps3)
# vp_init.show(vps, 'Projection', axes=1, viewup="z", interactive = True)

# print('camber surface mesh shape: ', camber_surface_mesh.shape)
''' Exporting mesh to txt file for Jiayao '''
# geo.assemble()
# geo.evaluate()
# a_file = open("points.txt", "w")

# for point in camber_surface_mesh.physical_coordinates:
#     np.savetxt(a_file, point)

# a_file.close()

''' Defining Points and directions for projections '''
up_direction = np.array([0., 0., 1.])
down_direction = np.array([0., 0., -1.])

left_lead_point = np.array([0.0, -9000., 2000.])/1000
left_trail_point = np.array([4000.0, -9000.0, 2000.])/1000
right_lead_point = np.array([0.0, 9000.0, 2000.])/1000
right_trail_point = np.array([4000.0, 9000.0, 2000.])/1000

''' Project points '''
wing_lead_left, wing_lead_left_coord = design.project_points(left_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
wing_trail_left ,wing_trail_left_coord = design.project_points(left_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

wing_lead_right ,wing_lead_right_coord = design.project_points(right_lead_point, projection_direction = down_direction, projection_targets_names=["wing"])
wing_trail_right ,wing_trail_right_coord = design.project_points(right_trail_point, projection_direction = down_direction, projection_targets_names=["wing"])

''' Creating a pointing vector across the wing for rotation '''
output_parameters = np.array([0.25])
quarter_left = design.perform_linear_interpolation(pointset_start=wing_lead_left, pointset_end=wing_trail_left, shape=(1,), output_parameters=output_parameters)
quarter_right = design.perform_linear_interpolation(pointset_start=wing_lead_right, pointset_end=wing_trail_right, shape=(1,), output_parameters=output_parameters)

''' Creating an actuation object and adding it to the simulation '''
# actuation_profile = np.linspace(0,1,100)
actuation_profile = np.array([0])
name = 'test_actuation'
actuation_obj = Actuation(name='wing_rot', actuation_profile=actuation_profile, origin=quarter_left, pointset2=quarter_right, actuating_components=[c])
cruise.add_actuations(design=design, actuation_list=[actuation_obj])

actuation_profile = np.linspace(0, np.pi/8, nt)
# actuation_profile = np.zeros((3,))
name = 'multiple_timesteps'
actuation_obj = Actuation(name='multi_wing_rot', actuation_profile=actuation_profile, origin=quarter_left, pointset2=quarter_right, actuating_components=[c])
test_sim.add_actuations(design=design, actuation_list=[actuation_obj])


''' Setting up Geometric Outputs and adding them to the design class '''
# geo.assemble()
geo_outputs = GeometricOuputs(geo=geo)
geo_outputs.compute_displacement('test_displacement', wing_lead_left, wing_trail_left)
geo_outputs.compute_magnitude('test_mag', wing_lead_left, wing_trail_left)
# print('Geometric Outputs Dictionary: ', geo_outputs.geometric_outputs_dict)

design.add_geometric_outputs(geo_outputs)
# print('DesignGeometry Geometric Outputs: ', design.design_geometry.geometric_outputs)

''' Setting up the VLM mesh and creating inputs necessary for the VLM Solver '''
vlm_mesh = VLMMesh('vlm_mesh', [camber_surface_mesh])
vlm_solver = VlmSolver('vlm_solver', mesh=vlm_mesh)
cruise.add_solver(vlm_solver)
test_sim.add_solver(vlm_solver)

''' Add Design to the Problem class and call the ProblemModel '''
geo.assemble()
geo.evaluate()

problem.set_design(design)
# problem.add_simulation(cruise)
problem.add_simulation(test_sim)

problem.assemble()
problem.run()
problem.sim.prob.model.list_outputs(prom_name=True, shape=True)


''' Printing of VLM results '''

for simulation_name, simulation_obj in problem.simulation_dict.items():
    solver_name = list(simulation_obj.solver_dict.keys())[0]
    for t in range(simulation_obj.nt):
        
        print(f'{simulation_name} LIFT t={t}:', problem.sim[f'simulation_models.{simulation_name}.{solver_name}.vlm_model_{t}.L'])
    print()

for simulation_name, simulation_obj in problem.simulation_dict.items():
    solver_name = list(simulation_obj.solver_dict.keys())[0]
    for t in range(simulation_obj.nt):
        
        print(f'{simulation_name} DRAG t={t}:', problem.sim[f'simulation_models.{simulation_name}.{solver_name}.vlm_model_{t}.D'])
    print()

for simulation_name, simulation_obj in problem.simulation_dict.items():
    solver_name = list(simulation_obj.solver_dict.keys())[0]
    for t in range(simulation_obj.nt):
        
        print(f'{simulation_name} C_L t={t}:', problem.sim[f'simulation_models.{simulation_name}.{solver_name}.vlm_model_{t}.C_L'])
    print()

for simulation_name, simulation_obj in problem.simulation_dict.items():
    solver_name = list(simulation_obj.solver_dict.keys())[0]
    for t in range(simulation_obj.nt):
        
        print(f'{simulation_name} C_D t={t}:', problem.sim[f'simulation_models.{simulation_name}.{solver_name}.vlm_model_{t}.C_D_i'])
    print()
