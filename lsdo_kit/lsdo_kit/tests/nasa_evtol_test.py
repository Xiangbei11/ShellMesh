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
from lsdo_kit.tests.thrust_vector_creation import generate_thrust_vector
from lsdo_kit.tests.center_mass_calculation import generate_CoM
from lsdo_kit.design.feature.point_mass import PointMass
from vedo import Points, Plotter, LegendBox
import matplotlib.pyplot as plt
import timeit

''' Creating File path for Rect Wing DesignGeometry'''
path_name = '../design/design_geometry/examples/CAD/'
file_name = 'fuse_fwing_rwing_tail.stp'

''' Creating instances of problem and design class '''
problem = Problem()

design = Design(path_name + file_name)
geo = design.design_geometry

nt = 1

velocity = np.tile(np.array([50, 0, 1]), (nt,1))
test_sim = CruiseSimulation('test_sim', velocity=velocity, density='pls_work') 
test_sim.nt = nt

''' Fuselage and Passengers and batteries'''
fuselage               = Component(stp_entity_names=['Fuselage'], name='fuselage')
design.add_component(fuselage)
fuselage_top, fuselage_bot, fuselage_CoM           = generate_CoM(geo, fuselage, np.array([12.77, 0.0, 5.197]))

passenger1_top, passenger1_bot, passenger1_CoM         = generate_CoM(geo, fuselage, np.array([7.014, -1.3, 4.907]))
passenger2_top, passenger2_bot, passenger2_CoM         = generate_CoM(geo, fuselage, np.array([7.014, 1.3, 4.907]))
passenger3_top, passenger3_bot, passenger3_CoM         = generate_CoM(geo, fuselage, np.array([9.848, 1.3, 4.907]))
passenger4_top, passenger4_bot, passenger4_CoM         = generate_CoM(geo, fuselage, np.array([9.848, -1.3, 4.907]))
passenger5_top, passenger5_bot, passenger5_CoM         = generate_CoM(geo, fuselage, np.array([12.681, 1.3, 4.907]))
passenger6_top, passenger6_bot, passenger6_CoM         = generate_CoM(geo, fuselage, np.array([12.681, -1.3, 4.907]))


battery1_top, battery1_bot, battery1_CoM               = generate_CoM(geo, fuselage, np.array([13.355, 0.0, 7.619]))
battery2_top, battery2_bot, battery2_CoM               = generate_CoM(geo, fuselage, np.array([15.72,  0.0, 3.923]))

''' Front and Rear Wing and Tail'''
fwing             = Component(stp_entity_names=['FrontWing'], name='front_wing')
design.add_component(fwing)
fwing_top, fwing_bot, fwing_CoM         = generate_CoM(geo, fwing, np.array([11.949, 0.0, 8.553]))

rwing             = Component(stp_entity_names=['RearWing'], name='Rear_wing')
design.add_component(rwing)
rwing_top, rwing_bot, rwing_CoM         = generate_CoM(geo, rwing, np.array([28.126, 0.0, 13.2]))

tail              = Component(stp_entity_names=['Tail'], name='tail')
design.add_component(tail)
tail_top, tail_bot, tail_CoM          = generate_CoM(geo, tail, np.array([26.11, 0.0, 9.927]))

''' Front wing nacelle '''
flnacelle1        = Component(stp_entity_names=['FrontLeftNacelle1'], name='flnacelle1')
design.add_component(flnacelle1)
flnacelle1_top, flnacelle1_bot, flnacelle1_CoM    = generate_CoM(geo, flnacelle1, np.array([8.953, 7.31, 8.144]))
flnacelle1_thrust = generate_thrust_vector(geo, flnacelle1)

flnacelle2          = Component(stp_entity_names=['FrontLeftNacelle2'], name='flnacelle2')
design.add_component(flnacelle2)
flnacelle2_top, flnacelle2_bot, flnacelle2_CoM      = generate_CoM(geo, flnacelle1, np.array([10.312, 14.624, 8.144]))
flnacelle2_thrust   = generate_thrust_vector(geo, flnacelle2)

flnacelle3          = Component(stp_entity_names=['FrontLeftNacelle3'], name='flnacelle3')
design.add_component(flnacelle3)
flnacelle3_top, flnacelle3_bot, flnacelle3_CoM      = generate_CoM(geo, flnacelle3, np.array([11.674, 21.859, 8.498]))
flnacelle3_thrust   = generate_thrust_vector(geo, flnacelle3)

frnacelle1          = Component(stp_entity_names=['FrontRightNacelle1'], name='frnacelle1')
design.add_component(frnacelle1)
frnacelle1_top, frnacelle1_bot, frnacelle1_CoM      = generate_CoM(geo, frnacelle1, np.array([8.953, -7.31, 8.144]))
frnacelle1_thrust   = generate_thrust_vector(geo, frnacelle1)

frnacelle2          = Component(stp_entity_names=['FrontRightNacelle2'], name='frnacelle2')
design.add_component(frnacelle2)
frnacelle2_top, frnacelle2_bot, frnacelle2_CoM      = generate_CoM(geo, frnacelle2, np.array([10.312, -14.624, 8.144]))
frnacelle2_thrust   = generate_thrust_vector(geo, frnacelle2)

frnacelle3          = Component(stp_entity_names=['FrontRightNacelle2'], name='frnacelle2')
design.add_component(frnacelle3)
frnacelle3_top, frnacelle3_bot, frnacelle3_CoM      = generate_CoM(geo, frnacelle3, np.array([11.674, -21.859, 8.498]))
frnacelle3_thrust   = generate_thrust_vector(geo, frnacelle3)

''' Front wing rotors '''

flrotor1_top, flrotor1_bot, flrotor1_CoM      = generate_CoM(geo, flnacelle1, np.array([7.548, 7.41, 8.139]))
flrotor2_top, flrotor2_bot, flrotor2_CoM      = generate_CoM(geo, flnacelle2, np.array([8.908, 14.624, 8.139]))
flrotor3_top, flrotor3_bot, flrotor3_CoM      = generate_CoM(geo, flnacelle3, np.array([10.27, 21.859, 8.5]))
frrotor1_top, frrotor1_bot, frrotor1_CoM      = generate_CoM(geo, frnacelle1, np.array([7.548, -7.41, 8.139]))
frrotor2_top, frrotor2_bot, frrotor2_CoM      = generate_CoM(geo, frnacelle2, np.array([8.908, -14.624, 8.139]))
frrotor3_top, frrotor3_bot, frrotor3_CoM      = generate_CoM(geo, frnacelle3, np.array([10.27, -21.859, 8.5]))


''' Front wing motor '''

flmotor1_top, flmotor1_bot, flmotor1_CoM        = generate_CoM(geo, flnacelle1, np.array([8.303, 7.41, 8.139]))
flmotor2_top, flmotor2_bot, flmotor2_CoM        = generate_CoM(geo, flnacelle2, np.array([9.662, 14.624, 8.139]))
flmotor3_top, flmotor3_bot, flmotor3_CoM        = generate_CoM(geo, flnacelle3, np.array([11.025, 21.859, 8.5]))
frmotor1_top, frmotor1_bot, frmotor1_CoM        = generate_CoM(geo, frnacelle1, np.array([8.303, -7.41, 8.139]))
frmotor2_top, frmotor2_bot, frmotor2_CoM        = generate_CoM(geo, frnacelle2, np.array([9.662, -14.624, 8.139]))
frmotor3_top, frmotor3_bot, frmotor3_CoM        = generate_CoM(geo, frnacelle3, np.array([11.025, -21.859, 8.5]))

''' Rear wing nacelle '''
rlnacelle1          = Component(stp_entity_names=['RearLeftNacelle1'], name='rlnacelle1')
design.add_component(rlnacelle1)
rlnacelle1_top, rlnacelle1_bot, rlnacelle1_CoM      = generate_CoM(geo, rlnacelle1, np.array([25.823, 5.04, 13.143]))
rlnacelle1_thrust   = generate_thrust_vector(geo, frnacelle1)

rrnacelle1          = Component(stp_entity_names=['RearRightNacelle1'], name='rrnacelle1')
design.add_component(rrnacelle1)
rrnacelle1_top, rrnacelle1_bot, rrnacelle1_CoM      = generate_CoM(geo, rrnacelle1, np.array([25.823, -5.04, 13.143]))
rrnacelle1_thrust   = generate_thrust_vector(geo, rrnacelle1)

''' Rear wing motor '''

rlmotor1_top, rlmotor1_bot, rlmotor1_CoM      = generate_CoM(geo, rlnacelle1, np.array([25.181, 5.04, 13.143]))
rrmotor1_top, rrmotor1_bot, rrmotor1_CoM      = generate_CoM(geo, rrnacelle1, np.array([25.181, -5.04, 13.143]))

''' Rear wing rotors '''
rlrotor1_top, rlrotor1_bot, rlrotor1_CoM      = generate_CoM(geo, rlnacelle1, np.array([24.426, 5.04, 13.143]))
rrrotor1_top, rrrotor1_bot, rrrotor1_CoM      = generate_CoM(geo, rrnacelle1, np.array([24.426, -5.04, 13.143]))

''' Point Mass Features ALL MASSES IN KG'''
pointmasses = PointMass('pointmasses')
pointmasses.add_pointmass(fuselage_CoM, 400)

pointmasses.add_pointmass(passenger1_CoM, 100)
pointmasses.add_pointmass(passenger2_CoM, 100)
pointmasses.add_pointmass(passenger3_CoM, 100)
pointmasses.add_pointmass(passenger4_CoM, 100)
pointmasses.add_pointmass(passenger5_CoM, 100)
pointmasses.add_pointmass(passenger6_CoM, 100)

pointmasses.add_pointmass(battery1_CoM, 20)
pointmasses.add_pointmass(battery2_CoM, 20)

pointmasses.add_pointmass(fwing_CoM, 150)
pointmasses.add_pointmass(rwing_CoM, 100)
pointmasses.add_pointmass(tail_CoM,  50)

pointmasses.add_pointmass(flnacelle1_CoM, 20)
pointmasses.add_pointmass(flnacelle2_CoM, 20)
pointmasses.add_pointmass(flnacelle3_CoM, 20)
pointmasses.add_pointmass(frnacelle1_CoM, 20)
pointmasses.add_pointmass(frnacelle2_CoM, 20)
pointmasses.add_pointmass(frnacelle3_CoM, 20)

pointmasses.add_pointmass(rlnacelle1_CoM, 20)
pointmasses.add_pointmass(rrnacelle1_CoM, 20)

pointmasses.add_pointmass(flrotor1_CoM, 10)
pointmasses.add_pointmass(flrotor2_CoM, 10)
pointmasses.add_pointmass(flrotor3_CoM, 10)
pointmasses.add_pointmass(frrotor1_CoM, 10)
pointmasses.add_pointmass(frrotor2_CoM, 10)
pointmasses.add_pointmass(frrotor3_CoM, 10)

pointmasses.add_pointmass(rlrotor1_CoM, 10)
pointmasses.add_pointmass(rrrotor1_CoM, 10)

pointmasses.add_pointmass(flmotor1_CoM, 10)
pointmasses.add_pointmass(flmotor2_CoM, 10)
pointmasses.add_pointmass(flmotor3_CoM, 10)
pointmasses.add_pointmass(frmotor1_CoM, 10)
pointmasses.add_pointmass(frmotor2_CoM, 10)
pointmasses.add_pointmass(frmotor3_CoM, 10)

pointmasses.add_pointmass(rlmotor1_CoM, 10)
pointmasses.add_pointmass(rrmotor1_CoM, 10)

design.add_feature(pointmasses)


# print('MAX NACELLE POS: ', flnacelle1.x_min)
# print('NACELLE EMBEDDED POINTS: ', flnacelle1.embedded_entities_control_points)
# print('MAX ROTOR POS: ', flrotor1.x_max)
# print('MIN ROTOR POS: ', flrotor1.x_min)

# for name in flnacelle1.embedded_entity_names:
#     surface = geo.input_bspline_entity_dict[name]
#     vp_init = Plotter()
#     vps1 = Points(surface.control_points, r=8, c = 'blue')
#     # vps.append(vps2)
#     vp_init.show(vps1, f'{surface.name}', axes=1, viewup="z", interactive = True)

# for name in flrotor1.embedded_entity_names:
#     surface = geo.input_bspline_entity_dict[name]
#     vp_init = Plotter()
#     vps1 = Points(surface.control_points, r=8, c = 'blue')
#     # vps.append(vps2)
#     vp_init.show(vps1, f'{surface.name}', axes=1, viewup="z", interactive = True)


''' Defining Points and directions for projections '''
up_direction = np.array([0., 0., 1.])
down_direction = np.array([0., 0., -1.])
mm2ft = 304.8

''' Front wing corner points '''
front_left_lead_point = np.array([3962.33, 6662.76, 2591])/mm2ft
front_left_trail_point = np.array([4673.88, 6662.76, 2590.8])/mm2ft
front_right_lead_point = np.array([3962.33, -6662.76, 2591])/mm2ft
front_right_trail_point = np.array([4673.88, -6662.76, 2590.8])/mm2ft

front_fuselage = np.array([3600.276, 453.0, 2600])/mm2ft
back_fuselage = np.array([3600.276, -453.0, 2600])/mm2ft


''' Project points '''
front_left_lead_point, front_left_lead_point_coord = design.project_points(front_left_lead_point, projection_direction = down_direction, projection_targets_names=["front_wing"])
front_left_trail_point ,front_left_trail_point_coord = design.project_points(front_left_trail_point, projection_direction = down_direction, projection_targets_names=["front_wing"])

front_right_lead_point ,front_right_lead_point_coord = design.project_points(front_right_lead_point, projection_direction = down_direction, projection_targets_names=["front_wing"])
front_right_trail_point ,front_right_trail_point_coord = design.project_points(front_right_trail_point, projection_direction = down_direction, projection_targets_names=["front_wing"])

front_fuse, front_fuse_coord = design.project_points(front_fuselage, projection_direction = down_direction, projection_targets_names=["front_wing"])
back_fuse, back_fuse_coord = design.project_points(back_fuselage, projection_direction = down_direction, projection_targets_names=["front_wing"])

''' Creating a pointing vector across the wing for rotation '''
output_parameters = np.array([0.75])
quarter_left = design.perform_linear_interpolation(pointset_start=front_left_lead_point, pointset_end=front_left_trail_point, shape=(1,), output_parameters=output_parameters)
quarter_right = design.perform_linear_interpolation(pointset_start=front_right_lead_point, pointset_end=front_right_trail_point, shape=(1,), output_parameters=output_parameters)

''' Define the actuation profile for the FRONT wings and rotors '''
actuation_profile = np.linspace(0, -np.pi/2, nt)
# actuation_profile = np.zeros((3,))
name = 'multiple_timesteps'
actuating_comps = [fwing, flnacelle1, flnacelle2, flnacelle3, frnacelle1, frnacelle2, frnacelle3]
actuation_obj = Actuation(name='multi_wing_rot', actuation_profile=actuation_profile, origin=front_fuse, pointset2=back_fuse, actuating_components=actuating_comps)
test_sim.add_actuations(design=design, actuation_list=[actuation_obj])
# test_sim.add_actuations(design=design, actuation_list=[actuation_obj1])


''' Add Design to the Problem class and call the ProblemModel '''

# geo.assemble()
geo.evaluate()

problem.set_design(design)
problem.add_simulation(test_sim)

problem.assemble()
problem.run()


# vp = Plotter()
# for t in range(180):
#     vps = Points(problem.sim['simulation_models.test_sim.actuation_model.actuated_control_points'][t,:,:], r=8, c = 'red')
#     vp.show(vps, f'Actuation_{t}', axes=1, viewup="z", interactive=False)
#     vp.screenshot(f'actuation{t}')



# ''' Rear wing corner points'''
# rear_left_lead_point = np.array([8277.012, 1536.289, 4005.94])/mm2ft
# rear_left_trail_point = np.array([9123.56, 1536.289, 4005.94])/mm2ft
# rear_right_lead_point = np.array([8277.012, -1536.289, 4005.94])/mm2ft
# rear_right_trail_point = np.array([9123.56, -1536.289, 4005.94])/mm2ft



# fig = plt.figure()
# ax = plt.axes(projection ='3d')
# x = geo.total_cntrl_pts_vector[:,0]
# y = geo.total_cntrl_pts_vector[:,1]
# z = geo.total_cntrl_pts_vector[:,2]
# ax.scatter(x, y, z, 'green')
# # ax.set_title(f'{surface.name}')
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# plt.show()


# print(geo.total_cntrl_pts_vector.shape)
# vp_init = Plotter()
# vps = []
# vps1 = Points(geo.total_cntrl_pts_vector, r=8, c = 'blue')
# # vps.append(vps2)
# vp_init.show(vps1, 'Total_control_points', axes=1, viewup="z", interactive = True)

# ''' Project points onto front wing '''
# front_wing_lead_left, front_wing_lead_left_coord = design.project_points(front_left_lead_point, projection_direction = down_direction, projection_targets_names=["front_wing"])
# front_wing_trail_left ,front_wing_trail_left_coord = design.project_points(front_left_trail_point, projection_direction = down_direction, projection_targets_names=["front_wing"])
# front_wing_lead_right ,front_wing_lead_right_coord = design.project_points(front_right_lead_point, projection_direction = down_direction, projection_targets_names=["front_wing"])
# front_wing_trail_right ,front_wing_trail_right_coord = design.project_points(front_right_trail_point, projection_direction = down_direction, projection_targets_names=["front_wing"])

# print(geo.total_cntrl_pts_vector.shape)

# for surface in geo.input_bspline_entity_dict.values():
# vp_init = Plotter()
# vps = []
# vps1 = Points(geo.total_cntrl_pts_vector[-1:3], r=8, c = 'blue')
# # vps.append(vps2)
# vp_init.show(vps1, 'Total_control_points', axes=1, viewup="z", interactive = True)




# for surface in geo.input_bspline_entity_dict.values():
#     fig = plt.figure()
#     ax = plt.axes(projection ='3d')
#     x = surface.control_points[:,0]
#     y = surface.control_points[:,1]
#     z = surface.control_points[:,2]
#     ax.scatter(x, y, z, 'green')
#     ax.set_title(f'{surface.name}')
#     ax.set_xlabel('X Label')
#     ax.set_ylabel('Y Label')
#     ax.set_zlabel('Z Label')
#     plt.show()


# for surface in geo.input_bspline_entity_dict.values():
#     vp_init = Plotter()
#     vps1 = Points(surface.control_points, r=8, c = 'blue')
#     # vps.append(vps2)
#     vp_init.show(vps1, f'{surface.name}', axes=1, viewup="z", interactive = True)


# print(geo.total_cntrl_pts_vector.shape)