from VLM_package.VLM_system.vlm_system import VLMSystemModel
from VLM_package.VLM_outputs.compute_force.compute_outputs_group import Outputs
import numpy as np

from VLM_package.VLM_preprocessing.generate_simple_mesh import *

# here nt is just a dummy variable that always equal to 2. since we are using a long wake panel,
# we can just make nt=2, delta_t=a large number.

# nx = 2
# ny = 50
# # ny = 100
# offset = 10

# frame_vel_val = np.array([1e-9, 0, -1])

# # multiple lifting surface
# # surface_names = ['wing', 'wing_1']
# # surface_shapes = [(nx, ny, 3), (nx, ny - 1, 3)]

# # single lifting surface
# surface_names = ['wing']
# surface_shapes = [(nx, ny, 3)]

# model_1 = csdl.Model()

# mesh_val = generate_simple_mesh(nx, ny).reshape(1, nx, ny, 3)



# frame_vel = model_1.create_input('frame_vel', val=frame_vel_val)


class VLMSolverModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('surface_names', types=list)
        self.parameters.declare('surface_shapes', types=list)
        # self.parameters.declare('free_stream_velocities', types=np.ndarray)

    def define(self):
        # add the mesh info
        surface_names = self.parameters['surface_names']
        surface_shapes = self.parameters['surface_shapes']

        print('SURFACE NAMES: ', surface_names)
        print('SURFACE SHAPES: ', surface_shapes)

        # free_stream_velocities = self.parameters['free_stream_velocities']
        # frame_vel_val = -free_stream_velocities
        # nx = 20
        # ny = 3
        # surface_shapes = [(nx, ny, 3)]

        # model = Model()
        # mesh_val = generate_simple_mesh(nx, ny).reshape(1, nx, ny, 3)
        # # surface_names = ['wing']
        # frame_vel_val = np.array([-50, 0, -1])
        # free_stream_velocities = -frame_vel_val

        # wing = self.create_input(surface_names[0], val=mesh_val)

        # frame_vel = self.declare_variable('frame_vel', shape=(3,))
        frame_vel = self.declare_variable('frame_vel', shape=(3,))
        '''
        The mesh would also have to be connected from outside this model into this model. 
        The shape of the mesh must be (1,nx,ny,3)
        nx - chordwise number of points 
        ny - spanwise number of points 
        NOTE: We want to keep the number of chordwise points low. Spanwise number of points can go up to 50. 
        '''


        wing = self.declare_variable('wing', shape=(1,) + surface_shapes[0])
        wing = wing + np.ones((1,) + surface_shapes[0])
        self.register_output('temp_wing', wing)

        # temp = self.create_input('temp', val=np.zeros((1,3,20,3)))
        # print('WING NAME: ', wing.name)
        # self.register_output('connect_wing', wing)
        # wing = wing + temp
        # The frame velocity will have to be connected to the declare variable in VLMSystemModel

        self.add(
            VLMSystemModel(
                surface_names=surface_names,
                surface_shapes=surface_shapes,
                # frame_vel=frame_vel_val,
            ),
            'ODE_system')

        eval_pts_names = [x + '_eval_pts_coords' for x in surface_names]
        eval_pts_shapes = [(x[0] - 1, x[1] - 1, 3) for x in surface_shapes]

        # compute lift and drag
        sub = Outputs(
            surface_names=surface_names,
            surface_shapes=surface_shapes,
            eval_pts_names=eval_pts_names,
            eval_pts_shapes=eval_pts_shapes,
        )
        self.add(sub, name='compute_lift_drag')


        '''
        Outputs:
        Total Lift: L
        Total Drag: D
        Lift Coefficient: C_L
        Drag Coefficient: C_D
        
        '''

if __name__ == "__main__":

    nx = 3
    ny = 20
    surface_shapes = [(nx, ny, 3)]

    f = np.loadtxt('../../points.txt')
    f = np.reshape(f, (1,nx,ny,3))

    model = Model()
    # mesh_val = generate_simple_mesh(nx, ny).reshape(1, nx, ny, 3)
    surface_names = ['wing']
    frame_vel_val = np.array([50, 0, 1])
    free_stream_velocities = -frame_vel_val
    print('TYPE MESH: ', type(f))
    model.create_input('vel', val=free_stream_velocities)
    # model.create_input('wing', val=f)
    model.create_input('wing', val=f)

    submodel = VLMSolverModel(
        surface_names=surface_names,
        surface_shapes=surface_shapes,
        # free_stream_velocities=free_stream_velocities,
    )
    model.add(submodel, 'VLMSolverModel', promotes=[])
    model.connect('wing', 'VLMSolverModel.wing')
    model.connect('vel', 'VLMSolverModel.frame_vel')
    sim = Simulator(model)

    sim.run()
    # sim.prob.model.list_inputs(prom_name=True)
    # sim.prob.model.list_outputs(prom_name=True)
    # print('WING NAME OUTSIDE: ', sim.prob['VLMSolverModel.connect_wing'])

    print('lift', sim.prob['VLMSolverModel.L'])
    print('drag', sim.prob['VLMSolverModel.D'])
    print('Cl', sim.prob['VLMSolverModel.C_L'])
    print('Cd', sim.prob['VLMSolverModel.C_D'])
    
    # sim.visualize_implementation()