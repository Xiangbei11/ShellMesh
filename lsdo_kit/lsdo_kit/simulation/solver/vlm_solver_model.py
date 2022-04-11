from VLM_package.VLM_system.vlm_system import VLMSystemModel
from VLM_package.VLM_outputs.compute_force.compute_outputs_group import Outputs
import numpy as np

from VLM_package.VLM_preprocessing.generate_simple_mesh import *

class VLMSolverModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('surface_names', types=list)
        self.parameters.declare('surface_shapes', types=list)
        self.parameters.declare('nt')
        # self.parameters.declare('free_stream_velocities', types=np.ndarray)

    def define(self):
        # add the mesh info
        surface_names = self.parameters['surface_names']
        surface_shapes = self.parameters['surface_shapes']
        nt = self.parameters['nt']

        for i in range(len(surface_names)):
           
            all_ts_mesh = self.declare_variable('all_timestep_' + surface_names[i], shape=(nt,) + surface_shapes[i])
            all_ts_vel  = self.declare_variable('all_timestep_frame_vel', shape=(nt,3))
            all_ts_vel  = -all_ts_vel

            for t in range(nt):
                self.register_output(f'mesh_{t}', all_ts_mesh[t,:,:,:])
                self.register_output(f'vel_{t}', all_ts_vel[t,:])
                
                model = Model()
                model.declare_variable('frame_vel', shape=(3,))
                model.declare_variable(surface_names[i], shape=(1,) + surface_shapes[i], val=np.ones((1,) + surface_shapes[i])) 

                model.add(
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
                model.add(sub, name='compute_lift_drag')
                self.add(model, f'vlm_model_{t}', promotes=[])

                self.connect(f'mesh_{t}', f'vlm_model_{t}.' + surface_names[i])
                self.connect(f'vel_{t}', f'vlm_model_{t}.frame_vel')

        # print('SURFACE NAMES: ', surface_names)
        # print('SURFACE SHAPES: ', surface_shapes)

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

        # '''
        # The mesh would also have to be connected from outside this model into this model. 
        # The shape of the mesh must be (1,nx,ny,3)
        # nx - chordwise number of points 
        # ny - spanwise number of points 
        # NOTE: We want to keep the number of chordwise points low. Spanwise number of points can go up to 50. 
        # '''

        # wing = self.declare_variable(surface_names[0], val=np.ones((1,) + surface_shapes[0])*10 ) 
        
        # self.add(
        #     VLMSystemModel(
        #         surface_names=surface_names,
        #         surface_shapes=surface_shapes,
        #         # frame_vel=frame_vel_val,
        #     ),
        #     'ODE_system')

        # eval_pts_names = [x + '_eval_pts_coords' for x in surface_names]
        # eval_pts_shapes = [(x[0] - 1, x[1] - 1, 3) for x in surface_shapes]

        # # compute lift and drag
        # sub = Outputs(
        #     surface_names=surface_names,
        #     surface_shapes=surface_shapes,
        #     eval_pts_names=eval_pts_names,
        #     eval_pts_shapes=eval_pts_shapes,
        # )
        # self.add(sub, name='compute_lift_drag')


        '''
        Outputs:
        Total Lift: L
        Total Drag: D
        Lift Coefficient: C_L
        Drag Coefficient: C_D
        
        '''

if __name__ == "__main__":

    nx = 20
    ny = 3
    surface_shapes = [(nx, ny, 3)]

    model = Model()
    mesh_val = generate_simple_mesh(nx, ny).reshape(1, nx, ny, 3)
    surface_names = ['wing']
    frame_vel_val = np.array([50, 0, 1])
    free_stream_velocities = -frame_vel_val

    wing = model.create_input('wing', val=mesh_val)

    submodel = VLMSolverModel(
        surface_names=surface_names,
        surface_shapes=surface_shapes,
        # free_stream_velocities=free_stream_velocities,
    )
    model.add(submodel, 'VLMSolverModel')
    sim = Simulator(model)

    sim.run()
    sim.prob.model.list_inputs(prom_name=True)
    # sim.prob.model.list_outputs(prom_name=True)
    print('lift', sim.prob['L'])
    print('drag', sim.prob['D'])
    # sim.visualize_implementation()