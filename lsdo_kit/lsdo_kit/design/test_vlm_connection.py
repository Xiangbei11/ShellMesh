import csdl
from csdl_om import Simulator

import numbers 
import numpy as np
# from VLM_package.VLM_system.vlm_system import VLMSystemModel
from lsdo_kit.simulation.solver.vlm_solver_model import VLMSolverModel

from VLM_package.VLM_system.vlm_system import VLMSystemModel
from VLM_package.VLM_outputs.compute_force.compute_outputs_group import Outputs

from VLM_package.VLM_preprocessing.generate_simple_mesh import *
from lsdo_kit.design.design_geometry.tests.organized_test.mesh_creation_script import generate_meshes

from lsdo_kit.design.design import Design
from lsdo_kit.simulation.cruise_simulation import CruiseSimulation
from lsdo_kit.design.design_geometry.core.component import Component

class test(csdl.Model):
    def define(self):

        nx = 20
        ny = 3
        surface_shapes = [(nx, ny, 3)]

        mesh_val = generate_simple_mesh(nx, ny).reshape(1, nx, ny, 3)
        surface_names = ['wing']
        frame_vel_val = np.array([-50, 0, -1])
        free_stream_velocities = -frame_vel_val
        wing = self.create_input('wing', val=mesh_val)
        
        submodel = VLMSolverModel(
            surface_names=surface_names,
            surface_shapes=surface_shapes,
            free_stream_velocities=free_stream_velocities,
        )
        self.add(submodel, 'VLMSolverModel')


if __name__ == "__main__":

    path_name = 'design_geometry/examples/CAD/'
    file_name = 'rect_wing.stp'

    ''' Creating instances of problem and design class '''
    design = Design(path_name + file_name)
    geo = design.design_geometry
 
    cruise = CruiseSimulation('cruise_sim', velocity=50, density='pls_work')

    ''' Creating a component and adding to design '''
    c = Component(stp_entity_names=['RectWing'], name='wing')
    c.add_shape_parameter(property_name = 'scale_y', parameter_name='linear', order=1, num_cp=2, dv=False, val=np.ones(2))
    design.add_component(c)

    camber_surface_mesh = generate_meshes(geo)

    print('camber surface mesh shape: ', camber_surface_mesh.shape)


    nx = 20
    ny = 3
    surface_shapes = [(nx, ny, 3)]

    model = Model()
    mesh_val = generate_simple_mesh(nx, ny).reshape(1, nx, ny, 3)
    surface_names = ['wing']
    frame_vel_val = np.array([-50, 0, -1])
    free_stream_velocities = -frame_vel_val

    wing = model.create_input('wing', val=mesh_val)

    submodel = VLMSolverModel(
        surface_names=surface_names,
        surface_shapes=surface_shapes,
        free_stream_velocities=free_stream_velocities,
    )
    model.add(submodel, 'VLMSolverModel')
    sim = Simulator(model)
    sim.run()
