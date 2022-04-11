import csdl
import numpy as np
from lsdo_kit.design.design_model import DesignModel
from lsdo_kit.simulation.top_simulation_model import TopSimulationModel
from lsdo_kit.simulation.solver.bem_solver import BemSolver


class ProblemModel(csdl.Model):
    def initialize(self):
        self.parameters.declare('design')
        self.parameters.declare('simulation_dict')

    def define(self):
        design = self.parameters['design']
        simulation_dict = self.parameters['simulation_dict']

        # design_model = DesignModel(design_obj=design)
        # self.add(design_model, name='design_model', promotes=[])

        simulation_model = TopSimulationModel(simulation_dict=simulation_dict, design=design)
        self.add(simulation_model, name='simulation_models')

        ''' CONNECTIONS BETWEEN DESIGN AND SIMULATION MODELS CAN ONLY BE DONE AT THE PROBLEM MODEL LEVEL'''
        # for simulation_name, simulation_obj in simulation_dict.items():
        #     if simulation_obj.actuation_dict:
        #         self.connect('design_model.design_geometry_model.ffd_model.application.total_geometry_control_points', f'simulation_models.{simulation_name}.actuation_model.total_geometry_control_points')

            # for solver_name, solver_obj in simulation_obj.solver_dict.items(): 
            #     if type(solver_obj) is BemSolver:
            #         for var_name, var_connection in solver_obj.design_solver_inputs_dict.items():
            #             self.connect(f'design_model.design_geometry_model.geometric_outputs.{var_name}', 
            #                     f'simulation_models.{simulation_name}.{solver_name}.{var_connection}')
        # for simulation_name, simulation_obj in simulation_dict.items():
        #     self.connect('design_model.design_geometry_model.ffd_model.application.total_geometry_control_points', f'simulation_models.{simulation_name}.actuation_model.total_geometry_control_points')



